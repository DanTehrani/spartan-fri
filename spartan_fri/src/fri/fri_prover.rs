use super::fft::fft;
use super::tree::MerkleTree;
use super::unipoly::UniPoly;
use super::utils::sample_indices;
use super::{FriEvalProof, FriProof, LayerProof};
use ff::PrimeField;
use pasta_curves::arithmetic::FieldExt;

use crate::transcript::Transcript;

pub struct FriProver<F: PrimeField> {
    domain: Vec<F>,
    // Number of colinearity checks per round
    num_colinearity_checks: usize,
}

impl<F> FriProver<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub fn new(max_degree: usize) -> Self {
        // TODO: Allow arbitrary degree
        assert!(max_degree.is_power_of_two());

        // Are these params OK?
        let expansion_factor = 2;
        let num_colinearity_checks = 2;

        let root_of_unity = F::root_of_unity();
        let root_of_unity_log_2 = F::S;

        let domain_order = (max_degree * expansion_factor).next_power_of_two();

        // Generator for the subgroup with order _subgroup_order_ in the field
        let domain_generator = root_of_unity.pow(&[
            2u32.pow(32 - ((domain_order as f64).log2() as u32)) as u64,
            0,
            0,
            0,
        ]);

        let domain = (0..domain_order)
            .map(|i| domain_generator.pow(&[i as u64, 0, 0, 0]))
            .collect();

        // Compute the domain generator from the root of unity
        Self {
            domain,
            num_colinearity_checks,
        }
    }

    fn num_rounds(&self) -> usize {
        let domain_order = self.domain.len();
        ((domain_order as f64).log2() as usize) - 3 // this `3` is just random
    }

    fn fold(&self, codeword: &[F], domain: &[F], alpha: F) -> Vec<F> {
        assert!(codeword.len() == domain.len());
        let two_inv = F::from(2).invert().unwrap();
        let one = F::from(1);

        let n = domain.len();

        let mut folded_codeword = vec![];
        for i in 0..(n / 2) {
            // f*(w^2i) = 1/2 * ((1 + alpha * w^-i) * f(w^i) + (1 - alpha * w^-i) * f(-w^i))
            // w^(n/2) = -1
            // -w^i = domain[i + n/2]

            //  let omega_pow_minus_i = self.omega.pow(&[i as u64, 0, 0, 0]).invert().unwrap();
            let omega_pow_minus_i = domain[n - 1 - i];

            let f_star_eval = two_inv
                * ((one + alpha * omega_pow_minus_i) * codeword[i]
                    + (one - alpha * omega_pow_minus_i) * codeword[i + (n / 2)]);
            folded_codeword.push(f_star_eval);
        }

        folded_codeword
    }

    fn eval_poly_over_domain(&self, poly: &UniPoly<F>) -> Vec<F> {
        let mut coeffs_expanded: Vec<F> = poly.coeffs.clone();
        coeffs_expanded.resize(self.domain.len(), F::zero());

        let evals = fft(&coeffs_expanded, &self.domain);
        evals
    }

    fn commit(
        &self,
        evals_d: Vec<F>,
        transcript: &mut Transcript<F>,
    ) -> (Vec<Vec<F>>, Vec<MerkleTree<F>>) {
        let mut codewords = vec![evals_d];

        let mut domain = self.domain.clone();
        let mut trees = vec![];

        for i in 0..self.num_rounds() {
            let current_codeword = &codewords[i];

            let mut tree = MerkleTree::new();
            let root = tree.commit(current_codeword);

            transcript.append_fe(&root);
            trees.push(tree);

            let alpha = transcript.challenge_fe();

            let next_codeword = self.fold(current_codeword, &domain, alpha);
            let mut domain_unique = vec![];
            domain.iter().map(|x| x.square()).for_each(|x| {
                if !domain_unique.contains(&x) {
                    domain_unique.push(x);
                }
            });
            domain = domain_unique;

            codewords.push(next_codeword.to_vec())
        }

        (codewords, trees)
    }

    fn query_with_f(
        &self,
        codewords: &[Vec<F>],
        trees: &[MerkleTree<F>],
        indices: &[usize],
        f_codewords: &[F],
        f_tree: &MerkleTree<F>,
        z: F,
        y: F,
    ) -> Vec<LayerProof<F>> {
        // A domain: w^i
        // B domain: w^{n/2 + i}
        // C domain: w^{2i}

        assert!(indices.len() == self.num_colinearity_checks);
        let mut indices = indices.to_vec();

        let mut queries = vec![];

        for (i, codeword) in codewords.iter().enumerate() {
            // Skip the last codeword since the verifier's gonna check it directly
            if i == codewords.len() - 1 {
                continue;
            }

            // Halve the range of the indices
            indices = indices
                .iter()
                .map(|index| {
                    if codeword.len() == 1 {
                        0
                    } else {
                        index % (codeword.len() / 2)
                    }
                })
                .collect::<Vec<usize>>();

            let a_indices = indices.clone();
            let b_indices = indices
                .iter()
                .map(|index| (codeword.len() / 2) + index)
                .collect::<Vec<usize>>();
            let c_indices = indices.clone();

            let mut q_openings = vec![];
            let mut f_openings = vec![];
            for j in 0..self.num_colinearity_checks {
                let a_y = codeword[a_indices[j]];
                let q_a_y_proof = trees[i].open(a_y);
                let f_a_y_proof = f_tree.open(f_codewords[a_indices[j]]);

                println!(
                    "p q(x) ?=  {:?} {:?}",
                    (f_codewords[a_indices[j]] - y)
                        * (self.domain[a_indices[j]] - z).invert().unwrap(),
                    a_y
                );

                let b_y = codeword[b_indices[j]];
                let q_b_y_proof = trees[i].open(b_y);
                let f_b_y_proof = f_tree.open(f_codewords[b_indices[j]]);

                /*
                println!(
                    "p: {:?} {:?}",
                    //                    (f_b_y_proof.leaf - y) * (self.domain[b_indices[j]] - z).invert().unwrap()
                    f_b_y_proof.leaf,
                    self.domain[b_indices[j]]
                );
                 */

                let c_y = codeword[c_indices[j]];
                let q_c_y_proof = trees[i].open(c_y);
                let f_c_y_proof = f_tree.open(f_codewords[c_indices[j]]);

                q_openings.push((q_a_y_proof, q_b_y_proof, q_c_y_proof));
                f_openings.push((f_a_y_proof, f_b_y_proof, f_c_y_proof));
            }

            queries.push(LayerProof {
                q_openings,
                f_openings,
            })
        }

        queries
    }

    /*
    pub fn prove_degree(&self, poly: &UniPoly<F>, transcript: &mut Transcript<F>) -> FriProof<F> {
        assert!(poly.degree().is_power_of_two());

        let evals_d = self.eval_poly_over_domain(poly);

        let (codewords, trees) = self.commit(evals_d, transcript);

        let indices = sample_indices(
            self.num_colinearity_checks,
            codewords[0].len(),                   // Length of the initial codeword
            codewords[codewords.len() - 2].len(), // Length of the reduced codeword
            transcript,
        );

        let queries = self.query_with_f(&codewords, &trees, &indices,);

        FriProof {
            reduced_codeword: codewords[codewords.len() - 1].clone(),
            queries,
        }
    }
     */

    pub fn prove_eval(
        &self,
        f: &UniPoly<F>,
        z: F,
        transcript: &mut Transcript<F>,
    ) -> FriEvalProof<F> {
        let y = f.eval(z);

        // Commit to `f`
        let f_evals_d = self.eval_poly_over_domain(f);
        let mut f_tree = MerkleTree::new();
        let f_comm = f_tree.commit(&f_evals_d);

        let q_evals_d = f_evals_d
            .iter()
            .zip(&self.domain)
            .map(|(e, d)| (*e - y) * (*d - z).invert().unwrap())
            .collect();

        // Prove proximity of q_evals_d to a low-degree polynomial

        let (q_codewords, q_trees) = self.commit(q_evals_d, transcript);

        let indices = sample_indices(
            self.num_colinearity_checks,
            q_codewords[0].len(), // Length of the initial codeword
            q_codewords[q_codewords.len() - 2].len(), // Length of the reduced codeword
            transcript,
        );

        let queries =
            self.query_with_f(&q_codewords, &q_trees, &indices, &f_evals_d, &f_tree, z, y);
        // Compute q_evals from the f_evals!
        // q(x) = (f(x) - y) (x  - z)
        // f(x) = z
        //        let q_evals =
        // Compute q(x)

        let q_ld_proof = FriProof {
            reduced_codeword: q_codewords[q_codewords.len() - 1].clone(),
            queries,
        };

        FriEvalProof {
            z,
            y,
            q_ld_proof,
            f_comm,
        }
    }
}
