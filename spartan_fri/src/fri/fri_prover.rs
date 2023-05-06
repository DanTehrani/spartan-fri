use super::fft::fft;
use super::tree::MerkleTree;
use super::unipoly::UniPoly;
use super::utils::sample_indices;
use super::{FriProof, LayerProof};
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

    fn commit(
        &self,
        codeword: &[F],
        transcript: &mut Transcript<F>,
    ) -> (Vec<Vec<F>>, Vec<MerkleTree<F>>) {
        let mut domain = self.domain.clone();

        let mut codewords = vec![codeword.to_vec()];
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

    fn query(
        &self,
        codewords: &[Vec<F>],
        trees: &[MerkleTree<F>],
        indices: &[usize],
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

            let mut openings = vec![];
            for j in 0..self.num_colinearity_checks {
                let a_y = codeword[a_indices[j]];
                let a_y_proof = trees[i].open(a_y);
                let b_y = codeword[b_indices[j]];
                let b_y_proof = trees[i].open(b_y);

                let c_y = codeword[c_indices[j]];
                let c_y_proof = trees[i].open(c_y);

                openings.push((a_y_proof, b_y_proof, c_y_proof));
            }

            queries.push(LayerProof { openings })
        }

        queries
    }

    pub fn prove_degree(&self, poly: &UniPoly<F>, transcript: &mut Transcript<F>) -> FriProof<F> {
        assert!(poly.degree().is_power_of_two());

        let mut coeffs_expanded: Vec<F> = poly.coeffs.clone();
        coeffs_expanded.resize(self.domain.len(), F::zero());

        let codewords = fft(&coeffs_expanded, &self.domain);

        let (codewords, trees) = self.commit(&codewords, transcript);

        let indices = sample_indices(
            self.num_colinearity_checks,
            codewords[0].len(),                   // Length of the initial codeword
            codewords[codewords.len() - 2].len(), // Length of the reduced codeword
            transcript,
        );

        let queries = self.query(&codewords, &trees, &indices);

        FriProof {
            reduced_codeword: codewords[codewords.len() - 1].clone(),
            queries,
        }
    }
}
