#![allow(non_snake_case)]

use super::tree::{MerkleProof, MerkleTree};
use super::unipoly::UniPoly;
use super::utils::{reduce_indices, sample_indices};
use super::{BatchedPolyEvalProof, FRIConfig};
use crate::transcript::Transcript;
use ark_std::{end_timer, start_timer};
use pasta_curves::arithmetic::FieldExt;

pub struct FRIPolyCommitProver<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub config: FRIConfig<F>,
}

impl<F: FieldExt> FRIPolyCommitProver<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub fn new(config: FRIConfig<F>) -> Self {
        Self { config }
    }

    fn fold(&self, codeword: &[F], domain: &[F], alpha: F) -> Vec<F> {
        assert!(codeword.len() == domain.len());
        let two_inv = F::from(2).invert().unwrap();
        let one = F::from(1);

        let n = domain.len();
        let mut folded_codeword = Vec::with_capacity(n / 2);
        for i in 0..(n / 2) {
            // f*(w^2i) = 1/2 * ((1 + alpha * w^-i) * f(w^i) + (1 - alpha * w^-i) * f(-w^i))
            // w^(n/2) = -1
            // -w^i = domain[i + n/2]

            let omega_pow_minus_i = domain[i].invert().unwrap();

            let f_star_eval = two_inv
                * ((one + alpha * omega_pow_minus_i) * codeword[i]
                    + (one - alpha * omega_pow_minus_i) * codeword[i + (n / 2)]);
            folded_codeword.push(f_star_eval);
        }

        folded_codeword
    }

    pub fn commit(
        &self,
        codeword: Vec<F>,
        transcript: &mut Transcript<F>,
    ) -> (Vec<Vec<F>>, Vec<MerkleTree<F>>) {
        let num_rounds = self.config.num_rounds();
        let mut oracles = Vec::<MerkleTree<F>>::with_capacity(num_rounds);
        let mut codewords = Vec::<Vec<F>>::with_capacity(num_rounds);

        let mut codeword = codeword;
        let mut challenges = Vec::with_capacity(num_rounds);
        for i in 0..self.config.num_rounds() {
            if i != 0 {
                // For i = 0, the verifier derives the opening from the oracle
                // of f(X), so we don't need to commit to it.
                let mut tree = MerkleTree::new();
                let comm = tree.commit(&codeword);
                oracles.push(tree);

                transcript.append_fe(&comm);
            }

            // "Fold" the codeword
            let x = transcript.challenge_fe();
            challenges.push(x);
            codeword = self.fold(&codeword, &self.config.L[i], x);
            codewords.push(codeword.clone());
        }

        (codewords, oracles)
    }

    pub fn query(
        &self,
        codewords: &[Vec<F>],
        trees: &[MerkleTree<F>],
        transcript: &mut Transcript<F>,
    ) -> Vec<Vec<(MerkleProof<F>, MerkleProof<F>, MerkleProof<F>)>> {
        let sample_indices_timer = start_timer!(|| "Sample indices");
        let mut indices =
            sample_indices(self.config.num_queries, self.config.L[0].len(), transcript);
        end_timer!(sample_indices_timer);

        let num_rounds = self.config.num_rounds();
        let mut openings = Vec::with_capacity(num_rounds);

        for i in 0..(num_rounds - 1) {
            let codeword = codewords[i].clone();
            let codeword_next = codewords[i + 1].clone();
            let tree = &trees[i];
            let tree_next = &trees[i + 1];

            reduce_indices(&mut indices, codeword_next.len());
            // Compute the degree-1 polynomial from the indices

            let mut openings_i = Vec::with_capacity(indices.len());
            for index in &indices {
                let alpha_0 = codeword[*index];
                let alpha_1 = codeword[index + (codeword.len() / 2)];
                let beta = codeword_next[*index];

                // Check that the folded polynomial equals the
                // random linear combination of the og polynomial.

                // Sanity check
                let s_0 = self.config.L[i][*index];
                let s_1 = self.config.L[i][index + (codeword.len() / 2)];
                let y = self.config.L[i + 1][*index];
                assert_eq!(y, s_0 * s_0);
                assert_eq!(y, s_1 * s_1);

                let alpha_0_opening = tree.open(alpha_0);
                let alpha_1_opening = tree.open(alpha_1);
                let beta_opening = tree_next.open(beta);

                let opening = (alpha_0_opening, alpha_1_opening, beta_opening);
                openings_i.push(opening);
            }
            openings.push(openings_i);
        }

        openings
    }

    pub fn batch_prove(
        &self,
        polys: &[UniPoly<F>],
        transcript: &mut Transcript<F>,
    ) -> BatchedPolyEvalProof<F> {
        let mut codeword = vec![F::zero(); self.config.L[0].len()];

        // Committed Merkle trees of the polynomials
        let mut f_trees = Vec::with_capacity(polys.len());

        // Store the codewords to open later
        let mut f_codewords = Vec::with_capacity(polys.len());

        let mut z = Vec::with_capacity(polys.len());
        let mut y = Vec::with_capacity(polys.len());
        // Compute the codeword of g(X)
        for i in 0..polys.len() {
            let poly = &polys[i];
            let k = (poly.degree() + 2).next_power_of_two();

            let evals = poly.eval_fft(&self.config.L[0]);
            f_codewords.push(evals.clone());

            let mut tree = MerkleTree::new();
            let comm = tree.commit(&evals);
            f_trees.push(tree);

            // Commit the codeword
            transcript.append_fe(&comm);

            // Get the challenge evaluation point
            let z_i = transcript.challenge_fe();
            let y_i = poly.eval(z_i);

            z.push(z_i);
            y.push(y_i);

            // Derive the codeword of the quotient polynomial
            let q_evals = evals
                .iter()
                .zip(self.config.L[0].iter())
                .map(|(f, x)| (*f - y_i) * (*x - z_i).invert().unwrap())
                .collect::<Vec<F>>();

            // Compute the challenges (i.e. random weights)
            let alpha = transcript.challenge_fe();
            let beta = transcript.challenge_fe();

            let pad_degree = k - poly.degree() - 1;

            let weighed_evals = q_evals
                .iter()
                .enumerate()
                .map(|(j, e)| {
                    let weighted = alpha * e
                        + beta * e * self.config.L[0][j].pow(&[pad_degree as u64, 0, 0, 0]);

                    weighted
                })
                .collect::<Vec<F>>();

            for (j, weighted_eval) in weighed_evals.iter().enumerate() {
                codeword[j] += weighted_eval;
            }
        }

        let commit_timer = start_timer!(|| "Commit");
        let (mut codewords, mut oracles) = self.commit(codeword.clone(), transcript);
        end_timer!(commit_timer);

        // TODO: We can do better!
        let mut unused_tree = MerkleTree::new();
        unused_tree.commit(&codeword);

        codewords.insert(0, codeword);
        oracles.insert(0, unused_tree);

        let query_timer = start_timer!(|| "Query");
        let queries = self.query(&codewords, &oracles, transcript);
        end_timer!(query_timer);

        let mut poly_openings = Vec::with_capacity(polys.len());

        // Open the polynomials at the indices
        for query in &queries[0] {
            let (alpha_0_opening, alpha_1_opening, _) = query;
            let mut openings = Vec::with_capacity(polys.len());
            for (i, tree) in f_trees.iter().enumerate() {
                let opening_alpha_0 = tree.open(f_codewords[i][alpha_0_opening.leaf_index]);
                let opening_alpha_1 = tree.open(f_codewords[i][alpha_1_opening.leaf_index]);
                openings.push((opening_alpha_0, opening_alpha_1));
            }
            poly_openings.push(openings);
        }

        BatchedPolyEvalProof {
            openings: queries,
            poly_openings,
            poly_degrees: polys.iter().map(|p| p.degree()).collect(),
            z,
            y,
            reduced_codeword: codewords[codewords.len() - 1].clone(),
        }
    }
}
