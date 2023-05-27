use crate::fri::tree::{CommittedMerkleTree, MerkleProof};
use crate::fri::unipoly::UniPoly;
use crate::fri::utils::{reduce_indices, sample_indices};
use crate::fri::{
    FRIConfig, FirstRoundBoundedPolyOpenings, MLPolyEvalProof, Query, QueryFirstRound, RoundQueries,
};
use crate::spartan::polynomial::ml_poly::MlPoly;
use crate::transcript::Transcript;
use crate::PolyCommitment;
use ark_std::{end_timer, start_timer};
use pasta_curves::arithmetic::FieldExt;

#[derive(Clone)]
pub struct FRIMLPolyCommitProver<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub config: FRIConfig<F>,
}

impl<F> PolyCommitment<F> for FRIMLPolyCommitProver<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    type Commitment = Vec<F>;
    type Opening = Vec<F>;
    fn new() -> Self {
        todo!();
    }

    fn commit(&self, evals: &[F]) -> Self::Commitment {
        todo!();
    }

    fn open(&self, point: &[F]) -> Self::Opening {
        todo!();
    }
}

impl<F> FRIMLPolyCommitProver<F>
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

    // Prove that the evaluation of a multilinear polynomial f(x_0, ..., x_{m-1})
    // at the point (x_0, ..., x_{m-1})
    pub fn prove_eval(
        &self,
        ml_poly: &MlPoly<F>,
        x: &[F],
        transcript: &mut Transcript<F>,
    ) -> MLPolyEvalProof<F> {
        assert!(ml_poly.coeffs.is_some());

        // Set the constants
        let m = x.len();

        // Evaluate the polynomial the point (x_0, ..., x_{m-1})
        let y = ml_poly.eval(x);

        // Commit to the univariate polynomial f(X) which has the
        // same coefficients as the given multilinear polynomial
        let f_0 = UniPoly::new(ml_poly.coeffs.as_ref().unwrap().clone());
        let f_0_comm = f_0.merkle_commit(&self.config.L[0], transcript);

        let mut f_comms = Vec::with_capacity(m);

        f_comms.push(f_0_comm.clone());

        let codeword = f_0_comm.committed_merkle_tree.leaves();
        let target_codeword_size = self.config.L[0].len();

        assert_eq!(codeword.len(), target_codeword_size);

        // ##########################################
        // Compute and commitments to f_1,...,f_m
        // ##########################################

        for i in 1..m {
            // Combine the even and odd coefficients by z;
            let (even_coeffs, odd_coeffs) = f_comms[i - 1].to_even_and_odd_coeffs();

            // Construct f_i(X) = g_{i-1}(X) * x[i] * h{i-1}(X)
            let coeffs = even_coeffs
                .iter()
                .zip(odd_coeffs.iter())
                .map(|(e, o)| *e + *o * x[i])
                .collect::<Vec<F>>();
            let f_i = UniPoly::new(coeffs);

            // Sanity check
            let f_i_prev = &f_comms[i - 1];
            assert_eq!(f_i_prev.eval(x[i]), f_i.eval(x[i].square()));

            // Commit to f_i(X)
            let f_i_comm = f_i.merkle_commit(&self.config.L[0], transcript);
            f_comms.push(f_i_comm);
        }

        // Get the challenge point to evaluate f_0,f_1,...,f_{m-1}
        let beta = transcript.challenge_fe();
        assert_eq!((-beta).square(), beta.square());

        // Compute quotient codewords for all evaluation challenges
        let mut quotient_codewords = Vec::with_capacity(f_comms.len() * 3);

        let mut f_i_evals_beta = Vec::with_capacity(m);
        let mut f_i_evals_beta_squared = Vec::with_capacity(m - 1);

        for i in 0..m {
            let f_i = &f_comms[i];

            // Evaluate the polynomial at beta, -beta, and beta^2
            let f_i_eval_beta = f_i.eval(beta);
            let f_i_eval_minus_beta = f_i.eval(-beta);

            // Derive the codewords of the quotient polynomials
            // about the evaluations above
            quotient_codewords.push(Self::to_quotient_codeword(
                &f_i.codeword(),
                &self.config.L[0],
                beta,
                f_i_eval_beta,
            ));

            quotient_codewords.push(Self::to_quotient_codeword(
                &f_i.codeword(),
                &self.config.L[0],
                -beta,
                f_i_eval_minus_beta,
            ));

            f_i_evals_beta.push((f_i_eval_beta, f_i_eval_minus_beta));

            if i != m - 1 {
                let f_i_next = &f_comms[i + 1];
                let f_i_next_eval_beta_squared = f_i_next.eval(beta.square());
                quotient_codewords.push(Self::to_quotient_codeword(
                    &f_i_next.codeword(),
                    &self.config.L[0],
                    beta.square(),
                    f_i_next_eval_beta_squared,
                ));
                f_i_evals_beta_squared.push(f_i_next_eval_beta_squared);
            }
        }

        // ##########################################
        // Compute the codeword of g(X),
        // which is the random linear combination of the quotient codewords
        // ##########################################
        let mut g_codeword = vec![F::zero(); self.config.L[0].len()];

        let quotient_poly_max_degree = ml_poly.evals.len() - 1;
        let k = quotient_poly_max_degree.next_power_of_two();

        let mut poly_degree = (quotient_codewords[0].len() / self.config.expansion_factor) - 1;

        for i in 0..quotient_codewords.len() {
            // Compute the challenges (i.e. random weights)
            let alpha = transcript.challenge_fe();
            let beta = transcript.challenge_fe();

            let pad_degree = k - poly_degree - 1;

            assert_eq!(quotient_codewords[i].len(), target_codeword_size);

            let weighed_evals = quotient_codewords[i]
                .iter()
                .enumerate()
                .map(|(j, e)| {
                    let weighted = alpha * e
                        + beta * e * self.config.L[0][j].pow(&[pad_degree as u64, 0, 0, 0]);

                    weighted
                })
                .collect::<Vec<F>>();

            for (j, weighted_eval) in weighed_evals.iter().enumerate() {
                g_codeword[j] += weighted_eval;
            }

            if i % 3 == 2 {
                poly_degree /= 2;
            }
        }

        // ##########################################
        // Commit phase for g(X)
        // ##########################################

        let commit_timer = start_timer!(|| "Commit");
        let (codewords, oracles) = self.commit(g_codeword.clone(), transcript);
        end_timer!(commit_timer);

        assert_eq!(
            codewords.len(),
            (g_codeword.len() as f64).log2() as usize - 1
        );

        // ##########################################
        // Query phase for g(X)
        // ##########################################

        let query_timer = start_timer!(|| "Query");
        let sample_indices_timer = start_timer!(|| "Sample indices");
        let mut indices = sample_indices(
            self.config.num_queries,
            self.config.L[0].len(),
            self.config.final_codeword_size(),
            transcript,
        );
        end_timer!(sample_indices_timer);
        reduce_indices(&mut indices, codewords[0].len());

        let queries_first_round = self.query_first_round(
            &mut indices,
            &codeword,
            &codewords[0],
            &f_comms
                .iter()
                .map(|f_comm| f_comm.committed_merkle_tree.clone())
                .collect::<Vec<CommittedMerkleTree<F>>>(),
            &oracles[0],
        );

        let queries = self.query(&mut indices, &codewords, &oracles, transcript);
        end_timer!(query_timer);

        MLPolyEvalProof {
            queries_first_round,
            queries,
            y,
            reduced_codeword: codewords[codewords.len() - 1].clone(),
            f_i_evals_beta,
            f_i_evals_beta_squared,
            x: x.to_vec(),
        }
    }

    fn to_quotient_codeword(codeword: &[F], domain: &[F], z: F, y: F) -> Vec<F> {
        codeword
            .iter()
            .zip(domain.iter())
            .map(|(f, x)| (*f - y) * (*x - z).invert().unwrap())
            .collect::<Vec<F>>()
    }

    pub fn commit(
        &self,
        codeword: Vec<F>,
        transcript: &mut Transcript<F>,
    ) -> (Vec<Vec<F>>, Vec<CommittedMerkleTree<F>>) {
        let num_rounds = self.config.num_rounds();
        let mut trees = Vec::<CommittedMerkleTree<F>>::with_capacity(num_rounds - 1);
        let mut codewords = Vec::<Vec<F>>::with_capacity(num_rounds);

        let mut codeword = codeword;
        let mut challenges = Vec::with_capacity(num_rounds);
        for i in 0..num_rounds {
            if i != 0 {
                // For i = 0, the verifier derives the opening from the oracle
                // of f(X), so we don't need to commit to it.
                let tree = CommittedMerkleTree::from_leaves(codeword.clone());
                transcript.append_fe(&tree.root());
                trees.push(tree);
            }

            // "Fold" the codeword
            let x = transcript.challenge_fe();
            challenges.push(x);
            codeword = self.fold(&codeword, &self.config.L[i], x);
            codewords.push(codeword.clone());
        }

        (codewords, trees)
    }

    pub fn query_first_round(
        &self,
        indices: &mut Vec<usize>,
        codeword: &[F],
        codeword_next: &[F],
        f_trees: &[CommittedMerkleTree<F>],
        tree_next: &CommittedMerkleTree<F>,
    ) -> Vec<QueryFirstRound<F>> {
        assert_eq!(codeword.len(), self.config.L[0].len());
        assert_eq!(codeword_next.len(), self.config.L[1].len());
        assert_eq!(self.config.num_queries, indices.len());

        let mut queries = Vec::with_capacity(self.config.num_queries);

        // Compute the degree-1 polynomial from the indices

        for index in indices {
            let e3 = codeword_next[*index];

            // Check that the folded polynomial equals the
            // random linear combination of the og polynomial.

            // Sanity check
            let s_0 = self.config.L[0][*index];
            let s_1 = self.config.L[0][*index + (codeword.len() / 2)];
            let y = self.config.L[1][*index];
            assert_eq!(y, s_0 * s_0);
            assert_eq!(y, s_1 * s_1);

            let mut openings_at_s0 = Vec::with_capacity(f_trees.len());
            let mut openings_at_s1 = Vec::with_capacity(f_trees.len());
            for f_tree in f_trees {
                let opening_at_s0 = f_tree.open_index(*index);
                let opening_at_s1 = f_tree.open_index(*index + (codeword.len() / 2));
                openings_at_s0.push(opening_at_s0);
                openings_at_s1.push(opening_at_s1);
            }

            let g_1_opening_at_t = tree_next.open(e3);

            let bounded_poly_openings = (openings_at_s0.iter().zip(openings_at_s1.iter()))
                .map(
                    |(opening_at_s0, opening_at_s1)| FirstRoundBoundedPolyOpenings {
                        opening_at_s0: opening_at_s0.clone(),
                        opening_at_s1: opening_at_s1.clone(),
                    },
                )
                .collect();

            queries.push(QueryFirstRound {
                bounded_poly_openings,
                g_1_opening_at_t,
            });
        }

        queries
    }

    pub fn query(
        &self,
        indices: &mut Vec<usize>,
        codewords: &[Vec<F>],
        trees: &[CommittedMerkleTree<F>],
        transcript: &mut Transcript<F>,
    ) -> Vec<RoundQueries<F>> {
        let num_rounds = self.config.num_rounds();
        let mut queries = Vec::with_capacity(num_rounds);

        // We skip the first round because the verifier derives the opening
        // from the oracle of the polynomials that consists of g(X)
        for i in 0..(num_rounds - 2) {
            let codeword = codewords[i].clone();
            let codeword_next = codewords[i + 1].clone();
            let tree = &trees[i];
            let tree_next = &trees[i + 1];

            reduce_indices(indices, codeword_next.len());
            // Compute the degree-1 polynomial from the indices

            let mut openings_i = Vec::with_capacity(indices.len());
            for index in indices.clone() {
                let e0 = codeword[index];
                let e1 = codeword[index + (codeword.len() / 2)];
                let e3 = codeword_next[index];

                // Check that the folded polynomial equals the
                // random linear combination of the og polynomial.

                // Sanity check
                let s_0 = self.config.L[i + 1][index];
                let s_1 = self.config.L[i + 1][index + (codeword.len() / 2)];
                let y = self.config.L[i + 2][index];

                assert_eq!(y, s_0 * s_0);
                assert_eq!(y, s_1 * s_1);

                let opening_at_s0 = tree.open(e0);
                let opening_at_s1 = tree.open(e1);
                let opening_at_t = tree_next.open(e3);

                openings_i.push(Query {
                    opening_at_s0,
                    opening_at_s1,
                    opening_at_t,
                });
            }
            queries.push(RoundQueries {
                queries: openings_i,
            });
        }

        queries
    }
}
