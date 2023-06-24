use crate::fri::sparse_unipoly::SparseUniPoly;
use crate::fri::tree::CommittedMerkleTree;
use crate::fri::unipoly::UniPoly;
use crate::fri::utils::{reduce_indices, sample_indices};
use crate::fri::{
    FirstRoundBoundedPolyOpenings, MLPolyEvalProof, Query, QueryFirstRound, RoundQueries,
    SparseFRIConfig,
};
use crate::spartan::polynomial::sparse_ml_poly::SparseMLPoly;
use crate::transcript::Transcript;
use crate::MultilinearPCS;
use crate::{FieldExt, SparseMultilinearPCS};
use ark_std::{end_timer, start_timer};
use ff::BatchInvert;

use super::sparse_unipoly::CommittedSparseUniPoly;

#[derive(Clone)]
pub struct SparseFRIMLPolyCommit<F>
where
    F: FieldExt,
{
    pub config: SparseFRIConfig<F>,
}

impl<F> SparseMultilinearPCS<F> for SparseFRIMLPolyCommit<F>
where
    F: FieldExt,
{
    type Commitment = [u8; 32];
    type Opening = MLPolyEvalProof<F>;
    type Config = SparseFRIConfig<F>;

    fn new(config: Self::Config) -> Self {
        Self { config }
    }

    // Commit to a multilinear polynomial
    fn commit(&self, poly: &SparseMLPoly<F>) -> Self::Commitment {
        // We commit to a univariate polynomials f(X) which has the same
        // coefficients as the given multilinear polynomial
        let mut val = self.config.domain_generator(0);
        let coeffs = (poly.coeffs.clone()).unwrap(); // TODO: Avoid cloning here
        let domain_size = coeffs.len() * self.config.expansion_factor;

        let mut evals = Vec::with_capacity(domain_size);
        for _ in 0..domain_size {
            let eval = poly.eval_as_uni_with_coeffs(&val);
            evals.push(eval);
            val *= val;
        }

        // TODO: Save this tree for later use?
        let tree = CommittedMerkleTree::from_leaves(evals);
        tree.root()
    }

    fn open(
        &self,
        poly: &SparseMLPoly<F>,
        point: &[F],
        transcript: &mut Transcript<F>,
    ) -> Self::Opening {
        let opening = self.prove_eval(poly, point, transcript);
        opening
    }

    fn verify(
        &self,
        opening: &Self::Opening,
        commitment: &Self::Commitment,
        transcript: &mut Transcript<F>,
    ) {
        self.verify_inner(commitment, opening, transcript);
    }
}

impl<F> SparseFRIMLPolyCommit<F>
where
    F: FieldExt,
{
    fn fold(&self, poly: &SparseUniPoly<F>, alpha: F) -> SparseUniPoly<F> {
        let (even_coeffs, odd_coeffs) = poly.to_even_and_odd_coeffs();

        // construct q_i(x) = (f_{i_e}(X) - a_i) / (X - \beta) + alpha * f_{i_o}(X)
        // Construct f_i(X) = g_{i-1}(X) + alpha * h{i-1}(X)
        let mut coeffs = vec![F::ZERO; poly.degree + 1];
        for e in even_coeffs {
            coeffs[e.0] += e.1;
        }

        for o in odd_coeffs {
            coeffs[o.0] += o.1 * alpha
        }

        let dense_coeffs = coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| (i, *c))
            .filter(|(_, c)| *c != F::ZERO)
            .collect();

        SparseUniPoly::new(dense_coeffs, poly.degree / 2)
    }

    // Prove that the evaluation of a multilinear polynomial f(x_0, ..., x_{m-1})
    // at the point (x_0, ..., x_{m-1})
    pub fn prove_eval(
        &self,
        ml_poly: &SparseMLPoly<F>,
        x: &[F],
        transcript: &mut Transcript<F>,
    ) -> MLPolyEvalProof<F> {
        // TODO: Avoid cloning here
        let mut ml_poly = ml_poly.clone();
        ml_poly.compute_coeffs();

        debug_assert!(ml_poly.coeffs.is_some());

        // Set the constants
        let m = x.len();

        // Evaluate the polynomial the point (x_0, ..., x_{m-1})
        let y = ml_poly.eval(x);

        // Commit to the univariate polynomial f(X) which has the
        // same coefficients as the given multilinear polynomial
        let compute_and_commit_f_timer = start_timer!(|| "Compute and commit f");
        let coeffs = ml_poly.coeffs.clone().unwrap();
        let f_0 = SparseUniPoly::new(coeffs, 2usize.pow(m as u32) - 1);
        let f_0_comm = f_0.merkle_commit(
            self.config.domain_generator(0),
            self.config.expansion_factor,
            transcript,
        );

        let mut f_comms = Vec::with_capacity(m);

        f_comms.push(f_0_comm.clone());

        let codeword = f_0_comm.committed_merkle_tree.leaves();

        // ##########################################
        // Compute and commitments to f_1,...,f_m
        // ##########################################

        let mut num_coeffs = 2usize.pow(m as u32);
        for i in 1..m {
            // Combine the even and odd coefficients by z;
            let (even_coeffs, odd_coeffs) = f_comms[i - 1].to_even_and_odd_coeffs();

            // Construct f_i(X) = g_{i-1}(X) + x[i] * h{i-1}(X)
            let mut coeffs = vec![F::ZERO; num_coeffs / 2];
            for e in even_coeffs {
                coeffs[e.0] += e.1;
            }

            for o in odd_coeffs {
                coeffs[o.0] += o.1 * x[i - 1];
            }

            let sparse_coeffs = coeffs
                .iter()
                .enumerate()
                .map(|(i, c)| (i, *c))
                .filter(|(_, c)| *c != F::ZERO)
                .collect::<Vec<(usize, F)>>();

            num_coeffs /= 2;
            let f_i = SparseUniPoly::new(sparse_coeffs, num_coeffs - 1);

            // Sanity check
            let f_i_prev = &f_comms[i - 1];
            debug_assert_eq!(f_i_prev.eval(x[i - 1]), f_i.eval(x[i - 1].square()));

            // Commit to f_i(X)
            let f_i_comm = f_i.merkle_commit(
                self.config.domain_generator(i),
                self.config.expansion_factor,
                transcript,
            );
            f_comms.push(f_i_comm);
        }
        end_timer!(compute_and_commit_f_timer);

        // Get the challenge point to evaluate f_0,f_1,...,f_{m-1}
        let beta = transcript.challenge_fe();
        debug_assert_eq!((-beta).square(), beta.square());
        let beta_squared = beta.square();
        let minus_beta = -beta;

        // evaluate f_0(x), f_1(x), ..., f_{m-1}(x) at \beta, -\beta

        let beta_powers_table = (0..(f_comms[0].degree + 1))
            .map(|i| beta.pow([i as u64, 0, 0, 0]))
            .collect::<Vec<F>>();

        let minus_beta_powers_table = (0..(f_comms[0].degree + 1))
            .map(|i| {
                if i & 1 == 0 {
                    beta_powers_table[i]
                } else {
                    -beta_powers_table[i]
                }
            })
            .collect::<Vec<F>>();

        let beta_squared_powers_table = (0..(f_comms[0].degree + 1))
            .map(|i| beta_squared.pow([i as u64, 0, 0, 0]))
            .collect::<Vec<F>>();

        let f_i_evals_beta = f_comms
            .iter()
            .map(|c| {
                (
                    c.eval_with_table(&beta_powers_table),
                    c.eval_with_table(&minus_beta_powers_table),
                )
            })
            .collect::<Vec<(F, F)>>();

        let f_i_evals_beta_squared = f_comms
            .pop()
            .iter()
            .map(|c| c.eval_with_table(&beta_squared_powers_table))
            .collect::<Vec<F>>();

        // ##########################################
        // Commit phase for g(X)
        // ##########################################

        let commit_timer = start_timer!(|| "Commit");
        let g_comms = self.commit_phase(&f_comms, transcript);
        end_timer!(commit_timer);

        // ##########################################
        // Query phase for g(X)
        // ##########################################

        let query_timer = start_timer!(|| "Query");
        let mut indices = sample_indices(
            self.config.num_queries,
            self.config.domain_order_in_round(1),
            self.config.final_codeword_size,
            transcript,
        );

        let f_trees = f_comms
            .iter()
            .map(|c| &c.committed_merkle_tree)
            .collect::<Vec<&CommittedMerkleTree<F>>>();
        let queries_first_round =
            self.query_first_round(&mut indices, &f_trees, &g_comms[0].committed_merkle_tree);

        let g_trees = g_comms
            .iter()
            .map(|c| &c.committed_merkle_tree)
            .collect::<Vec<&CommittedMerkleTree<F>>>();
        let queries = self.query(&mut indices, &g_trees, transcript);
        end_timer!(query_timer);

        let reduced_codeword = g_trees[g_trees.len() - 1].leaves();

        MLPolyEvalProof {
            queries_first_round,
            queries,
            y,
            reduced_codeword,
            f_i_evals_beta,
            f_i_evals_beta_squared,
            x: x.to_vec(),
        }
    }

    fn commit_phase(
        &self,
        committed_polys: &[CommittedSparseUniPoly<F>],
        transcript: &mut Transcript<F>,
    ) -> Vec<CommittedSparseUniPoly<F>> {
        let num_rounds = self.config.num_rounds();
        let mut committed_round_polys = Vec::<CommittedSparseUniPoly<F>>::with_capacity(num_rounds);

        let weight_alpha = transcript.challenge_fe();
        let weight_beta = transcript.challenge_fe();

        // This is the degree of f_0
        let poly_max_degree = committed_polys[0].degree;

        // g(x) will be the degree-corrected random linear combination of the committed polynomials
        let g_poly_degree = poly_max_degree.next_power_of_two();

        // TODO adjust degree of g(x) to be power of two

        // Compute the degree-corrected random  combination of the committed polynomials
        let mut g_poly_coeffs = vec![F::ZERO; poly_max_degree + 1];

        // Compute the DEEP polynomials for all the committed polynomials

        for committed_poly in committed_polys {
            for coeff in &committed_poly.coeffs {
                // TODO: Random linear combination
                g_poly_coeffs[coeff.0] += coeff.1;
            }
        }

        let g_poly_coeffs_dense = g_poly_coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| (i, *c))
            .filter(|(_, c)| *c != F::ZERO)
            .collect::<Vec<(usize, F)>>();

        let mut g_poly = SparseUniPoly::new(g_poly_coeffs_dense, poly_max_degree);

        for i in 0..num_rounds {
            let tree = g_poly.merkle_commit(
                self.config.domain_generator(i),
                self.config.expansion_factor,
                transcript,
            );
            committed_round_polys.push(tree);

            // "Fold" the polynomial
            let x = transcript.challenge_fe();
            g_poly = self.fold(&g_poly, x);
        }

        committed_round_polys
    }

    // query f_0(x), f_1(x), f_2(x), ..., f_m-1(x) at random indices
    pub fn query_first_round(
        &self,
        indices: &mut Vec<usize>,
        f_trees: &[&CommittedMerkleTree<F>],
        g_1: &CommittedMerkleTree<F>,
    ) -> Vec<QueryFirstRound<F>> {
        debug_assert_eq!(self.config.num_queries, indices.len());

        let mut queries = Vec::with_capacity(self.config.num_queries);

        for index in indices {
            let domain_generator = self.config.domain_generator(0);

            // s_0 = w^i
            let s_0_index = *index;
            let s_0 = domain_generator.pow(&[s_0_index as u64, 0, 0, 0]); // w^i

            // s_1 = w^{i + n} where n is the domain order of the next round
            let s_1_index = s_0_index + self.config.domain_order_in_round(1);
            let s_1 = domain_generator.pow(&[s_1_index as u64, 0, 0, 0]);

            // y = s_0^2 = s_1^2
            let y = s_0 * s_0;
            debug_assert_eq!(y, s_1 * s_1);

            // open f_0(x), f_1(x), ..., f_{m-1}(x) at s_0
            let f_openings = f_trees
                .iter()
                .map(|f_tree| FirstRoundBoundedPolyOpenings {
                    opening_at_s0: f_tree.open_index(s_0_index),
                    opening_at_s1: f_tree.open_index(s_1_index),
                })
                .collect::<Vec<FirstRoundBoundedPolyOpenings<F>>>();

            // open g_1(x) at y
            let g_opening = g_1.open(y);

            queries.push(QueryFirstRound {
                bounded_poly_openings: f_openings,
                g_1_opening_at_t: g_opening,
            })
        }

        queries
    }

    pub fn query(
        &self,
        indices: &mut Vec<usize>,
        g_trees: &[&CommittedMerkleTree<F>],
        transcript: &mut Transcript<F>,
    ) -> Vec<RoundQueries<F>> {
        let num_rounds = self.config.num_rounds();
        let mut queries = Vec::with_capacity(num_rounds);

        // We skip the first round because the verifier derives the opening
        // from the oracle of the polynomials that consists of g(X)
        for i in 0..(num_rounds - 2) {
            let round = i + 1; // Round 0 is done in query_first_round
            let round_domain_order = self.config.domain_order_in_round(round);

            let tree = &g_trees[i];
            let tree_next = &g_trees[i + 1];

            reduce_indices(indices, round_domain_order / 2);

            let openings_i = indices
                .iter()
                .map(|index| Query {
                    opening_at_s0: tree.open_index(*index),
                    opening_at_s1: tree.open_index(*index + (round_domain_order / 2)),
                    opening_at_t: tree_next.open_index(*index),
                })
                .collect();

            queries.push(RoundQueries {
                queries: openings_i,
            });
        }

        queries
    }

    pub fn verify_folding(&self, proof: &MLPolyEvalProof<F>, beta_challenge: F) {
        assert_eq!(proof.x.len(), proof.f_i_evals_beta.len());
        let two_inv = F::TWO_INV;
        let two_inv_beta = two_inv * beta_challenge.invert().unwrap();
        let n = proof.f_i_evals_beta.len();
        for i in 0..n {
            let (f_i_eval_beta, f_i_eval_minus_beta) = proof.f_i_evals_beta[i];

            let rhs = (f_i_eval_beta + f_i_eval_minus_beta) * two_inv
                + proof.x[i] * (f_i_eval_beta - f_i_eval_minus_beta) * two_inv_beta;

            if i != n - 1 {
                let f_i_eval_beta_squared = proof.f_i_evals_beta_squared[i];
                assert_eq!(f_i_eval_beta_squared, rhs);
            } else {
                assert_eq!(proof.y, rhs);
            }
        }
    }

    pub fn get_quotient_evaluations(y: F, f_eval: (F, F), x: F, f_x: (F, F)) -> (F, F) {
        let eval_1 = (f_eval.0 - y) * (f_x.0 - x).invert().unwrap();
        let eval_2 = (f_eval.1 - y) * (f_x.1 - x).invert().unwrap();
        (eval_1, eval_2)
    }

    pub fn get_padded_rlc_eval(y: F, x_powered: F, weight: (F, F)) -> F {
        weight.0 * y + weight.1 * y * x_powered
    }

    // Interpolate the first degree-1 round polynomial
    // from the Merkle opening proofs of the quotient polynomials
    // that represents the evaluation of f_0, f_1, ..., f_{m-1}
    // at \beta, -\beta, \beta^2.
    // The quotient polynomials are evaluated at \e0, \e1, and \beta.
    pub fn interpolate_round_0_poly_at_index(
        &self,
        f_i_evals_beta: &[(F, F)],
        f_i_evals_beta_squared: &[F],
        query_first_round: &QueryFirstRound<F>,
        index: usize,
        challenge_weights: &[(F, F)],
        beta_challenge: F,
    ) -> (F, F) {
        let bounded_poly_openings = &query_first_round.bounded_poly_openings;

        let round_0_domain_generator = self.config.domain_generator(0);
        let round_0_domain_order = self.config.domain_order_in_round(0);

        let mut poly_degree = round_0_domain_order / self.config.expansion_factor - 1;
        let k = (poly_degree + 1).next_power_of_two();

        // round_0_openings consists of all openings for all the queries

        let s_0 = round_0_domain_generator.pow(&[index as u64, 0, 0, 0]);
        let s_1 =
            round_0_domain_generator.pow(&[(index + (round_0_domain_order / 2)) as u64, 0, 0, 0]);

        // \e0 and e1 openings consists of all openings of the batched polynomials,
        // whereas \beta opening is a single opening to the batched polynomial.
        let mut weighted_alpha_0 = F::ZERO;
        let mut weighted_alpha_1 = F::ZERO;
        // First, we take the random linear combination of \e0 and \e1 openings

        assert_eq!(
            f_i_evals_beta.len() * 2 + f_i_evals_beta_squared.len(),
            challenge_weights.len()
        );
        assert_eq!(
            query_first_round.bounded_poly_openings.len(),
            f_i_evals_beta.len()
        );
        assert_eq!(f_i_evals_beta.len(), f_i_evals_beta_squared.len() + 1);

        let num_polys = bounded_poly_openings.len();

        for i in 0..num_polys {
            let opening_at_s0 = &bounded_poly_openings[i].opening_at_s0;
            let opening_at_s1 = &bounded_poly_openings[i].opening_at_s1;
            let (eval_beta, eval_minus_beta) = f_i_evals_beta[i];

            // e0 and e1 are the supposed evaluations of f_i(s_0) and f_i(s_1)
            // respectively.
            let e0 = opening_at_s0.leaf_fe();
            let e1 = opening_at_s1.leaf_fe();

            // ##########################################
            // Derive the evaluation of the quotient polynomial
            // ##########################################

            // Two evaluations for f_i(\beta)

            // q_beta_0 = f_i(\s_0) - a_eval / (\s_0 - \beta)
            // q_beta_1 = f_i(\s_1) - a_eval / (\s_1 - \beta)
            let (q_beta_0, q_beta_1) =
                Self::get_quotient_evaluations(eval_beta, (e0, e1), beta_challenge, (s_0, s_1));

            // Two evaluations for f_i(-\beta)
            // q_minus_beta_0 = f(\s_0) - b_eval / (\s_0 - \-beta)
            // q_minus_beta_1 = f(\s_1) - b_eval / (\s_1 - \-beta)
            let (q_minus_beta_0, q_minus_beta_1) = Self::get_quotient_evaluations(
                eval_minus_beta,
                (e0, e1),
                -beta_challenge,
                (s_0, s_1),
            );

            // Adjust the degree of the polynomial
            let pad_degree = k - poly_degree - 1;

            // x_pow_0 = w^{index}^{k-d-1}
            let x_pow_0 = round_0_domain_generator.pow(&[(index + pad_degree) as u64, 0, 0, 0]);
            // Compute w^{index + (domain.len() / 2)}^{k-d-1}
            let x_pow_1 = round_0_domain_generator.pow(&[
                (index + (round_0_domain_order / 2) + pad_degree) as u64,
                0,
                0,
                0,
            ]);

            let q_beta_weight = challenge_weights[i * 3];

            weighted_alpha_0 += Self::get_padded_rlc_eval(q_beta_0, x_pow_0, q_beta_weight);
            weighted_alpha_1 += Self::get_padded_rlc_eval(q_beta_1, x_pow_1, q_beta_weight);

            let q_minus_beta_weight = challenge_weights[i * 3 + 1];

            weighted_alpha_0 +=
                Self::get_padded_rlc_eval(q_minus_beta_0, x_pow_0, q_minus_beta_weight);

            weighted_alpha_1 +=
                Self::get_padded_rlc_eval(q_minus_beta_1, x_pow_1, q_minus_beta_weight);

            if i != num_polys - 1 {
                let eval_beta_squared = f_i_evals_beta_squared[i];

                // Two evaluations for f_{i+1}(\beta^2)
                let (q_beta_squared_0, q_beta_squared_1) = Self::get_quotient_evaluations(
                    eval_beta_squared,
                    (
                        bounded_poly_openings[i + 1].opening_at_s0.leaf_fe(),
                        bounded_poly_openings[i + 1].opening_at_s1.leaf_fe(),
                    ),
                    beta_challenge.square(),
                    (s_0, s_1),
                );

                let q_beta_squared_weight = challenge_weights[i * 3 + 2];

                weighted_alpha_0 +=
                    Self::get_padded_rlc_eval(q_beta_squared_0, x_pow_0, q_beta_squared_weight);

                weighted_alpha_1 +=
                    Self::get_padded_rlc_eval(q_beta_squared_1, x_pow_1, q_beta_squared_weight);
            }

            poly_degree /= 2;
        }

        // Interpolate the points (s_0, e0) (s_1, e1)
        let coeff = (weighted_alpha_0 - weighted_alpha_1) * (s_0 - s_1).invert().unwrap();
        let intercept = weighted_alpha_0 - coeff * s_0;

        (coeff, intercept)
    }

    pub fn verify_inner(
        &self,
        comm: &[u8; 32],
        proof: &MLPolyEvalProof<F>,
        transcript: &mut Transcript<F>,
    ) {
        let num_rounds = self.config.num_rounds();

        // Append all commitments to f_0,...f_{n-1} polynomials to the transcript
        for openings in &proof.queries_first_round[0].bounded_poly_openings {
            // Append the commitment to the polynomial to the transcript
            transcript.append_bytes(&openings.opening_at_s0.root);
        }

        let beta_challenge = transcript.challenge_fe();
        self.verify_folding(proof, beta_challenge);

        // Compute the challenges of the weighted sum
        let challenge_weights = (0..(proof.f_i_evals_beta.len() * 2
            + proof.f_i_evals_beta_squared.len()))
            .map(|_| (transcript.challenge_fe(), transcript.challenge_fe()))
            .collect::<Vec<(F, F)>>();

        let mut challenges = Vec::with_capacity(num_rounds);
        challenges.push(transcript.challenge_fe());

        for i in 0..(num_rounds - 2) {
            // Append the commitment to the batched quotient polynomial to the transcript,
            // and get the challenge value
            transcript.append_bytes(&proof.queries[i].queries[0].opening_at_s0.root);
            challenges.push(transcript.challenge_fe());

            if i == num_rounds - 3 {
                transcript.append_bytes(&proof.queries[i].queries[0].opening_at_t.root);
                challenges.push(transcript.challenge_fe());
            }
        }

        let final_domain_generator = self.config.domain_generator(num_rounds - 1);
        let reduced_domain = (0..self.config.final_codeword_size)
            .map(|i| final_domain_generator.pow(&[i as u64, 0, 0, 0]))
            .collect::<Vec<F>>();
        debug_assert_eq!(reduced_domain.len(), proof.reduced_codeword.len());

        // Check that the final codeword corresponds to a low-degree polynomial
        let interpolant = UniPoly::interpolate(&reduced_domain, &proof.reduced_codeword);
        let target_degree = self.config.final_codeword_size / self.config.expansion_factor - 1;
        assert_eq!(interpolant.degree(), target_degree);

        let mut indices = sample_indices(
            self.config.num_queries,
            self.config.domain_order_in_round(1),
            self.config.final_codeword_size,
            transcript,
        );

        reduce_indices(&mut indices, self.config.domain_order_in_round(1));

        for (index, query) in indices.iter().zip(&proof.queries_first_round) {
            let (coeff, intercept) = self.interpolate_round_0_poly_at_index(
                &proof.f_i_evals_beta,
                &proof.f_i_evals_beta_squared,
                &query,
                *index,
                &challenge_weights,
                beta_challenge,
            );

            let g_i_opening_at_t = &query.g_1_opening_at_t;
            let beta = g_i_opening_at_t.leaf_fe();

            assert_eq!(
                beta,
                coeff * challenges[0] + intercept,
                "failed at round 0 ",
            );

            for (i, openings) in query.bounded_poly_openings.iter().enumerate() {
                let opening_at_s0 = &openings.opening_at_s0;
                let opening_at_s1 = &openings.opening_at_s1;
                assert!(opening_at_s0.verify());
                assert!(opening_at_s1.verify());

                if i == 0 {
                    assert_eq!(&opening_at_s0.root, comm);
                    assert_eq!(&opening_at_s1.root, comm);
                } else {
                    // Check that the e0 and e1 come from the same oracle
                    assert_eq!(opening_at_s0.root, opening_at_s1.root);
                }
            }

            assert!(g_i_opening_at_t.verify());
        }

        for i in 0..(num_rounds - 2) {
            let round = i + 1;
            let domain_generator = self.config.domain_generator(round);
            let domain_order = self.config.domain_order_in_round(round);
            reduce_indices(&mut indices, domain_order / 2);
            for (j, openings) in proof.queries[i].queries.iter().enumerate() {
                let opening_at_s0 = &openings.opening_at_s0;
                let opening_at_s1 = &openings.opening_at_s1;
                let opening_at_t = &openings.opening_at_t;

                let index = indices[j];
                let s_0 = domain_generator.pow(&[index as u64, 0, 0, 0]);
                let s_1 = domain_generator.pow(&[(index + domain_order / 2) as u64, 1, 0, 0]);

                let e0 = opening_at_s0.leaf_fe();
                let e1 = opening_at_s1.leaf_fe();
                let e3 = opening_at_t.leaf_fe();

                // Interpolate the points (s_0, e0) (s_1, e1)
                let coeff = (e0 - e1) * (s_0 - s_1).invert().unwrap();
                let intercept = e0 - coeff * s_0;

                assert_eq!(
                    e3,
                    coeff * challenges[i + 1] + intercept,
                    "failed at round = {}",
                    i
                );

                assert!(opening_at_s0.verify());
                assert!(opening_at_s1.verify());
                assert!(opening_at_t.verify());

                // Check that the e0 and e1 come from the same oracle
                assert_eq!(opening_at_s0.root, opening_at_s1.root);
            }
        }
    }
}
