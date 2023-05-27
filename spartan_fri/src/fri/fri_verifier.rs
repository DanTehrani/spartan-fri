#![allow(non_snake_case)]

use crate::fri::unipoly::UniPoly;
use crate::fri::utils::{reduce_indices, sample_indices};
use crate::fri::{FRIConfig, MLPolyEvalProof, QueryFirstRound};
use crate::transcript::Transcript;
use pasta_curves::arithmetic::FieldExt;

pub struct FRIMLPolyCommitVerifier<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub config: FRIConfig<F>,
}

impl<F> FRIMLPolyCommitVerifier<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub fn new(config: FRIConfig<F>) -> Self {
        Self { config }
    }

    pub fn verify_folding(&self, proof: &MLPolyEvalProof<F>, beta_challenge: F) {
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
        let f_1_beta_squared_opening = &query_first_round.g_1_opening_at_t;

        let mut poly_degree = self.config.L[0].len() / self.config.expansion_factor - 1;
        let k = (poly_degree + 1).next_power_of_two();

        // round_0_openings consists of all openings for all the queries

        let domain = &self.config.L[0];

        let s_0 = domain[index];
        let s_1 = domain[index + (domain.len() / 2)];

        // Sanity check
        assert_eq!(self.config.L[1][index], s_0 * s_0);

        // \e0 and e1 openings consists of all openings of the batched polynomials,
        // whereas \beta opening is a single opening to the batched polynomial.
        let mut weighted_alpha_0 = F::zero();
        let mut weighted_alpha_1 = F::zero();
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
            let e0 = opening_at_s0.leaf;
            let e1 = opening_at_s1.leaf;

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
            let x_pow_0 = domain[index].pow(&[pad_degree as u64, 0, 0, 0]);
            // Compute w^{index + (domain.len() / 2)}^{k-d-1}
            let x_pow_1 = domain[index + (domain.len() / 2)].pow(&[pad_degree as u64, 0, 0, 0]);

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
                        bounded_poly_openings[i + 1].opening_at_s0.leaf,
                        bounded_poly_openings[i + 1].opening_at_s1.leaf,
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

    pub fn verify(&self, proof: &MLPolyEvalProof<F>, transcript: &mut Transcript<F>) {
        let num_rounds = self.config.num_rounds();

        // Append all commitments to f_0,...f_{n-1} polynomials to the transcript
        for openings in &proof.queries_first_round[0].bounded_poly_openings {
            // Append the commitment to the polynomial to the transcript
            transcript.append_fe(&openings.opening_at_s0.root);
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
            transcript.append_fe(&proof.queries[i].queries[0].opening_at_s0.root);
            challenges.push(transcript.challenge_fe());

            if i == num_rounds - 3 {
                transcript.append_fe(&proof.queries[i].queries[0].opening_at_t.root);
                challenges.push(transcript.challenge_fe());
            }
        }

        let reduced_domain = &self.config.L[self.config.L.len() - 1];
        assert_eq!(reduced_domain.len(), proof.reduced_codeword.len());

        // Check that the final codeword corresponds to a low-degree polynomial
        let interpolant = UniPoly::interpolate(&reduced_domain, &proof.reduced_codeword);
        let target_degree = self.config.final_codeword_size() / self.config.expansion_factor - 1;
        assert_eq!(interpolant.degree(), target_degree);

        let mut indices = sample_indices(
            self.config.num_queries,
            self.config.L[0].len(),
            self.config.final_codeword_size(),
            transcript,
        );

        reduce_indices(&mut indices, self.config.L[1].len());

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
            let beta = g_i_opening_at_t.leaf;

            assert_eq!(
                beta,
                coeff * challenges[0] + intercept,
                "failed at round 0 ",
            );

            for openings in &query.bounded_poly_openings {
                let opening_at_s0 = &openings.opening_at_s0;
                let opening_at_s1 = &openings.opening_at_s1;
                assert!(opening_at_s0.verify());
                assert!(opening_at_s1.verify());
                // Check that the e0 and e1 come from the same oracle
                assert_eq!(opening_at_s0.root, opening_at_s1.root);
            }

            assert!(g_i_opening_at_t.verify());
        }

        for i in 0..(num_rounds - 2) {
            reduce_indices(&mut indices, self.config.L[i + 2].len());
            for (j, openings) in proof.queries[i].queries.iter().enumerate() {
                let opening_at_s0 = &openings.opening_at_s0;
                let opening_at_s1 = &openings.opening_at_s1;
                let opening_at_t = &openings.opening_at_t;

                let L_i = &self.config.L[i + 1];
                let L_i_next = &self.config.L[i + 2];
                let index = indices[j];
                let s_0 = L_i[index];
                let s_1 = L_i[index + (L_i.len() / 2)];
                let y = L_i_next[index];
                assert_eq!(y, s_0 * s_0);
                assert_eq!(y, s_1 * s_1);

                let e0 = opening_at_s0.leaf;
                let e1 = opening_at_s1.leaf;
                let e3 = opening_at_t.leaf;

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
