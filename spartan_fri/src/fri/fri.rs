use crate::fri::tree::CommittedMerkleTree;
use crate::fri::unipoly::UniPoly;
use crate::fri::utils::{reduce_indices, sample_indices};
use crate::fri::{
    FRIConfig, FirstRoundBoundedPolyOpenings, MLPolyEvalProof, Query, QueryFirstRound, RoundQueries,
};
use crate::spartan::polynomial::ml_poly::MlPoly;
use crate::transcript::Transcript;
use crate::FieldExt;
use crate::MultilinearPCS;
use ark_std::{end_timer, start_timer};
use ff::BatchInvert;

#[derive(Clone)]
pub struct FRIMLPolyCommit<F>
where
    F: FieldExt,
{
    pub config: FRIConfig<F>,
}

impl<F> MultilinearPCS<F> for FRIMLPolyCommit<F>
where
    F: FieldExt,
{
    type Commitment = [u8; 32];
    type Opening = MLPolyEvalProof<F>;
    type Config = FRIConfig<F>;

    fn new(config: Self::Config) -> Self {
        Self { config }
    }

    fn commit(&self, poly: &MlPoly<F>) -> Self::Commitment {
        // Compute the evaluations of the polynomial at the domain
        let poly = UniPoly::new(poly.to_coeffs());
        let evals = poly.eval_fft(&self.config.L[0]);
        // TODO: Save this tree for later use
        let tree = CommittedMerkleTree::from_leaves(evals);

        tree.root()
    }

    fn open(&self, poly: &MlPoly<F>, point: &[F], transcript: &mut Transcript<F>) -> Self::Opening {
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

impl<F> FRIMLPolyCommit<F>
where
    F: FieldExt,
{
    fn fold(&self, codeword: &[F], domain: &[F], alpha: F) -> Vec<F> {
        debug_assert!(codeword.len() == domain.len());
        let two_inv = F::from(2).invert().unwrap();
        let one = F::from(1);

        let n = domain.len();

        let folded_codeword_size = n / self.config.folding_factor;
        let mut folded_codeword = Vec::with_capacity(folded_codeword_size);
        for i in 0..folded_codeword_size {
            // f*(w^2i) = 1/2 * ((1 + alpha * w^-i) * f(w^i) + (1 - alpha * w^-i) * f(-w^i))
            // w^(n/2) = -1
            // -w^i = domain[i + n/2]

            let omega_pow_minus_i = domain[i].invert().unwrap();

            let f_star_eval = two_inv
                * ((one + alpha * omega_pow_minus_i) * codeword[i]
                    + (one - alpha * omega_pow_minus_i) * codeword[i + folded_codeword_size]);
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
        // TODO: Avoid cloning here
        let mut ml_poly = ml_poly.clone();
        ml_poly.compute_coeffs();

        debug_assert!(ml_poly.coeffs.is_some());

        // Set the constants
        let m = x.len();

        // Evaluate the polynomial the point (x_0, ..., x_{m-1})
        let eval_timer = start_timer!(|| "Eval");
        let y = ml_poly.eval(x);
        end_timer!(eval_timer);

        // Commit to the univariate polynomial f(X) which has the
        // same coefficients as the given multilinear polynomial
        let compute_and_commit_f_timer = start_timer!(|| "Compute and commit f");
        let f_0 = UniPoly::new(ml_poly.coeffs.as_ref().unwrap().clone());
        let f_0_comm = f_0.merkle_commit(&self.config.L[0], transcript);

        let mut f_comms = Vec::with_capacity(m);

        f_comms.push(f_0_comm.clone());

        let codeword = f_0_comm.committed_merkle_tree.leaves();
        let target_codeword_size = self.config.L[0].len();

        debug_assert_eq!(codeword.len(), target_codeword_size);

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
                .map(|(e, o)| *e + *o * x[i - 1])
                .collect::<Vec<F>>();
            let f_i = UniPoly::new(coeffs);

            // Sanity check
            let f_i_prev = &f_comms[i - 1];
            debug_assert_eq!(f_i_prev.eval(x[i - 1]), f_i.eval(x[i - 1].square()));

            // Commit to f_i(X)
            let f_i_comm = f_i.merkle_commit(&self.config.L[0], transcript);
            f_comms.push(f_i_comm);
        }
        end_timer!(compute_and_commit_f_timer);

        // Get the challenge point to evaluate f_0,f_1,...,f_{m-1}
        let beta = transcript.challenge_fe();
        debug_assert_eq!((-beta).square(), beta.square());
        let beta_squared = beta.square();
        let minus_beta = -beta;

        // Compute quotient codewords for all evaluation challenges
        let mut quotient_codewords = Vec::with_capacity(f_comms.len() * 3);

        let mut f_i_evals_beta = Vec::with_capacity(m);
        let mut f_i_evals_beta_squared = Vec::with_capacity(m - 1);

        let compute_quotient_codewords_timer = start_timer!(|| "Compute quotient codewords");
        for i in 0..m {
            let f_i = &f_comms[i];

            // Evaluate the polynomial at beta, -beta, and beta^2
            let f_i_eval_beta = f_i.eval(beta);
            let f_i_eval_minus_beta = f_i.eval(minus_beta);

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
                minus_beta,
                f_i_eval_minus_beta,
            ));

            f_i_evals_beta.push((f_i_eval_beta, f_i_eval_minus_beta));

            if i != m - 1 {
                let f_i_next = &f_comms[i + 1];
                let f_i_next_eval_beta_squared = f_i_next.eval(beta_squared);
                quotient_codewords.push(Self::to_quotient_codeword(
                    &f_i_next.codeword(),
                    &self.config.L[0],
                    beta_squared,
                    f_i_next_eval_beta_squared,
                ));
                f_i_evals_beta_squared.push(f_i_next_eval_beta_squared);
            }
        }
        end_timer!(compute_quotient_codewords_timer);

        // ##########################################
        // Compute the codeword of g(X),
        // which is the random linear combination of the quotient codewords
        // ##########################################
        let mut g_codeword = vec![F::ZERO; self.config.L[0].len()];

        let quotient_poly_max_degree = ml_poly.evals.len() - 1;
        let k = quotient_poly_max_degree.next_power_of_two();

        let mut poly_degree = (quotient_codewords[0].len() / self.config.expansion_factor) - 1;

        let combine_quotient_codewords_timer = start_timer!(|| "Combine quotient codewords");
        for i in 0..quotient_codewords.len() {
            // Compute the challenges (i.e. random weights)
            let alpha = transcript.challenge_fe();
            let beta = transcript.challenge_fe();

            let pad_degree = k - poly_degree - 1;

            debug_assert_eq!(quotient_codewords[i].len(), target_codeword_size);

            let weighed_evals = quotient_codewords[i]
                .iter()
                .enumerate()
                .map(|(j, e)| {
                    // TODO: Optimize this
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
        end_timer!(combine_quotient_codewords_timer);

        // ##########################################
        // Commit phase for g(X)
        // ##########################################

        let commit_timer = start_timer!(|| "Commit");
        let (codewords, oracles) = self.commit_phase(g_codeword.clone(), transcript);
        end_timer!(commit_timer);

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
        let mut denominators = domain.iter().map(|x| *x - z).collect::<Vec<F>>();
        denominators.iter_mut().batch_invert();

        codeword
            .iter()
            .enumerate()
            .map(|(i, f)| (*f - y) * denominators[i])
            .collect::<Vec<F>>()
    }

    fn commit_phase(
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
                let tree = CommittedMerkleTree::<F>::from_leaves(codeword.clone());
                transcript.append_bytes(&tree.root());
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
        debug_assert_eq!(codeword.len(), self.config.L[0].len());
        debug_assert_eq!(codeword_next.len(), self.config.L[1].len());
        debug_assert_eq!(self.config.num_queries, indices.len());

        let mut queries = Vec::with_capacity(self.config.num_queries);

        // Compute the degree-1 polynomial from the indices

        let next_layer_size = codeword_next.len();
        let folding_factor = self.config.folding_factor as u64;

        for index in indices {
            let e3 = codeword_next[*index];

            // Check that the folded polynomial equals the
            // random linear combination of the og polynomial.

            // Sanity check
            let s_0 = self.config.L[0][*index];
            let s_1 = self.config.L[0][*index + next_layer_size];
            let y = self.config.L[1][*index];
            debug_assert_eq!(y, s_0.pow([folding_factor, 0, 0, 0]));
            debug_assert_eq!(y, s_1.pow([folding_factor, 0, 0, 0]));

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
        let folding_factor = self.config.folding_factor as u64;
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

                debug_assert_eq!(y, s_0.pow([folding_factor, 0, 0, 0]));
                debug_assert_eq!(y, s_1.pow([folding_factor, 0, 0, 0]));

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

        let reduced_domain = &self.config.L[self.config.L.len() - 1];
        debug_assert_eq!(reduced_domain.len(), proof.reduced_codeword.len());

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
