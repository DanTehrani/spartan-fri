#![allow(non_snake_case)]

use super::unipoly::UniPoly;
use super::utils::{reduce_indices, sample_indices};
use super::{BatchedPolyEvalProof, FRIConfig};
use crate::transcript::Transcript;
use pasta_curves::arithmetic::FieldExt;

pub struct FRIPolyCommitVerifier<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    config: FRIConfig<F>,
}

impl<F> FRIPolyCommitVerifier<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub fn new(config: FRIConfig<F>) -> Self {
        Self { config }
    }

    pub fn batch_verify(&self, proof: &BatchedPolyEvalProof<F>, transcript: &mut Transcript<F>) {
        // Compute the challenges of the weighted sum
        let num_rounds = self.config.num_rounds();
        let mut challenge_weights = Vec::with_capacity(proof.poly_degrees.len());

        for i in 0..proof.poly_openings[0].len() {
            // Append the commitment to the polynomial to the transcript
            transcript.append_fe(&proof.poly_openings[0][i].0.root);

            // Get the evaluation challenges
            let _z = transcript.challenge_fe();
            assert_eq!(proof.z[i], _z);

            // Get the random weights
            challenge_weights.push((transcript.challenge_fe(), transcript.challenge_fe()))
        }

        let mut challenges = Vec::with_capacity(num_rounds);
        for i in 0..(num_rounds - 1) {
            // Append the commitment to the batched quotient polynomial to the transcript,
            // and get the challenge value
            if i != 0 {
                transcript.append_fe(&proof.openings[i][0].1.root);
            }

            challenges.push(transcript.challenge_fe());
        }

        transcript.append_fe(&proof.openings[num_rounds - 2][0].2.root);
        let _x = transcript.challenge_fe();

        let reduced_domain = &self.config.L[self.config.L.len() - 1];
        assert_eq!(reduced_domain.len(), proof.reduced_codeword.len());

        // Check that the final codeword corresponds to a low-degree polynomial
        let interpolant = UniPoly::interpolate(&reduced_domain, &proof.reduced_codeword);
        let degree = (reduced_domain.len() / self.config.R) - 1;
        assert_eq!(interpolant.degree(), degree);

        let mut indices =
            sample_indices(self.config.num_queries, self.config.L[0].len(), transcript);

        for i in 0..(num_rounds - 1) {
            reduce_indices(&mut indices, self.config.L[i + 1].len());
            for (j, openings) in proof.openings[i].iter().enumerate() {
                let (alpha_0_opening, alpha_1_opening, beta_opening) = openings;
                let L_i = &self.config.L[i];
                let L_i_next = &self.config.L[i + 1];
                let index = indices[j];
                let s_0 = L_i[index];
                let s_1 = L_i[index + (L_i.len() / 2)];
                let y = L_i_next[index];
                assert_eq!(y, s_0 * s_0);
                assert_eq!(y, s_1 * s_1);

                let alpha_0 = alpha_0_opening.leaf;
                let alpha_1 = alpha_1_opening.leaf;
                let beta = beta_opening.leaf;

                let (coeff, intercept) =
                    if i == 0 {
                        // Interpolate by evaluating the weighted sum
                        let mut weighted_alpha_0 = F::zero();
                        let mut weighted_alpha_1 = F::zero();

                        let max_poly_degree = proof.poly_degrees.iter().max().unwrap();
                        let k = (max_poly_degree + 1).next_power_of_two();

                        for (l, poly_opening) in proof.poly_openings[j].iter().enumerate() {
                            poly_opening.0.verify();
                            poly_opening.1.verify();

                            assert_eq!(poly_opening.0.root, poly_opening.1.root);

                            let alpha_0_x = self.config.L[0][index];
                            let alpha_1_x = self.config.L[0][index + (L_i.len() / 2)];

                            let alpha_0 = (poly_opening.0.leaf - proof.y[l])
                                * (alpha_0_x - proof.z[l]).invert().unwrap();
                            let alpha_1 = (poly_opening.1.leaf - proof.y[l])
                                * (alpha_1_x - proof.z[l]).invert().unwrap();

                            let challenge_weight = &challenge_weights[l];

                            let pad_degree = k - proof.poly_degrees[l] - 1;
                            // Compute w^{alpha_i}^{k-d-1}
                            let alpha_0_x_powered =
                                self.config.L[0][index].pow(&[pad_degree as u64, 0, 0, 0]);
                            // Compute w^{alpha_i}^{k-d-1}
                            let alpha_1_x_powered = self.config.L[0][index + (L_i.len() / 2)]
                                .pow(&[pad_degree as u64, 0, 0, 0]);

                            weighted_alpha_0 += alpha_0 * challenge_weight.0
                                + alpha_0 * alpha_0_x_powered * challenge_weight.1;

                            weighted_alpha_1 += alpha_1 * challenge_weight.0
                                + alpha_1 * alpha_1_x_powered * challenge_weight.1;
                        }

                        // Interpolate the points (s_0, alpha_0) (s_1, alpha_1)
                        let coeff =
                            (weighted_alpha_0 - weighted_alpha_1) * (s_0 - s_1).invert().unwrap();
                        let intercept = weighted_alpha_0 - coeff * s_0;

                        (coeff, intercept)
                    } else {
                        // Interpolate the points (s_0, alpha_0) (s_1, alpha_1)
                        let coeff = (alpha_0 - alpha_1) * (s_0 - s_1).invert().unwrap();
                        let intercept = alpha_0 - coeff * s_0;
                        (coeff, intercept)
                    };

                assert_eq!(
                    beta,
                    coeff * challenges[i] + intercept,
                    "failed at round = {}",
                    i
                );

                assert!(alpha_0_opening.verify());
                assert!(alpha_1_opening.verify());
                assert!(beta_opening.verify());

                // Check that the alpha_0 and alpha_1 come from the same oracle
                assert_eq!(alpha_0_opening.root, alpha_1_opening.root);
            }
        }
    }
}
