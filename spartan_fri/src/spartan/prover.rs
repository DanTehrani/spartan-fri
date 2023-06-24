#![allow(non_snake_case)]
use std::marker::PhantomData;
use std::time::Instant;

use ark_std::iterable::Iterable;
use ark_std::{end_timer, start_timer};

use crate::spartan::polynomial::eq_poly::EqPoly;
use crate::spartan::polynomial::ml_poly::MlPoly;
use crate::spartan::polynomial::sparse_ml_poly::SparseMLPoly;
use crate::spartan::sumcheck::{SumCheckPhase1, SumCheckPhase2};
use crate::spartan::PartialSpartanProof;
use crate::transcript::{AppendToTranscript, Transcript};
use crate::FieldExt;
use crate::MultilinearPCS;
use crate::R1CS;

use super::utils::boolean_hypercube;

pub struct SpartanProver<F: FieldExt, PCS: MultilinearPCS<F>> {
    r1cs: R1CS<F>,
    pcs_witness: PCS,
}

impl<F: FieldExt, PCS: MultilinearPCS<F>> SpartanProver<F, PCS> {
    pub fn new(r1cs: R1CS<F>, pcs_witness: PCS) -> Self {
        Self { r1cs, pcs_witness }
    }

    pub fn prove(
        &self,
        witness: &[F],
        transcript: &mut Transcript<F>,
    ) -> (PartialSpartanProof<F, PCS>, Vec<F>) {
        // Compute the multilinear extension of the witness
        let mut witness_poly = MlPoly::new(witness.to_vec());
        witness_poly.compute_coeffs();

        // Commit the witness polynomial
        let comm_witness_timer = start_timer!(|| "Commit witness");
        let witness_comm = self.pcs_witness.commit(&witness_poly);
        end_timer!(comm_witness_timer);

        witness_comm.append_to_transcript(transcript);

        // ############################
        // Phase 1: The sum-checks
        // ###################

        let m = (self.r1cs.num_vars as f64).log2() as usize;
        let tau = transcript.challenge_vec(m);

        // First
        // Compute the multilinear extension of the R1CS matrices.
        // Prove that he Q_poly is a zero-polynomial

        // Q_poly is a zero-polynomial iff F_io evaluates to zero
        // over the m-dimensional boolean hypercube..

        // We prove using the sum-check protocol.

        // G_poly = A_poly * B_poly - C_poly

        let num_rows = self.r1cs.num_cons;
        let Az_poly = self.r1cs.A.mul_vector(num_rows, witness);
        let Bz_poly = self.r1cs.B.mul_vector(num_rows, witness);
        let Cz_poly = self.r1cs.C.mul_vector(num_rows, witness);

        // Prove that the polynomial Q(t)
        // \sum_{x \in {0, 1}^m} (Az_poly(x) * Bz_poly(x) - Cz_poly(x)) eq(tau, x)
        // is a zero-polynomial using the sum-check protocol.

        let rx = transcript.challenge_vec(m);
        let mut rx_rev = rx.clone();
        rx_rev.reverse();

        let sc_phase_1_timer = start_timer!(|| "Sumcheck phase 1");

        let sc_phase_1 = SumCheckPhase1::new(
            Az_poly.clone(),
            Bz_poly.clone(),
            Cz_poly.clone(),
            tau.clone(),
            rx.clone(),
        );
        let (sc_proof_1, (v_A, v_B, v_C)) = sc_phase_1.prove(transcript);
        end_timer!(sc_phase_1_timer);

        transcript.append_fe(&v_A);
        transcript.append_fe(&v_B);
        transcript.append_fe(&v_C);

        // Phase 2
        let r = transcript.challenge_vec(3);

        // T_2 should equal teh evaluations of the random linear combined polynomials

        let ry = transcript.challenge_vec(m);
        let sc_phase_2_timer = start_timer!(|| "Sumcheck phase 2");
        let sc_phase_2 = SumCheckPhase2::new(
            self.r1cs.A.clone(),
            self.r1cs.B.clone(),
            self.r1cs.C.clone(),
            witness.to_vec(),
            rx.clone(),
            r.as_slice().try_into().unwrap(),
            ry.clone(),
        );

        let sc_proof_2 = sc_phase_2.prove(transcript);
        end_timer!(sc_phase_2_timer);

        let mut ry_rev = ry.clone();
        ry_rev.reverse();
        let z_open_timer = start_timer!(|| "Open witness poly");
        // Prove the evaluation of the polynomial Z(y) at ry
        let z_eval_proof = self.pcs_witness.open(&witness_poly, &ry_rev, transcript);
        end_timer!(z_open_timer);

        // Prove the evaluation of the polynomials A(y), B(y), C(y) at ry

        let rx_ry = vec![ry_rev, rx_rev].concat();
        (
            PartialSpartanProof {
                z_comm: witness_comm,
                sc_proof_1,
                sc_proof_2,
                z_eval_proof,
                v_A,
                v_B,
                v_C,
            },
            rx_ry,
        )
    }
}
