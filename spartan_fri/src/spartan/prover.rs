#![allow(non_snake_case)]
use std::marker::PhantomData;

use crate::spartan::polynomial::ml_poly::MlPoly;
use crate::spartan::sumcheck::{SumCheckPhase1, SumCheckPhase2};
use crate::spartan::SpartanProof;
use crate::transcript::{AppendToTranscript, Transcript};
use crate::FieldExt;
use crate::MultilinearPCS;
use crate::R1CS;

pub struct SpartanProver<F: FieldExt, PCS: MultilinearPCS<F>> {
    r1cs: R1CS<F>,
    pcs_witness: PCS,
    pcs_index: PCS,
}

impl<F: FieldExt, PCS: MultilinearPCS<F>> SpartanProver<F, PCS> {
    pub fn new(r1cs: R1CS<F>, pcs_witness: PCS, pcs_index: PCS) -> Self {
        Self {
            r1cs,
            pcs_witness,
            pcs_index,
        }
    }

    pub fn prove(&self, witness: &[F], transcript: &Transcript<F>) -> SpartanProof<F, PCS> {
        let mut transcript: Transcript<F> = transcript.clone();

        // Compute the multilinear extension of the witness
        let witness_poly = MlPoly::new(witness.to_vec());

        // Commit the witness polynomial
        let witness_comm = self.pcs_witness.commit(&witness_poly);
        witness_comm.append_to_transcript(&mut transcript);

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

        // Represent the terms as polynomiales multiplied
        let Az = self.r1cs.A.mul_vector(self.r1cs.num_cons, witness);
        let Bz = self.r1cs.B.mul_vector(self.r1cs.num_cons, witness);
        let Cz = self.r1cs.C.mul_vector(self.r1cs.num_cons, witness);

        let Az_poly = MlPoly::new(Az);
        let Bz_poly = MlPoly::new(Bz);
        let Cz_poly = MlPoly::new(Cz);

        // Prove that the polynomial Q(t)
        // \sum_{x \in {0, 1}^m} (Az_poly(x) * Bz_poly(x) - Cz_poly(x)) eq(tau, x)
        // is a zero-polynomial using the sum-check protocol.

        let rx = transcript.challenge_vec(m);

        let sc_phase_1 = SumCheckPhase1::new(
            Az_poly.clone(),
            Bz_poly.clone(),
            Cz_poly.clone(),
            tau.clone(),
            rx.clone(),
        );
        let sc_proof_1 = sc_phase_1.prove();

        let v_A = Az_poly.eval(&rx);
        let v_B = Bz_poly.eval(&rx);
        let v_C = Cz_poly.eval(&rx);

        transcript.append_fe(&v_A);
        transcript.append_fe(&v_B);
        transcript.append_fe(&v_C);

        // Phase 2
        let r = transcript.challenge_vec(3);
        // let T_2 = r_A * v_A + r_B * v_B + r_C * v_C;

        let r_A = r[0];
        let r_B = r[1];
        let r_C = r[2];

        let ry = transcript.challenge_vec(m);
        let sc_phase_2 = SumCheckPhase2::new(
            Az_poly,
            Bz_poly,
            Cz_poly,
            r.as_slice().try_into().unwrap(),
            ry.clone(),
        );

        let sc_proof_2 = sc_phase_2.prove();

        // Prove the evaluation of the polynomial Z(y) at ry
        let z_eval_proof = self.pcs_witness.open(&witness_poly, &ry, &mut transcript);

        // Prove the evaluation of the polynomials A(y), B(y), C(y) at ry
        // TODO: Pre-compute the MLEs

        // This part can be ported to the remote server
        let s = (self.r1cs.num_vars as f64).log2() as usize;
        let A_mle = self.r1cs.A.to_ml_extension(s);
        let B_mle = self.r1cs.B.to_ml_extension(s);
        let C_mle = self.r1cs.C.to_ml_extension(s);

        let rx_ry = vec![rx, ry].concat();
        let A_eval_proof = self.pcs_index.open(&A_mle, &rx_ry, &mut transcript);
        let B_eval_proof = self.pcs_index.open(&B_mle, &rx_ry, &mut transcript);
        let C_eval_proof = self.pcs_index.open(&C_mle, &rx_ry, &mut transcript);

        SpartanProof {
            z_comm: witness_comm,
            sc_proof_1,
            sc_proof_2,
            z_eval_proof,
            v_A,
            v_B,
            v_C,
            A_eval_proof,
            B_eval_proof,
            C_eval_proof,
            _marker: PhantomData,
        }
    }
}
