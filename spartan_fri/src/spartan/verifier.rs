use crate::spartan::polynomial::eq_poly::EqPoly;
use crate::{FieldExt, MultilinearPCSOpening};

use crate::spartan::indexer::IndexedR1CS;
use crate::spartan::sumcheck::{SumCheckPhase1, SumCheckPhase2};
use crate::spartan::SpartanProof;
use crate::transcript::{AppendToTranscript, Transcript};
use crate::MultilinearPCS;

pub struct SpartanVerifier<F: FieldExt, ML_PCS: MultilinearPCS<F>> {
    index_r1cs: IndexedR1CS<F, ML_PCS>,
    pcs_witness: ML_PCS,
    pcs_index: ML_PCS,
}

impl<F: FieldExt, ML_PCS: MultilinearPCS<F>> SpartanVerifier<F, ML_PCS> {
    pub fn new(index_r1cs: IndexedR1CS<F, ML_PCS>, pcs_witness: ML_PCS, pcs_index: ML_PCS) -> Self {
        Self {
            index_r1cs,
            pcs_witness,
            pcs_index,
        }
    }

    pub fn verify(&self, proof: &SpartanProof<F, ML_PCS>, transcript: &Transcript<F>) -> bool {
        let mut transcript = transcript.clone();
        proof.z_comm.append_to_transcript(&mut transcript);
        let m = (self.index_r1cs.num_vars as f64).log2() as usize;
        let tau = transcript.challenge_vec(m);
        let rx = transcript.challenge_vec(m);
        let mut rx_rev = rx.clone();
        rx_rev.reverse();

        let ex = SumCheckPhase1::verify_round_polys(&proof.sc_proof_1, &rx);

        // The final eval should equal
        let v_A = proof.v_A;
        let v_B = proof.v_B;
        let v_C = proof.v_C;

        let T_1_eq = EqPoly::new(tau);
        let T_1 = (v_A * v_B - v_C) * T_1_eq.eval(&rx_rev);
        assert_eq!(T_1, ex);

        transcript.append_fe(&v_A);
        transcript.append_fe(&v_B);
        transcript.append_fe(&v_C);

        let r = transcript.challenge_vec(3);
        let r_A = r[0];
        let r_B = r[1];
        let r_C = r[2];

        let T_2 = r_A * v_A + r_B * v_B + r_C * v_C;

        let ry = transcript.challenge_vec(m);
        let final_poly_eval = SumCheckPhase2::verify_round_polys(T_2, &proof.sc_proof_2, &ry);

        let mut ry_rev = ry.clone();
        ry_rev.reverse();

        let rx_ry = [ry_rev.clone(), rx_rev].concat();
        assert_eq!(proof.z_eval_proof.x(), ry_rev);
        assert_eq!(proof.A_eval_proof.x(), rx_ry);
        assert_eq!(proof.B_eval_proof.x(), rx_ry);
        assert_eq!(proof.C_eval_proof.x(), rx_ry);

        let z_eval = proof.z_eval_proof.y();
        let A_eval = proof.A_eval_proof.y();
        let B_eval = proof.B_eval_proof.y();
        let C_eval = proof.C_eval_proof.y();

        self.pcs_witness
            .verify(&proof.z_eval_proof, &proof.z_comm, &mut transcript);

        self.pcs_index.verify(
            &proof.A_eval_proof,
            &self.index_r1cs.a_comm,
            &mut transcript,
        );
        self.pcs_index.verify(
            &proof.B_eval_proof,
            &self.index_r1cs.b_comm,
            &mut transcript,
        );
        self.pcs_index.verify(
            &proof.C_eval_proof,
            &self.index_r1cs.c_comm,
            &mut transcript,
        );

        let T_opened = (r_A * A_eval + r_B * B_eval + r_C * C_eval) * z_eval;
        assert_eq!(T_opened, final_poly_eval);

        true
    }
}
