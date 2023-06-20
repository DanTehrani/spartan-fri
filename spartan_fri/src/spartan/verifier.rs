use crate::spartan::polynomial::eq_poly::EqPoly;
use crate::{FieldExt, MultilinearPCSOpening};

use crate::spartan::indexer::IndexedR1CS;
use crate::spartan::sumcheck::{SumCheckPhase1, SumCheckPhase2};
use crate::spartan::SpartanProof;
use crate::transcript::{AppendToTranscript, Transcript};
use crate::MultilinearPCS;

pub struct SpartanVerifier<F: FieldExt, ML_PCS: MultilinearPCS<F>> {
    index_r1cs: IndexedR1CS<F, ML_PCS>,
    pcs: ML_PCS,
}

impl<F: FieldExt, ML_PCS: MultilinearPCS<F>> SpartanVerifier<F, ML_PCS> {
    pub fn new(index_r1cs: IndexedR1CS<F, ML_PCS>, pcs: ML_PCS) -> Self {
        Self { index_r1cs, pcs }
    }

    pub fn verify(&self, proof: &SpartanProof<F, ML_PCS>, transcript: &Transcript<F>) -> bool {
        let mut transcript = transcript.clone();
        proof.z_comm.append_to_transcript(&mut transcript);
        let m = (self.index_r1cs.num_vars as f64).log2() as usize;
        let tau = transcript.challenge_vec(m);
        let rx = transcript.challenge_vec(m);

        let (ex, sum) = SumCheckPhase1::verify_round_polys(&proof.sc_proof_1, &rx);
        // TODO: The sum should be zero.

        // The final eval should equal
        let v_A = proof.v_A;
        let v_B = proof.v_B;
        let v_C = proof.v_C;

        let T_1_eq = EqPoly::new(tau);
        let T_1 = (v_A + v_B - v_C) * T_1_eq.eval(&rx);
        assert_eq!(T_1, ex);

        transcript.append_fe(&v_A);
        transcript.append_fe(&v_B);
        transcript.append_fe(&v_C);

        let r = transcript.challenge_vec(3);
        let r_A = r[0];
        let r_B = r[1];
        let r_C = r[2];

        let ry = transcript.challenge_vec(m);
        let (sum_claim, final_poly_eval) =
            SumCheckPhase2::verify_round_polys(&proof.sc_proof_2, &ry);

        assert_eq!(proof.z_eval_proof.x(), ry);

        let z_eval = proof.z_eval_proof.y();
        let A_eval = proof.A_eval_proof.y();
        let B_eval = proof.B_eval_proof.y();
        let C_eval = proof.C_eval_proof.y();

        let T_opened = (r_A * A_eval + r_B * B_eval + r_C * C_eval) * z_eval;
        assert_eq!(T_opened, final_poly_eval);

        true
    }
}
