use crate::spartan::polynomial::eq_poly::EqPoly;
use crate::{FieldExt, MultilinearPCSOpening};

use crate::spartan::indexer::IndexedR1CS;
use crate::spartan::sumcheck::{SumCheckPhase1, SumCheckPhase2};
use crate::spartan::PartialSpartanProof;
use crate::transcript::{AppendToTranscript, Transcript};
use crate::MultilinearPCS;

use super::FullSpartanProof;

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

    pub fn verify(&self, proof: &FullSpartanProof<F, ML_PCS>, transcript: &Transcript<F>) -> bool {
        let mut transcript = transcript.clone();
        let partial_proof = &proof.partial_proof;
        partial_proof.z_comm.append_to_transcript(&mut transcript);
        let m = (self.index_r1cs.num_vars as f64).log2() as usize;
        let tau = transcript.challenge_vec(m);
        let rx = transcript.challenge_vec(m);
        let mut rx_rev = rx.clone();
        rx_rev.reverse();

        transcript.append_fe(&partial_proof.sc_proof_1.blinder_poly_sum);
        let rho = transcript.challenge_fe();

        let ex = SumCheckPhase1::verify_round_polys(&partial_proof.sc_proof_1, &rx, rho);

        // The final eval should equal
        let v_A = partial_proof.v_A;
        let v_B = partial_proof.v_B;
        let v_C = partial_proof.v_C;

        let T_1_eq = EqPoly::new(tau);
        let T_1 = (v_A * v_B - v_C) * T_1_eq.eval(&rx_rev)
            + rho * partial_proof.sc_proof_1.blinder_poly_eval_claim;
        assert_eq!(T_1, ex);

        transcript.append_fe(&v_A);
        transcript.append_fe(&v_B);
        transcript.append_fe(&v_C);

        let r = transcript.challenge_vec(3);
        let r_A = r[0];
        let r_B = r[1];
        let r_C = r[2];

        let ry = transcript.challenge_vec(m);

        transcript.append_fe(&partial_proof.sc_proof_2.blinder_poly_sum);
        let rho_2 = transcript.challenge_fe();

        let T_2 =
            (r_A * v_A + r_B * v_B + r_C * v_C) + rho_2 * partial_proof.sc_proof_2.blinder_poly_sum;
        let final_poly_eval =
            SumCheckPhase2::verify_round_polys(T_2, &partial_proof.sc_proof_2, &ry);

        let mut ry_rev = ry.clone();
        ry_rev.reverse();

        let rx_ry = [ry_rev.clone(), rx_rev].concat();
        assert_eq!(partial_proof.z_eval_proof.x(), ry_rev);
        assert_eq!(proof.A_eval_proof.x(), rx_ry);
        assert_eq!(proof.B_eval_proof.x(), rx_ry);
        assert_eq!(proof.C_eval_proof.x(), rx_ry);

        let z_eval = partial_proof.z_eval_proof.y();
        let A_eval = proof.A_eval_proof.y();
        let B_eval = proof.B_eval_proof.y();
        let C_eval = proof.C_eval_proof.y();

        self.pcs_witness.verify(
            &partial_proof.z_eval_proof,
            &partial_proof.z_comm,
            &mut transcript,
        );

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

        let T_opened = (r_A * A_eval + r_B * B_eval + r_C * C_eval) * z_eval
            + rho_2 * partial_proof.sc_proof_2.blinder_poly_eval_claim;
        assert_eq!(T_opened, final_poly_eval);

        true
    }
}
