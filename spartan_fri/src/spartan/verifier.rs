use crate::FieldExt;
use ff::Field;

use crate::spartan::sumcheck::{SumCheckPhase1, SumCheckPhase2};
use crate::spartan::{SpartanPP, SpartanProof};
use crate::PolyCommitment;

pub struct SpartanVerifier<F: FieldExt, PCS: PolyCommitment<F>> {
    pp: SpartanPP<F, PCS>,
}

impl<F: FieldExt, PCS: PolyCommitment<F>> SpartanVerifier<F, PCS> {
    pub fn new(pp: SpartanPP<F, PCS>) -> Self {
        Self { pp }
    }

    pub fn verify(&self, proof: &SpartanProof<F, PCS>) -> bool {
        let mut transcript = self.pp.transcript.clone();
        let m = (self.pp.r1cs.num_vars as f64).log2() as usize;
        let rx = transcript.challenge_vec(m);

        let (eval_claim, sum) = SumCheckPhase1::verify_round_polys(&proof.sc_proof_1, &rx);

        let r = transcript.challenge_vec(3);
        let r_A = r[0];
        let r_B = r[1];
        let r_C = r[2];

        let ry = transcript.challenge_vec(m);
        let (sum_claim, final_poly_eval) =
            SumCheckPhase2::verify_round_polys(&proof.sc_proof_2, &ry);

        true
    }
}
