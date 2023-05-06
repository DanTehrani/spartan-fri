use pasta_curves::arithmetic::FieldExt;

use crate::r1cs::R1CS;
use crate::sc_phase_1::SumCheckPhase1;
use crate::sc_phase_2::SumCheckPhase2;
use crate::transcript::Transcript;
use crate::SpartanFRIProof;

pub struct SpartanFRIVerifier<F: FieldExt> {
    r1cs: R1CS<F>,
    transcript: Transcript<F>,
}

impl<F: FieldExt<Repr = [u8; 32]>> SpartanFRIVerifier<F> {
    pub fn new(r1cs: R1CS<F>, transcript: Transcript<F>) -> Self {
        Self { r1cs, transcript }
    }

    pub fn verify(&mut self, proof: &SpartanFRIProof<F>) -> bool {
        let m = (self.r1cs.num_vars as f64).log2() as usize;
        let rx = self.transcript.challenge_vec(m);

        let (eval_claim, sum) = SumCheckPhase1::verify_round_polys(&proof.sc_proof_1, &rx);

        let r = self.transcript.challenge_vec(3);
        let r_A = r[0];
        let r_B = r[1];
        let r_C = r[2];

        let ry = self.transcript.challenge_vec(m);
        let (sum_claim, final_poly_eval) =
            SumCheckPhase2::verify_round_polys(&proof.sc_proof_2, &ry);

        true
    }
}
