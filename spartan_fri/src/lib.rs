#![allow(non_snake_case)]
mod eq_poly;
mod fri;
mod ml_poly;
mod prover;
mod r1cs;
mod sc_phase_1;
mod sc_phase_2;
mod transcript;
mod utils;
mod verifier;

use pasta_curves::arithmetic::FieldExt;
use sc_phase_1::SCPhase1Proof;
use sc_phase_2::SCPhase2Proof;

pub use prover::SpartanFRIProver;
pub use r1cs::R1CS;
pub use transcript::Transcript;

pub struct SpartanFRIProof<F: FieldExt> {
    pub sc_proof_1: SCPhase1Proof<F>,
    pub sc_proof_2: SCPhase2Proof<F>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{prover::SpartanFRIProver, r1cs::R1CS, verifier::SpartanFRIVerifier};

    use pasta_curves::Fp;
    use transcript::Transcript;
    type F = Fp;

    #[test]
    fn test_spartan_fri() {
        let num_cons = 2usize.pow(4);
        let num_vars = num_cons;
        let num_input = 5;

        let r1cs = R1CS::<F>::produce_synthetic_r1cs(num_cons, num_vars, num_input);
        let prover_transcript = Transcript::<F>::new(b"test_spartan_fri");

        let mut prover = SpartanFRIProver::<F>::new(r1cs.clone(), prover_transcript);
        let proof = prover.prove();

        let verifier_transcript = Transcript::<F>::new(b"test_spartan_fri");
        let mut verifier = SpartanFRIVerifier::<F>::new(r1cs.clone(), verifier_transcript);
        assert!(verifier.verify(&proof));
    }
}
