pub mod polynomial;
pub mod prover;
pub mod sumcheck;
pub mod utils;
pub mod verifier;

use crate::r1cs::R1CS;
use crate::spartan::sumcheck::{SCPhase1Proof, SCPhase2Proof};
use crate::transcript::Transcript;
use crate::FieldExt;
use crate::PolyCommitment;
use std::marker::PhantomData;

// Public parameters for the Spartan IOP
#[derive(Clone)]
pub struct SpartanPP<F: FieldExt, PCS: PolyCommitment<F>> {
    r1cs: R1CS<F>,
    transcript: Transcript<F>,
    pcs: PCS,
}

impl<F: FieldExt, PCS: PolyCommitment<F>> SpartanPP<F, PCS> {
    pub fn new(r1cs: R1CS<F>, label: &'static [u8]) -> Self {
        let transcript = Transcript::<F>::new(label);
        let pcs = PCS::new();
        Self {
            transcript,
            pcs,
            r1cs,
        }
    }
}

pub struct SpartanProof<F: FieldExt, PCS: PolyCommitment<F>> {
    pub sc_proof_1: SCPhase1Proof<F>,
    pub sc_proof_2: SCPhase2Proof<F>,
    _marker: PhantomData<PCS>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{PolyCommitment, SpartanProver, SpartanVerifier, R1CS};
    use ff::Field;
    use pasta_curves::Fp;

    type F = Fp;

    #[derive(Clone)]
    pub struct MockPCS {}
    impl PolyCommitment<F> for MockPCS {
        type Commitment = F;
        type Opening = F;
        fn new() -> Self {
            Self {}
        }

        fn commit(&self, evals: &[F]) -> Self::Commitment {
            F::ZERO
        }

        fn open(&self, point: &[F]) -> Self::Opening {
            F::ZERO
        }
    }

    #[test]
    fn test_spartan() {
        let num_cons = 2usize.pow(5);
        let num_vars = num_cons;
        let num_input = 5;

        let (r1cs, witness) = R1CS::<F>::produce_synthetic_r1cs(num_cons, num_vars, num_input);
        let pp = SpartanPP::<F, MockPCS>::new(r1cs, b"test_spartan");
        let prover = SpartanProver::new(pp.clone());
        let mut proof = prover.prove(&witness);

        let verifier = SpartanVerifier::new(pp);
        assert!(verifier.verify(&proof), "Verification failed");

        // Should assert invalid round-1 sum-check proof
        proof.sc_proof_1.round_polys[0].coeffs[0] =
            proof.sc_proof_1.round_polys[0].coeffs[0] + Fp::one();
        assert!(
            std::panic::catch_unwind(|| verifier.verify(&proof)).is_err(),
            "Verification should fail"
        );
        proof.sc_proof_1.round_polys[0].coeffs[0] =
            proof.sc_proof_1.round_polys[0].coeffs[0] - Fp::one();

        // Should assert invalid round-2 sum-check proof
        proof.sc_proof_2.round_polys[0].coeffs[0] =
            proof.sc_proof_2.round_polys[0].coeffs[0] + Fp::one();
        assert!(
            std::panic::catch_unwind(|| verifier.verify(&proof)).is_err(),
            "Verification should fail"
        );
    }
}
