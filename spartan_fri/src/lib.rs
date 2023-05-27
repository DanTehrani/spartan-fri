#![allow(non_snake_case)]
use fri::FRIMLPolyCommitProver;
use pasta_curves::arithmetic::FieldExt;
use spartan::prover::SpartanProver;
use spartan::verifier::SpartanVerifier;

mod fri;
mod r1cs;
mod spartan;
mod transcript;

use spartan::SpartanProof;

// Commitment scheme for multilinear polynomials
pub trait PolyCommitment<F: FieldExt>: Clone {
    type Commitment;
    type Opening;
    fn new() -> Self;
    fn commit(&self, evals: &[F]) -> Self::Commitment;
    fn open(&self, point: &[F]) -> Self::Opening;
}

// #######################
// Re-exports
// #######################

pub type SpartanFRIProver<F> = SpartanProver<F, FRIMLPolyCommitProver<F>>;
pub type SpartanFRIProof<F> = SpartanProof<F, FRIMLPolyCommitProver<F>>;
pub use r1cs::R1CS;
pub use spartan::SpartanPP;

/*
#[cfg(test)]
mod tests {
    use super::*;

    use r1cs::R1CS;
    use pasta_curves::Fp;
    use transcript::Transcript;
    use spartan::prover

    type F = Fp;

    #[test]
    fn test_spartan_fri() {
        let num_cons = 2usize.pow(4);
        let num_vars = num_cons;
        let num_input = 5;

        let r1cs = R1CS::<F>::produce_synthetic_r1cs(num_cons, num_vars, num_input);
        let prover_transcript = Transcript::<F>::new(b"test_spartan_fri");

        let mut prover = SpartanProver::<F>::new(r1cs.clone(), prover_transcript);
        let proof = prover.prove();

        let verifier_transcript = Transcript::<F>::new(b"test_spartan_fri");
        let mut verifier = SpartanVerifier::<F>::new(r1cs.clone(), verifier_transcript);
        assert!(verifier.verify(&proof));
    }
}
 */
