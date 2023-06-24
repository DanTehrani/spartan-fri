#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
use ff::{FromUniformBytes, PrimeField};
use fri::FRIMLPolyCommit;
//use fri::{FRIMLPolyCommit, UniPoly};
use pasta_curves::Fp;
use serde::Serialize;
use spartan::polynomial::{ml_poly::MlPoly, sparse_ml_poly::SparseMLPoly};
use spartan::prover::SpartanProver;
use spartan::verifier::SpartanVerifier;
use spartan::PartialSpartanProof;

pub mod fri;
mod r1cs;
pub mod spartan;
pub mod transcript;

pub trait FieldExt: FromUniformBytes<64, Repr = [u8; 32]> + Serialize {}

impl FieldExt for Fp {}

pub trait MultilinearPCS<F: FieldExt>: Clone {
    type Commitment: AppendToTranscript<F>;
    type Opening: MultilinearPCSOpening<F>;
    type Config;
    fn new(config: Self::Config) -> Self;
    fn commit(&self, poly: &MlPoly<F>) -> Self::Commitment;
    fn open(&self, poly: &MlPoly<F>, point: &[F], transcript: &mut Transcript<F>) -> Self::Opening;
    fn verify(
        &self,
        opening: &Self::Opening,
        commitment: &Self::Commitment,
        transcript: &mut Transcript<F>,
    );
}
// Commitment scheme for multilinear polynomials
pub trait SparseMultilinearPCS<F: FieldExt>: Clone {
    type Commitment: AppendToTranscript<F>;
    type Opening: MultilinearPCSOpening<F>;
    type Config;
    fn new(config: Self::Config) -> Self;
    fn commit(&self, poly: &SparseMLPoly<F>) -> Self::Commitment;
    fn open(
        &self,
        poly: &SparseMLPoly<F>,
        point: &[F],
        transcript: &mut Transcript<F>,
    ) -> Self::Opening;
    fn verify(
        &self,
        opening: &Self::Opening,
        commitment: &Self::Commitment,
        transcript: &mut Transcript<F>,
    );
}

pub trait MultilinearPCSOpening<F: FieldExt>: Clone {
    fn x(&self) -> Vec<F>;
    fn y(&self) -> F;
}

// #######################
// Re-exports
// #######################

pub type SpartanFRIProver<F> = SpartanProver<F, FRIMLPolyCommit<F>>;
pub type SpartanFRIProof<F> = PartialSpartanProof<F, FRIMLPolyCommit<F>>;
pub use r1cs::R1CS;
use transcript::{AppendToTranscript, Transcript};

#[cfg(test)]
mod tests {
    use crate::fri::FRIConfig;
    use crate::fri::FRIMLPolyCommit;
    use crate::spartan::assistant::Assistant;
    use crate::spartan::indexer::Indexer;
    use crate::transcript::Transcript;
    use crate::{MultilinearPCS, SpartanProver, SpartanVerifier, R1CS};
    use ark_std::end_timer;
    use ark_std::start_timer;
    use pasta_curves::Fp;

    type F = Fp;

    #[test]
    fn test_spartan() {
        let num_cons = 2usize.pow(4);
        let num_vars = num_cons;
        let num_input = 0;

        let (r1cs, witness) = R1CS::<F>::produce_synthetic_r1cs(num_cons, num_vars, num_input);

        let witness_poly_degree = num_cons;
        let num_queries = 30;
        let expansion_factor = 2;
        let folding_factor = 2;
        let final_codeword_size = 1;

        let fri_config = FRIConfig::<F>::new(
            witness_poly_degree,
            expansion_factor,
            folding_factor,
            num_queries,
            final_codeword_size,
        );
        let pcs_witness = FRIMLPolyCommit::<F>::new(fri_config);

        // #####################
        // Prover
        // #####################

        let proving_timer = start_timer!(|| "Proving");
        let prover = SpartanProver::<F, FRIMLPolyCommit<F>>::new(r1cs.clone(), pcs_witness.clone());
        let mut prover_transcript = Transcript::<F>::new(b"test_spartan");
        let (partial_proof, rx_ry) = prover.prove(&witness, &mut prover_transcript);
        end_timer!(proving_timer);

        // #####################
        // Assistant
        // #####################
        let indexer_fri_config = FRIConfig::<F>::new(
            num_cons * num_vars,
            expansion_factor,
            folding_factor,
            num_queries,
            final_codeword_size,
        );
        let pcs_indexer = FRIMLPolyCommit::<F>::new(indexer_fri_config);

        // Pre-process the R1CS
        let indexer = Indexer::new(r1cs.clone(), pcs_indexer.clone());
        let index_r1cs = indexer.pre_process();

        let assistant = Assistant::new(r1cs.clone(), pcs_indexer.clone());
        let full_proof = assistant.complete_proof(partial_proof, &rx_ry, &mut prover_transcript);
        let verifier = SpartanVerifier::new(index_r1cs, pcs_witness.clone(), pcs_indexer.clone());
        let verifier_transcript = Transcript::<F>::new(b"test_spartan");
        assert!(
            verifier.verify(&full_proof, &verifier_transcript),
            "Verification failed"
        );

        /*
        // Should assert invalid round-1 sum-check proof
        proof.sc_proof_1.round_polys[0].coeffs[0] =
            proof.sc_proof_1.round_polys[0].coeffs[0] + Fp::one();
        assert!(
            std::panic::catch_unwind(|| verifier.verify(&proof, &verifier_transcript)).is_err(),
            "Verification should fail"
        );
        proof.sc_proof_1.round_polys[0].coeffs[0] =
            proof.sc_proof_1.round_polys[0].coeffs[0] - Fp::one();

        // Should assert invalid round-2 sum-check proof
        proof.sc_proof_2.round_polys[0].coeffs[0] =
            proof.sc_proof_2.round_polys[0].coeffs[0] + Fp::one();
        assert!(
            std::panic::catch_unwind(|| verifier.verify(&proof, &verifier_transcript)).is_err(),
            "Verification should fail"
        );
         */
    }
}
