pub mod polynomial;
pub mod prover;
pub mod sumcheck;
pub mod utils;
pub mod verifier;

use crate::r1cs::R1CS;
use crate::spartan::sumcheck::{SCPhase1Proof, SCPhase2Proof};
use crate::transcript::Transcript;
use crate::PolyCommitment;
use pasta_curves::arithmetic::FieldExt;
use std::marker::PhantomData;

// Public parameters for the Spartan IOP
pub struct SpartanPP<F: FieldExt, PCS: PolyCommitment<F>> {
    r1cs: R1CS<F>,
    transcript: Transcript<F>,
    pcs: PCS,
}

impl<F: FieldExt<Repr = [u8; 32]>, PCS: PolyCommitment<F>> SpartanPP<F, PCS> {
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
