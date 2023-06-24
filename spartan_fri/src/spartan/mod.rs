pub mod assistant;
pub mod indexer;
pub mod polynomial;
pub mod prover;
pub mod sumcheck;
pub mod utils;
pub mod verifier;

use ff::Field;

use crate::spartan::sumcheck::{SCPhase1Proof, SCPhase2Proof};
use crate::FieldExt;
use crate::MultilinearPCS;
use std::marker::PhantomData;

pub struct PartialSpartanProof<F: FieldExt, PCS: MultilinearPCS<F>> {
    pub z_comm: PCS::Commitment,
    pub sc_proof_1: SCPhase1Proof<F>,
    pub sc_proof_2: SCPhase2Proof<F>,
    pub z_eval_proof: PCS::Opening,
    pub v_A: F,
    pub v_B: F,
    pub v_C: F,
}

pub struct FullSpartanProof<F: FieldExt, PCS: MultilinearPCS<F>> {
    pub partial_proof: PartialSpartanProof<F, PCS>,
    pub A_eval_proof: PCS::Opening,
    pub B_eval_proof: PCS::Opening,
    pub C_eval_proof: PCS::Opening,
}
