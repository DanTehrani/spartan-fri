pub mod indexer;
pub mod polynomial;
pub mod prover;
pub mod sumcheck;
pub mod utils;
pub mod verifier;

use crate::spartan::sumcheck::{SCPhase1Proof, SCPhase2Proof};
use crate::FieldExt;
use crate::MultilinearPCS;
use std::marker::PhantomData;

pub struct SpartanProof<F: FieldExt, PCS: MultilinearPCS<F>> {
    pub z_comm: PCS::Commitment,
    pub sc_proof_1: SCPhase1Proof<F>,
    pub sc_proof_2: SCPhase2Proof<F>,
    pub z_eval_proof: PCS::Opening,
    pub v_A: F,
    pub v_B: F,
    pub v_C: F,
    pub A_eval_proof: PCS::Opening,
    pub B_eval_proof: PCS::Opening,
    pub C_eval_proof: PCS::Opening,
    _marker: PhantomData<PCS>,
}
