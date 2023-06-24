use ark_std::{end_timer, start_timer};

use crate::{fri::MLPolyEvalProof, transcript::Transcript, FieldExt, MultilinearPCS, R1CS};

use super::{FullSpartanProof, PartialSpartanProof};

pub struct Assistant<F: FieldExt, PCS: MultilinearPCS<F>> {
    pub r1cs: R1CS<F>,
    pub pcs: PCS,
}

impl<F: FieldExt, PCS: MultilinearPCS<F>> Assistant<F, PCS> {
    pub fn new(r1cs: R1CS<F>, pcs: PCS) -> Self {
        Self { r1cs, pcs }
    }

    pub fn complete_proof(
        &self,
        partial_proof: PartialSpartanProof<F, PCS>,
        x: &[F],
        transcript: &mut Transcript<F>,
    ) -> FullSpartanProof<F, PCS> {
        let compute_ml_timer = start_timer!(|| "Compute MLE");
        // This part can be ported to the remote server
        let s = (self.r1cs.num_vars as f64).log2() as usize;
        let A_mle = self.r1cs.A.to_ml_extension(s);
        let B_mle = self.r1cs.B.to_ml_extension(s);
        let C_mle = self.r1cs.C.to_ml_extension(s);
        end_timer!(compute_ml_timer);

        let open_matrices_timer = start_timer!(|| "Open matrices");

        let A_eval_proof = self.pcs.open(&A_mle, x, transcript);
        let B_eval_proof = self.pcs.open(&B_mle, x, transcript);
        let C_eval_proof = self.pcs.open(&C_mle, x, transcript);
        end_timer!(open_matrices_timer);

        FullSpartanProof {
            partial_proof,
            A_eval_proof,
            B_eval_proof,
            C_eval_proof,
        }
    }
}
