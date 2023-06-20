use crate::{fri::UniPoly, FieldExt, MultilinearPCS, R1CS};

pub struct Indexer<F: FieldExt, PCS: MultilinearPCS<F>> {
    pub r1cs: R1CS<F>,
    pub pcs: PCS,
}

pub struct IndexedR1CS<F: FieldExt, PCS: MultilinearPCS<F>> {
    pub a_comm: PCS::Commitment,
    pub b_comm: PCS::Commitment,
    pub c_comm: PCS::Commitment,
    pub num_cons: usize,
    pub num_vars: usize,
}

impl<F: FieldExt, PCS: MultilinearPCS<F>> Indexer<F, PCS> {
    pub fn new(r1cs: R1CS<F>, pcs: PCS) -> Self {
        Self { r1cs, pcs }
    }

    pub fn pre_process(&self) -> IndexedR1CS<F, PCS> {
        // Compute the multilinear extension of the R1CS matrices.
        // Commit the evaluations
        let s = (self.r1cs.num_vars as f64).log2() as usize;
        let a_mle = self.r1cs.A.to_ml_extension(s);
        let b_mle = self.r1cs.B.to_ml_extension(s);
        let c_mle = self.r1cs.C.to_ml_extension(s);

        let a_comm = self.pcs.commit(&a_mle);
        let b_comm = self.pcs.commit(&b_mle);
        let c_comm = self.pcs.commit(&c_mle);

        IndexedR1CS {
            a_comm,
            b_comm,
            c_comm,
            num_cons: self.r1cs.num_cons,
            num_vars: self.r1cs.num_vars,
        }
    }
}
