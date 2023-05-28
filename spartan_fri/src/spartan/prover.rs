#![allow(non_snake_case)]
use std::marker::PhantomData;

use crate::spartan::polynomial::ml_poly::MlPoly;
use crate::spartan::sumcheck::{SumCheckPhase1, SumCheckPhase2};
use crate::spartan::{SpartanPP, SpartanProof};
use crate::FieldExt;
use crate::PolyCommitment;
use ff::Field;

pub struct SpartanProver<F: FieldExt, PCS: PolyCommitment<F>> {
    pp: SpartanPP<F, PCS>,
    _marker: PhantomData<F>,
}

impl<F: FieldExt, PCS: PolyCommitment<F>> SpartanProver<F, PCS> {
    pub fn new(pp: SpartanPP<F, PCS>) -> Self {
        Self {
            pp,
            _marker: PhantomData,
        }
    }

    pub fn prove(&self, witness: &[F]) -> SpartanProof<F, PCS> {
        let mut transcript = self.pp.transcript.clone();
        // Commit the witness polynomial

        // Derive the commitment to the quotient polynomial
        // Q(ry) = w(ry) - v / (X - ry)

        // How to derive the commitment to the quotient polynomial?
        // 1. Just query at some point x_i in the domain.
        // 2. Compute (f(x_i) - y) / (x_i - ry)

        let m = (self.pp.r1cs.num_vars as f64).log2() as usize;
        // First
        // Compute the multilinear extension of the R1CS matrices.
        // Prove that he Q_poly is a zero-polynomial

        // Q_poly is a zero-polynomial iff F_io evaluates to zero
        // over the m-dimensional boolean hypercube..

        // We prove using the sum-check protocol.

        // G_poly = A_poly * B_poly - C_poly

        // Represent the terms as polynomiales multiplied
        let Az = self.pp.r1cs.A.mul_vector(self.pp.r1cs.num_cons, witness);
        let Bz = self.pp.r1cs.B.mul_vector(self.pp.r1cs.num_cons, witness);
        let Cz = self.pp.r1cs.C.mul_vector(self.pp.r1cs.num_cons, witness);

        let Az_poly = MlPoly::new(Az);
        let Bz_poly = MlPoly::new(Bz);
        let Cz_poly = MlPoly::new(Cz);

        // Prove that the polynomial Q(t)
        // \sum_{x \in {0, 1}^m} (Az_poly(x) * Bz_poly(x) - Cz_poly(x)) eq(tau, x)
        // is a zero-polynomial using the sum-check protocol.

        let rx = transcript.challenge_vec(m);

        let sc_phase_1 = SumCheckPhase1::new(Az_poly.clone(), Bz_poly.clone(), Cz_poly.clone(), rx);
        let sc_proof_1 = sc_phase_1.prove();

        let r = transcript.challenge_vec(3);
        let r_A = r[0];
        let r_B = r[1];
        let r_C = r[2];

        let ry = transcript.challenge_vec(m);

        let v_A = Az_poly.eval(&ry);
        let v_B = Bz_poly.eval(&ry);
        let v_C = Cz_poly.eval(&ry);

        let T_2 = r_A * v_A + r_B * v_B + r_C * v_C;

        let sc_phase_2 = SumCheckPhase2::new(
            Az_poly,
            Bz_poly,
            Cz_poly,
            r.as_slice().try_into().unwrap(),
            ry,
        );

        let sc_proof_2 = sc_phase_2.prove();

        SpartanProof {
            sc_proof_1,
            sc_proof_2,
            _marker: PhantomData,
        }
    }
}
