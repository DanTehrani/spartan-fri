// Phase 2 sum-check of Spartan.

use crate::spartan::polynomial::eq_poly::EqPoly;
use crate::spartan::polynomial::ml_poly::MlPoly;
use crate::spartan::sumcheck::unipoly::UniPoly;
use crate::spartan::utils::boolean_hypercube;

use crate::FieldExt;

pub struct SCPhase1Proof<F: FieldExt> {
    pub round_polys: Vec<UniPoly<F>>,
}

pub struct SumCheckPhase1<F: FieldExt> {
    A_mle: MlPoly<F>,
    B_mle: MlPoly<F>,
    C_mle: MlPoly<F>,
    Z: Vec<F>,
    bound_eq_poly: EqPoly<F>,
    challenge: Vec<F>,
}

impl<F: FieldExt> SumCheckPhase1<F> {
    pub fn new(
        A_mle: MlPoly<F>,
        B_mle: MlPoly<F>,
        C_mle: MlPoly<F>,
        Z: Vec<F>,
        tau: Vec<F>,
        challenge: Vec<F>,
    ) -> Self {
        let bound_eq_poly = EqPoly::new(tau);
        Self {
            A_mle,
            B_mle,
            C_mle,
            Z,
            bound_eq_poly,
            challenge,
        }
    }

    fn eval_poly(&self, vars: &[F]) -> F {
        let s = self.A_mle.num_vars / 2;
        debug_assert_eq!(vars.len(), s);

        let mut a_eval = F::ZERO;
        let mut b_eval = F::ZERO;
        let mut c_eval = F::ZERO;

        // Create an eq_poly evaluation (for the mle) table for the first s variables.
        // Run a dot product for the last s - j variables

        for (i, b) in boolean_hypercube(s).iter().enumerate() {
            let z_eval = self.Z[i];
            let eval_at = [b.as_slice(), vars].concat();
            // Precompute the eq evaluation table and run dot products?
            a_eval += self.A_mle.eval(&eval_at) * z_eval;
            b_eval += self.B_mle.eval(&eval_at) * z_eval;
            c_eval += self.C_mle.eval(&eval_at) * z_eval;
        }

        (a_eval * b_eval - c_eval) * self.bound_eq_poly.eval(vars)
    }

    pub fn round(&self, j: usize) -> UniPoly<F> {
        // evaluate at points 0, 1, 2, 3

        let zero = F::ZERO;
        let one = F::ONE;

        let m = self.A_mle.num_vars / 2;
        let mut evals = [F::ZERO; 4];

        for (i, vars) in boolean_hypercube(m - j - 1).iter().enumerate() {
            let mut eval_at = vec![];

            // Eval at 0
            eval_at.extend_from_slice(&self.challenge[..j]);
            eval_at.push(zero);
            eval_at.extend_from_slice(vars);

            evals[0] += self.eval_poly(&eval_at);

            // Eval at 1
            eval_at[j] = one;
            evals[1] += self.eval_poly(&eval_at);

            // Eval at 2
            eval_at[j] = F::from(2u64);
            evals[2] += self.eval_poly(&eval_at);

            // Eval at 3
            eval_at[j] = F::from(3u64);
            evals[3] += self.eval_poly(&eval_at);
        }

        let round_poly = UniPoly::interpolate(&evals);
        round_poly
    }

    pub fn prove(&self) -> SCPhase1Proof<F> {
        let num_vars = self.A_mle.num_vars / 2;
        let mut round_polys = Vec::<UniPoly<F>>::with_capacity(num_vars - 1);

        for i in 0..(num_vars - 1) {
            let round_poly = self.round(i);
            round_polys.push(round_poly);
        }

        SCPhase1Proof { round_polys }
    }

    pub fn verify_round_polys(proof: &SCPhase1Proof<F>, challenge: &[F]) -> F {
        debug_assert_eq!(proof.round_polys.len(), challenge.len() - 1);
        let m = challenge.len();

        let zero = F::ZERO;
        let one = F::ONE;

        let round_polys = &proof.round_polys;

        let mut target = zero;
        for (i, round_poly) in proof.round_polys.iter().enumerate() {
            assert_eq!(round_poly.eval(zero) + round_poly.eval(one), target);
            target = round_poly.eval(challenge[i]);
        }

        let final_poly_eval = round_polys[m - 2].eval(challenge[m - 1]);
        final_poly_eval
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pasta_curves::Fp;
    type F = Fp;
    use ff::Field;

    #[test]
    fn test_unipoly_3() {
        let coeffs = [F::from(1u64), F::from(2u64), F::from(3u64), F::from(4u64)];
        let eval_at = Fp::from(33);

        let mut expected_eval = F::ZERO;
        for i in 0..coeffs.len() {
            expected_eval += coeffs[i] * eval_at.pow(&[3 - i as u64, 0, 0, 0]);
        }

        let mut evals = [F::ZERO; 4];
        for i in 0..4 {
            let eval_at = F::from(i as u64);
            let mut eval_i = F::ZERO;
            for j in 0..coeffs.len() {
                eval_i += coeffs[j] * eval_at.pow(&[3 - j as u64, 0, 0, 0]);
            }
            evals[i] = eval_i;
        }

        let uni_poly = UniPoly::interpolate(&evals);
        let eval = uni_poly.eval(eval_at);
        assert_eq!(eval, expected_eval);
    }
}
