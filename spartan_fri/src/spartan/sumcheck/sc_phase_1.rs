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
    Az_poly: MlPoly<F>,
    Bz_poly: MlPoly<F>,
    Cz_poly: MlPoly<F>,
    bound_eq_poly: EqPoly<F>,
    challenge: Vec<F>,
}

impl<F: FieldExt> SumCheckPhase1<F> {
    pub fn new(
        Az_poly: MlPoly<F>,
        Bz_poly: MlPoly<F>,
        Cz_poly: MlPoly<F>,
        challenge: Vec<F>,
    ) -> Self {
        let bound_eq_poly = EqPoly::new(challenge.clone());
        Self {
            Az_poly,
            Bz_poly,
            Cz_poly,
            bound_eq_poly,
            challenge,
        }
    }

    fn eval_poly(&self, vars: &[F]) -> F {
        assert_eq!(vars.len(), self.Az_poly.num_vars);

        let eval = (self.Az_poly.eval(vars) + self.Bz_poly.eval(vars) - self.Cz_poly.eval(vars))
            * self.bound_eq_poly.eval(vars);

        eval
    }

    pub fn round(&self, j: usize) -> UniPoly<F> {
        // evaluate at points 0, 1, 2, 3

        let zero = F::ZERO;
        let one = F::ONE;

        let m = self.Az_poly.num_vars;
        let mut evals = [F::ZERO; 4];

        for vars in &boolean_hypercube(m - j - 1) {
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
        let num_vars = self.Az_poly.num_vars;
        let mut round_polys = Vec::<UniPoly<F>>::with_capacity(num_vars);

        for i in 0..num_vars {
            let round_poly = self.round(i);
            round_polys.push(round_poly);
        }

        SCPhase1Proof { round_polys }
    }

    pub fn verify_round_polys(proof: &SCPhase1Proof<F>, challenge: &[F]) -> (F, F) {
        let m = challenge.len();

        let zero = F::ZERO;
        let one = F::ONE;

        let round_polys = &proof.round_polys;
        let sum_claim = round_polys[0].eval(zero) + round_polys[0].eval(one);

        for i in 0..(m - 1) {
            let round_poly = &round_polys[i];
            let round_poly_eval = round_poly.eval(challenge[i]);
            let next_round_poly = &round_polys[i + 1];
            let eval_0 = next_round_poly.eval(zero);
            let eval_1 = next_round_poly.eval(one);

            assert_eq!(round_poly_eval, eval_0 + eval_1);
        }

        let final_poly_eval = round_polys[m - 1].eval(challenge[m - 1]);

        (sum_claim, final_poly_eval)
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
