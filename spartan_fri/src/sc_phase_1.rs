use crate::utils::boolean_hypercube;
use crate::{eq_poly::EqPoly, ml_poly::MlPoly};
use pasta_curves::arithmetic::FieldExt;

pub struct UniPoly3<F: FieldExt> {
    pub coeffs: [F; 4],
}

impl<F: FieldExt> UniPoly3<F> {
    pub fn eval(&self, x: F) -> F {
        // ax^3 + bx^2 + cx + d
        let x_sq = x.square();
        let x_cub = x_sq * x;

        let a = self.coeffs[0];
        let b = self.coeffs[1];
        let c = self.coeffs[2];
        let d = self.coeffs[3];

        a * x_cub + b * x_sq + c * x + d
    }

    pub fn interpolate(evals: &[F; 4]) -> Self {
        // ax^3 + bx^2 + cx + d
        let two_inv = F::from(2u64).invert().unwrap();
        let six_inv = F::from(6u64).invert().unwrap();

        let d = evals[0];
        let a = six_inv
            * (evals[3] - evals[2] - evals[2] - evals[2] + evals[1] + evals[1] + evals[1]
                - evals[0]);
        let b = two_inv
            * (evals[0] + evals[0] - evals[1] - evals[1] - evals[1] - evals[1] - evals[1]
                + evals[2]
                + evals[2]
                + evals[2]
                + evals[2]
                - evals[3]);

        let c = evals[1] - d - a - b;

        Self {
            coeffs: [a, b, c, d],
        }
    }
}

pub struct SCPhase1Proof<F: FieldExt> {
    pub round_polys: Vec<UniPoly3<F>>,
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

    pub fn round(&self, j: usize) -> UniPoly3<F> {
        // evaluate at points 0, 1, 2, 3

        let zero = F::zero();
        let one = F::one();

        let m = self.Az_poly.num_vars;
        let mut evals = [F::zero(); 4];

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

        let round_poly = UniPoly3::interpolate(&evals);
        round_poly
    }

    pub fn prove(&self) -> SCPhase1Proof<F> {
        let num_vars = self.Az_poly.num_vars;
        let mut round_polys = Vec::<UniPoly3<F>>::with_capacity(num_vars);

        for i in 0..num_vars {
            let round_poly = self.round(i);
            round_polys.push(round_poly);
        }

        SCPhase1Proof { round_polys }
    }

    pub fn verify_round_polys(proof: &SCPhase1Proof<F>, challenge: &[F]) -> (F, F) {
        let m = challenge.len();

        let zero = F::zero();
        let one = F::one();

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

    #[test]
    fn test_unipoly_3() {
        let coeffs = [F::from(1u64), F::from(2u64), F::from(3u64), F::from(4u64)];
        let eval_at = Fp::from(33);

        let mut expected_eval = F::zero();
        for i in 0..coeffs.len() {
            expected_eval += coeffs[i] * eval_at.pow(&[3 - i as u64, 0, 0, 0]);
        }

        let mut evals = [F::zero(); 4];
        for i in 0..4 {
            let eval_at = F::from(i as u64);
            let mut eval_i = F::zero();
            for j in 0..coeffs.len() {
                eval_i += coeffs[j] * eval_at.pow(&[3 - j as u64, 0, 0, 0]);
            }
            evals[i] = eval_i;
        }

        let uni_poly = UniPoly3::interpolate(&evals);
        let eval = uni_poly.eval(eval_at);
        assert_eq!(eval, expected_eval);
    }
}
