use crate::utils::boolean_hypercube;
use crate::{eq_poly::EqPoly, ml_poly::MlPoly};
use pasta_curves::{arithmetic::FieldExt, Eq};

pub struct UniPoly2<F: FieldExt> {
    pub coeffs: [F; 3],
}

impl<F: FieldExt> UniPoly2<F> {
    pub fn eval(&self, x: F) -> F {
        // ax^3 + bx^2 + cx + d
        let x_sq = x.square();

        let a = self.coeffs[0];
        let b = self.coeffs[1];
        let c = self.coeffs[2];

        a * x_sq + b * x + c
    }

    pub fn interpolate(evals: &[F; 3]) -> Self {
        // ax^2 + bx + c
        let two_inv = F::TWO_INV;

        // y = ax^2 + bx + c

        // eval at 0
        // y_0 = c

        // eval at 1
        // y_1 = a + b + c
        // b = y_1 - a - c

        // eval at 2
        // y_2 = 4a + 2b + c

        // 4a + 2b + c - (a + b + c) = y_2  - y1
        // -> 3a + b = y_2  - y1

        // a + b + c - c = y_1 - y_0
        // -> a + b = y_1 - y_0

        // 3a + b - (y_2  - y1) = a + b - (y_1 - y_0)
        // 3a - (y_2  - y1) = a - (y_1 - y_0)
        // 2a = (y_2  - y1) - (y_1 - y_0)
        // a = (y_2 - 2y1) + y_0 / 2

        let c = evals[0];
        let a = (evals[2] - evals[1] - evals[1] + evals[0]) * two_inv;
        let b = evals[1] - a - c;

        Self { coeffs: [a, b, c] }
    }
}

pub struct SCPhase2Proof<F: FieldExt> {
    pub round_polys: Vec<UniPoly2<F>>,
}

pub struct SumCheckPhase2<F: FieldExt> {
    Az_poly: MlPoly<F>,
    Bz_poly: MlPoly<F>,
    Cz_poly: MlPoly<F>,
    r: [F; 3],
    bound_eq_poly: EqPoly<F>,
    challenge: Vec<F>,
}

impl<F: FieldExt> SumCheckPhase2<F> {
    pub fn new(
        Az_poly: MlPoly<F>,
        Bz_poly: MlPoly<F>,
        Cz_poly: MlPoly<F>,
        r: [F; 3],
        challenge: Vec<F>,
    ) -> Self {
        let bound_eq_poly = EqPoly::new(challenge.clone());
        Self {
            Az_poly,
            Bz_poly,
            Cz_poly,
            r,
            bound_eq_poly,
            challenge,
        }
    }

    fn eval_poly(&self, vars: &[F]) -> F {
        assert_eq!(vars.len(), self.Az_poly.num_vars);

        let r_A = self.r[0];
        let r_B = self.r[1];
        let r_C = self.r[2];

        let eval = self.Az_poly.eval(vars) * r_A
            + self.Bz_poly.eval(vars) * r_B
            + self.Cz_poly.eval(vars) * r_C;

        eval
    }

    pub fn round(&self, j: usize) -> UniPoly2<F> {
        // evaluate at points 0, 1, 2, 3

        let zero = F::zero();
        let one = F::one();

        let m = self.Az_poly.num_vars;
        let mut evals = [F::zero(); 3];

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
        }

        let round_poly = UniPoly2::interpolate(&evals);
        round_poly
    }

    pub fn prove(&self) -> SCPhase2Proof<F> {
        let num_vars = self.Az_poly.num_vars;
        let mut round_polys = Vec::<UniPoly2<F>>::with_capacity(num_vars);

        for i in 0..num_vars {
            let round_poly = self.round(i);
            round_polys.push(round_poly);
        }

        SCPhase2Proof { round_polys }
    }

    pub fn verify_round_polys(proof: &SCPhase2Proof<F>, challenge: &[F]) -> (F, F) {
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
    fn test_unipoly_2() {
        let coeffs = [F::from(1u64), F::from(2u64), F::from(3u64)];
        let eval_at = Fp::from(33);

        let mut expected_eval = F::zero();
        for i in 0..coeffs.len() {
            expected_eval += coeffs[i] * eval_at.pow(&[i as u64, 0, 0, 0]);
        }

        let mut evals = [F::zero(); 3];
        for i in 0..3 {
            let eval_at = F::from(i as u64);
            let mut eval_i = F::zero();
            for j in 0..coeffs.len() {
                eval_i += coeffs[j] * eval_at.pow(&[j as u64, 0, 0, 0]);
            }
            evals[i] = eval_i;
        }

        let uni_poly = UniPoly2::interpolate(&evals);
        let eval = uni_poly.eval(eval_at);
        assert_eq!(eval, expected_eval);
    }
}
