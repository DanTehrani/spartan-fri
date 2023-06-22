// Phase 2 sum-check of Spartan.

use crate::spartan::polynomial::ml_poly::MlPoly;
use crate::spartan::sumcheck::unipoly::UniPoly;
use crate::spartan::utils::boolean_hypercube;
use crate::FieldExt;

pub struct SCPhase2Proof<F: FieldExt> {
    pub round_polys: Vec<UniPoly<F>>,
}

pub struct SumCheckPhase2<F: FieldExt> {
    A_mle: MlPoly<F>,
    B_mle: MlPoly<F>,
    C_mle: MlPoly<F>,
    Z_mle: MlPoly<F>,
    rx: Vec<F>,
    r: [F; 3],
    challenge: Vec<F>,
}

impl<F: FieldExt> SumCheckPhase2<F> {
    pub fn new(
        A_mle: MlPoly<F>,
        B_mle: MlPoly<F>,
        C_mle: MlPoly<F>,
        Z_mle: MlPoly<F>,
        rx: Vec<F>,
        r: [F; 3],
        challenge: Vec<F>,
    ) -> Self {
        Self {
            A_mle,
            B_mle,
            C_mle,
            Z_mle,
            rx,
            r,
            challenge,
        }
    }

    fn eval_poly(&self, vars: &[F]) -> F {
        debug_assert_eq!(vars.len(), self.A_mle.num_vars / 2);

        let r_A = self.r[0];
        let r_B = self.r[1];
        let r_C = self.r[2];

        let z_eval = self.Z_mle.eval(&vars);
        let mut r_x = self.rx.clone();
        r_x.reverse();
        let eval_at = [&vars, r_x.as_slice()].concat();
        let a_eval = self.A_mle.eval(&eval_at) * z_eval;
        let b_eval = self.B_mle.eval(&eval_at) * z_eval;
        let c_eval = self.C_mle.eval(&eval_at) * z_eval;

        a_eval * r_A + b_eval * r_B + c_eval * r_C
    }

    pub fn round(&self, j: usize) -> UniPoly<F> {
        // evaluate at points 0, 1, 2, 3

        let zero = F::ZERO;
        let one = F::ONE;

        let m = self.A_mle.num_vars / 2;
        let mut evals = [F::ZERO; 3];

        for vars in &boolean_hypercube(m - j - 1) {
            let mut c = self.challenge[..j].to_vec();
            c.reverse();
            let mut eval_at = vec![vars.as_slice(), &[zero], &c].concat();

            // Eval at 0
            evals[0] += self.eval_poly(&eval_at);

            // Eval at 1
            eval_at[vars.len()] = one;
            evals[1] += self.eval_poly(&eval_at);

            // Eval at 2
            eval_at[vars.len()] = F::from(2u64);
            evals[2] += self.eval_poly(&eval_at);
        }

        let round_poly = UniPoly::interpolate(&evals);
        round_poly
    }

    pub fn prove(&self) -> SCPhase2Proof<F> {
        let num_vars = self.A_mle.num_vars / 2;
        let mut round_polys = Vec::<UniPoly<F>>::with_capacity(num_vars - 1);

        for i in 0..num_vars {
            let round_poly = self.round(i);
            round_polys.push(round_poly);
        }

        SCPhase2Proof { round_polys }
    }

    pub fn verify_round_polys(sum_target: F, proof: &SCPhase2Proof<F>, challenge: &[F]) -> F {
        debug_assert_eq!(proof.round_polys.len(), challenge.len());

        let zero = F::ZERO;
        let one = F::ONE;

        let mut target = sum_target;
        for (i, round_poly) in proof.round_polys.iter().enumerate() {
            assert_eq!(
                round_poly.eval(zero) + round_poly.eval(one),
                target,
                "i = {}",
                i
            );
            target = round_poly.eval(challenge[i]);
        }

        target
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff::Field;
    use pasta_curves::Fp;
    type F = Fp;

    #[test]
    fn test_unipoly_2() {
        let coeffs = [F::from(1u64), F::from(2u64), F::from(3u64)];
        let eval_at = Fp::from(33);

        let mut expected_eval = F::ZERO;
        for i in 0..coeffs.len() {
            expected_eval += coeffs[i] * eval_at.pow(&[i as u64, 0, 0, 0]);
        }

        let mut evals = [F::ZERO; 3];
        for i in 0..3 {
            let eval_at = F::from(i as u64);
            let mut eval_i = F::ZERO;
            for j in 0..coeffs.len() {
                eval_i += coeffs[j] * eval_at.pow(&[j as u64, 0, 0, 0]);
            }
            evals[i] = eval_i;
        }

        let uni_poly = UniPoly::interpolate(&evals);
        let eval = uni_poly.eval(eval_at);
        assert_eq!(eval, expected_eval);
    }
}
