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
        tau: Vec<F>,
        challenge: Vec<F>,
    ) -> Self {
        let bound_eq_poly = EqPoly::new(tau);
        Self {
            Az_poly,
            Bz_poly,
            Cz_poly,
            bound_eq_poly,
            challenge,
        }
    }

    pub fn prove(&self) -> SCPhase1Proof<F> {
        let num_vars = self.Az_poly.num_vars;
        let mut round_polys = Vec::<UniPoly<F>>::with_capacity(num_vars - 1);

        let mut A_table = self.Az_poly.evals.clone();
        let mut B_table = self.Bz_poly.evals.clone();
        let mut C_table = self.Cz_poly.evals.clone();
        let mut eq_table = boolean_hypercube(num_vars)
            .iter()
            .map(|b| self.bound_eq_poly.eval(b))
            .collect::<Vec<F>>();

        for j in 0..num_vars {
            let zero = F::ZERO;
            let one = F::ONE;

            let m = self.Az_poly.num_vars;

            let high_index = 2usize.pow((m - j - 1) as u32);

            let mut evals = [F::ZERO; 4];

            // https://eprint.iacr.org/2019/317.pdf#subsection.3.2
            for b in 0..high_index {
                for (i, eval_at) in [zero, one, F::from(2), F::from(3)].iter().enumerate() {
                    let a_eval = A_table[b] + (A_table[b + high_index] - A_table[b]) * eval_at;
                    let b_eval = B_table[b] + (B_table[b + high_index] - B_table[b]) * eval_at;
                    let c_eval = C_table[b] + (C_table[b + high_index] - C_table[b]) * eval_at;
                    let eq_eval = eq_table[b] + (eq_table[b + high_index] - eq_table[b]) * eval_at;
                    evals[i] += (a_eval * b_eval - c_eval) * eq_eval;
                }

                let r_i = self.challenge[j];
                A_table[b] = A_table[b] + (A_table[b + high_index] - A_table[b]) * r_i;
                B_table[b] = B_table[b] + (B_table[b + high_index] - B_table[b]) * r_i;
                C_table[b] = C_table[b] + (C_table[b + high_index] - C_table[b]) * r_i;
                eq_table[b] = eq_table[b] + (eq_table[b + high_index] - eq_table[b]) * r_i;
            }

            let round_poly = UniPoly::interpolate(&evals);

            round_polys.push(round_poly);
        }

        SCPhase1Proof { round_polys }
    }

    pub fn verify_round_polys(proof: &SCPhase1Proof<F>, challenge: &[F]) -> F {
        debug_assert_eq!(proof.round_polys.len(), challenge.len());

        let zero = F::ZERO;
        let one = F::ONE;

        let mut target = zero;
        for (i, round_poly) in proof.round_polys.iter().enumerate() {
            assert_eq!(
                round_poly.eval(zero) + round_poly.eval(one),
                target,
                "round poly {} failed",
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
