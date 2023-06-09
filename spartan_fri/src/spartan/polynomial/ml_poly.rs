use crate::{spartan::utils::boolean_hypercube, FieldExt};
use ff::Field;

#[derive(Clone, Debug)]
pub struct MlPoly<F> {
    pub evals: Vec<F>,
    pub coeffs: Option<Vec<F>>,
    pub num_vars: usize,
}

impl<F: FieldExt> MlPoly<F> {
    pub fn new(evals: Vec<F>) -> Self {
        assert!(evals.len().is_power_of_two());
        let num_vars = (evals.len() as f64).log2() as usize;
        Self {
            evals,
            num_vars,
            coeffs: None,
        }
    }

    pub fn bound_high(&self, x: &[F]) -> MlPoly<F> {
        let new_n = 2usize.pow((self.num_vars - x.len()) as u32);
        let m = self.num_vars;
        let new_m = self.num_vars - x.len();

        let mut evals = self.evals[..new_n].to_vec();
        for i in 0..new_n {
            for j in new_m..m {
                if (i >> j & 1) == 1 {
                    evals[i] *= x[j - new_m];
                }
            }
        }
        MlPoly::new(evals)
    }

    pub fn compute_coeffs(&mut self) {
        let coeffs = self.to_coeffs();
        self.coeffs = Some(coeffs);
    }

    pub fn eval_coeffs(&self, coeffs: &[F], vals: &[F]) -> F {
        let n = coeffs.len();
        let m = vals.len();
        debug_assert_eq!(2usize.pow(m as u32), n);

        let mut result = F::ZERO;
        for i in 0..n {
            let mut term = coeffs[i];
            for j in 0..m {
                if (i >> j & 1) == 1 {
                    term *= vals[j];
                }
            }
            result += term;
        }

        result
    }

    pub fn from_coeffs(coeffs: Vec<F>) -> Self {
        let n = coeffs.len();
        assert!(n.is_power_of_two());
        let num_vars = (n as f64).log2() as usize;
        let domain = boolean_hypercube::<F>(num_vars);

        let evals = domain
            .iter()
            .map(|x| {
                let mut result = F::ZERO;
                for i in 0..n {
                    let mut term = coeffs[i];
                    for j in 0..num_vars {
                        if (i >> j & 1) == 1 {
                            term *= x[j];
                        }
                    }
                    result += term;
                }
                result
            })
            .collect::<Vec<F>>();

        Self {
            evals,
            coeffs: Some(coeffs),
            num_vars,
        }
    }

    pub fn to_coeffs(&self) -> Vec<F> {
        let n = self.evals.len();
        let mut coeffs = Vec::with_capacity(n);
        for i in 0..n {
            let mut coeff = self.evals[i];
            for j in 0..coeffs.len() {
                if j & i == j {
                    coeff -= coeffs[j];
                }
            }
            coeffs.push(coeff);
        }

        coeffs
    }

    pub fn eval_with_coeffs(&self, vals: &[F]) -> F {
        let n = self.evals.len();
        let m = vals.len();
        debug_assert_eq!(2usize.pow(m as u32), n);

        let mut result = F::ZERO;
        let coeffs = self.coeffs.as_ref().unwrap();
        for i in 0..n {
            let mut term = coeffs[i];
            for j in 0..m {
                if (i >> j & 1) == 1 {
                    term *= vals[j];
                }
            }
            result += term;
        }

        result
    }

    // Evaluate this multilinear polynomial as a univariate polynomial
    // with the same coefficients.
    pub fn eval_as_uni_with_coeffs(&self, val: &F) -> F {
        let mut result = F::ZERO;

        let coeffs = self.coeffs.as_ref().unwrap();
        let n = coeffs.len();
        for i in 0..n {
            if coeffs[i] == F::ZERO {
                continue;
            }
            result += coeffs[i] * val.pow(&[i as u64, 0, 0, 0]);
        }

        result
    }

    pub fn eval_as_uni_with_evals(&self, val: &F) -> F {
        let n = self.evals.len();

        let coeffs = self.to_coeffs();

        let mut result = F::ZERO;
        for i in 0..n {
            result += coeffs[i] * val.pow(&[i as u64, 0, 0, 0]);
        }

        result

        /*
        let mut term_evals = Vec::with_capacity(n);

        // Evaluate per coefficient
        for i in 0..n {
            let mut vars = Vec::with_capacity(self.num_vars);
            let mut val_set = false;
            for j in 0..m {
                if (i >> j & 1) == 1 && !val_set {
                    vars.insert(j, val.pow(&[i as u64, 0, 0, 0]));
                    val_set = true;
                } else if (i >> j & 1) == 1 {
                    vars.insert(j, F::ONE);
                } else {
                    vars.insert(j, F::ZERO;);
                }
            }

            let mut eval_i = self.eval(&vars);
            for l in 0..i {
                let mask = l & i;
                for j in 0..m {
                    if (mask >> j & 1) == 1 {
                        eval_i -= term_evals[l];
                    }
                }
            }

            term_evals.push(eval_i);
        }

        let mut eval = F::ZERO;;
        for term_eval in term_evals {
            eval += term_eval;
        }
        eval
        */
    }

    // Evaluate the multilinear extension of the polynomial `a`, at point `t`.
    // `a` is in evaluation form.
    pub fn eval(&self, t: &[F]) -> F {
        let n = self.evals.len();
        debug_assert_eq!((n as f64).log2() as usize, t.len());
        // Evaluate the multilinear extension of the polynomial `a`,
        // over the boolean hypercube

        let m = t.len();

        let mut result = F::ZERO;
        let one = F::ONE;

        for i in 0..n {
            let mut result_i = F::ONE;
            for j in 0..m {
                let i_b = F::from((i >> j & 1) as u64);
                result_i *= t[j] * i_b + (one - t[j]) * (one - i_b);
            }

            result += self.evals[i] * result_i;
        }

        result

        /*
        // Borrowed from https://github.com/microsoft/Spartan/blob/1e431e2bbfe74b8a53488c43adf32de1dd974777/src/dense_mlpoly.rs#L67
        let mut evals = vec![F::ONE; n];

        let mut size = 1;
        for j in 0..m {
            // in each iteration, we double the size of chis
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                // copy each element from the prior iteration twice
                let scalar = evals[i / 2];
                evals[i] = scalar * t[j];
                evals[i - 1] = scalar - evals[i];
            }
        }

        // TODO: Make this faster!
        let mut result = F::ZERO;;
        for i in 0..n {
            result += self.evals[i];
        }

        result
         */
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type F = pasta_curves::Fp;

    #[test]
    pub fn test_ml_eval_uni() {
        let m = 4;
        let n = 2usize.pow(m as u32);
        let coeffs = (0..n)
            .map(|i| F::from((i + 123) as u64))
            .collect::<Vec<F>>();

        let ml_poly = MlPoly::from_coeffs(coeffs.clone());

        let vals = (0..m).map(|i| F::from(i + 33 as u64)).collect::<Vec<F>>();
        let eval_coeffs = ml_poly.eval_with_coeffs(&vals);
        let eval = ml_poly.eval(&vals);
        assert_eq!(eval_coeffs, eval);

        let uni_eval_at = F::from(33);
        let eval_as_uni_with_coeffs = ml_poly.eval_as_uni_with_coeffs(&uni_eval_at);
        let eval_as_uni_with_evals = ml_poly.eval_as_uni_with_evals(&uni_eval_at);
        assert_eq!(eval_as_uni_with_coeffs, eval_as_uni_with_evals);
    }

    #[test]
    pub fn test_bound_high() {
        let m = 4;
        let n = 2usize.pow(m as u32);
        let evals = (0..n)
            .map(|i| F::from((i + 123) as u64))
            .collect::<Vec<F>>();

        let ml_poly = MlPoly::new(evals);

        let bound_vars = [F::from(333), F::from(333)];
        let bounded_poly = ml_poly.bound_high(&bound_vars);

        let eval_bounded_at = [F::from(333), F::from(333)];
        let eval_at = [eval_bounded_at, bound_vars].concat();
        assert_eq!(bounded_poly.eval(&eval_bounded_at), ml_poly.eval(&eval_at));
    }
}
