use crate::{spartan::utils::boolean_hypercube, FieldExt};
use ff::Field;

#[derive(Clone, Debug)]
pub struct SparseMLPoly<F> {
    pub evals: Vec<(usize, F)>,
    pub coeffs: Option<Vec<(usize, F)>>,
    pub num_vars: usize,
}

impl<F: FieldExt> SparseMLPoly<F> {
    pub fn new(evals: Vec<(usize, F)>, num_vars: usize) -> Self {
        Self {
            evals,
            coeffs: None,
            num_vars,
        }
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

    pub fn from_coeffs(coeffs: Vec<(usize, F)>) -> Self {
        let n = coeffs.len();
        assert!(n.is_power_of_two());
        let num_vars = (n as f64).log2() as usize;
        let domain: Vec<Vec<F>> = boolean_hypercube::<F>(num_vars);

        let evals = domain
            .iter()
            .enumerate()
            .map(|(i, x)| {
                let mut result = F::ZERO;
                for i in 0..n {
                    let mut term = coeffs[i].1;
                    for j in 0..num_vars {
                        if (i >> j & 1) == 1 {
                            term *= x[j];
                        }
                    }
                    result += term;
                }
                (i, result)
            })
            .collect::<Vec<(usize, F)>>();

        Self {
            evals,
            coeffs: Some(coeffs),
            num_vars,
        }
    }

    pub fn to_coeffs(&self) -> Vec<(usize, F)> {
        let mut coeffs: Vec<(usize, F)> = vec![];
        for eval in &self.evals {
            let mut coeff = eval.1;
            for j in 0..coeffs.len() {
                if j & eval.0 == j {
                    coeff -= coeffs[j].1;
                }
            }
            coeffs.push((eval.0, coeff));
        }

        coeffs
    }

    pub fn eval_with_coeffs(&self, vals: &[F]) -> F {
        let mut result = F::ZERO;
        for coeff in self.coeffs.as_ref().unwrap() {
            let mut term = coeff.1;
            for j in 0..self.num_vars {
                if (coeff.0 >> j & 1) == 1 {
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
        for coeff in coeffs {
            result += coeff.1 * val.pow(&[coeff.0 as u64, 0, 0, 0]);
        }

        result
    }

    pub fn eval_as_uni_with_evals(&self, val: &F) -> F {
        let n = self.evals.len();

        let coeffs = self.to_coeffs();

        let mut result = F::ZERO;
        for i in 0..n {
            result += coeffs[i].1 * val.pow(&[i as u64, 0, 0, 0]);
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

    pub fn eval(&self, t: &[F]) -> F {
        // Evaluate the multilinear extension of the polynomial `a`,
        // over the boolean hypercube

        let m = t.len();

        let mut result = F::ZERO;
        let one = F::ONE;

        for eval in &self.evals {
            let mut result_i = F::ONE;
            for j in 0..m {
                let i_b = F::from((eval.0 >> j & 1) as u64);
                result_i *= t[j] * i_b + (one - t[j]) * (one - i_b);
            }

            result += eval.1 * result_i;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MlPoly;
    type F = pasta_curves::Fp;

    #[test]
    pub fn test_ml_eval_uni() {
        let m = 4;
        let n = 2usize.pow(m as u32);
        let coeffs = (0..n)
            .map(|i| (i, F::from((i + 123) as u64)))
            .collect::<Vec<(usize, F)>>();

        let ml_poly = SparseMLPoly::from_coeffs(coeffs.clone());

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
    pub fn test_to_coeffs() {
        let num_vars = 3;
        let dense_evals = vec![3, 33, 0, 333, 0, 0, 0, 3333]
            .iter()
            .map(|e| F::from(*e as u64))
            .collect::<Vec<F>>();
        let sparse_evals = dense_evals
            .iter()
            .enumerate()
            .map(|(i, e)| (i, *e))
            .filter(|(_, e)| *e != F::ZERO)
            .collect::<Vec<(usize, F)>>();

        let mut sparse_ml_poly = SparseMLPoly::new(sparse_evals, num_vars);
        let mut dense_ml_poly = MlPoly::new(dense_evals);

        sparse_ml_poly.compute_coeffs();
        dense_ml_poly.compute_coeffs();

        let eval_at = vec![F::from(3), F::from(33), F::from(333)];
        let sparse_eval = sparse_ml_poly.eval_with_coeffs(&eval_at);
        let dense_eval = dense_ml_poly.eval_with_coeffs(&eval_at);

        println!("dense coeffs: {:?}", dense_ml_poly.coeffs.unwrap());
        println!("sparse eval: {:?}", sparse_eval);
        println!("dense eval: {:?}", dense_eval);
    }
}
