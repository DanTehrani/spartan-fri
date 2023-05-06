use pasta_curves::{arithmetic::FieldExt, Eq};

#[derive(Clone, Debug)]
pub struct MlPoly<F> {
    pub evals: Vec<F>,
    pub num_vars: usize,
}

impl<F: FieldExt> MlPoly<F> {
    pub fn new(evals: Vec<F>) -> Self {
        assert!(evals.len().is_power_of_two());
        let num_vars = (evals.len() as f64).log2() as usize;
        Self { evals, num_vars }
    }

    // Evaluate the multilinear extension of the polynomial `a`, at point `t`.
    // `a` is in evaluation form.
    pub fn eval(&self, t: &[F]) -> F {
        let n = self.evals.len();
        assert_eq!((n as f64).log2() as usize, t.len());
        // Evaluate the multilinear extension of the polynomial `a`,
        // over the boolean hypercube

        let m = t.len();

        /*
        let mut result = F::zero();
        let one = F::one();

        for i in 0..n {
            let mut result_i = F::one();
            for j in 0..m {
                let i_b = F::from((i >> j & 1) as u64);
                result_i *= t[j] * i_b + (one - t[j]) * (one - i_b);
            }

            result += self.evals[i] * result_i;
        }

        result
         */

        // Borrowed from https://github.com/microsoft/Spartan/blob/1e431e2bbfe74b8a53488c43adf32de1dd974777/src/dense_mlpoly.rs#L67
        let mut evals = vec![F::one(); n];

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
        let mut result = F::zero();
        for i in 0..n {
            result += self.evals[i];
        }

        result
    }
}
