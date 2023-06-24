use crate::FieldExt;

pub struct EqPoly<F: FieldExt> {
    t: Vec<F>,
}

impl<F: FieldExt> EqPoly<F> {
    pub fn new(t: Vec<F>) -> Self {
        Self { t }
    }

    pub fn eval(&self, x: &[F]) -> F {
        let mut result = F::ONE;
        let one = F::ONE;

        for i in 0..x.len() {
            result *= self.t[i] * x[i] + (one - self.t[i]) * (one - x[i]);
        }
        result
    }

    // Copied from microsoft/Spartan
    pub fn evals(&self) -> Vec<F> {
        let ell = self.t.len();

        let mut evals: Vec<F> = vec![F::ONE; 2usize.pow(ell as u32)];
        let mut size = 1;
        for j in 0..ell {
            // in each iteration, we double the size of chis
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                // copy each element from the prior iteration twice
                let scalar = evals[i / 2];
                evals[i] = scalar * self.t[j];
                evals[i - 1] = scalar - evals[i];
            }
        }
        evals
    }
}
