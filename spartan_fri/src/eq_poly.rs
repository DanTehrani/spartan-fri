use pasta_curves::arithmetic::FieldExt;

pub struct EqPoly<F: FieldExt> {
    t: Vec<F>,
}

impl<F: FieldExt> EqPoly<F> {
    pub fn new(t: Vec<F>) -> Self {
        Self { t }
    }

    pub fn eval(&self, x: &[F]) -> F {
        let mut result = F::one();
        let one = F::one();

        for i in 0..x.len() {
            result *= self.t[i] * x[i] + (one - self.t[i]) * (one - x[i]);
        }
        result
    }
}
