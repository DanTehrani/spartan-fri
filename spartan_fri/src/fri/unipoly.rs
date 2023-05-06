use super::fft::{fft, ifft};
use pasta_curves::arithmetic::FieldExt;
use pasta_curves::group::ff::PrimeField;
use std::ops::Div;
use std::process::Output;

pub struct UniPoly<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub coeffs: Vec<F>,
}

impl<F> UniPoly<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub fn new(coeffs: Vec<F>) -> Self {
        Self { coeffs } // [x^0, x^1, x^2, x^3...]
    }

    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    pub fn eval(&self, x: F) -> F {
        let mut result = F::zero();
        for (i, coeff) in self.coeffs.iter().enumerate() {
            result += *coeff * x.pow(&[i as u64, 0, 0, 0]);
        }

        result
    }

    pub fn interpolate(domain: &[F], evals: &[F]) -> Self {
        assert!(domain.len() == evals.len());
        let coeffs = ifft(&domain, &evals);
        let mut degree = 0;

        for i in 0..coeffs.len() {
            if coeffs[i] != F::zero() {
                degree = i;
            }
        }

        Self {
            coeffs: coeffs[..(degree + 1)].to_vec(),
        }
    }
}

/*
impl<F: FieldExt<Repr = [> ; 32]>> Div for UniPoly<F {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let evals = fft(domain, self.coeffs.clone());
    }
}
 */

#[cfg(test)]
mod tests {
    use super::*;
    use pasta_curves::Fp;

    #[test]
    fn test_interpolate() {
        let coeffs = vec![
            Fp::from(1),
            Fp::from(2),
            Fp::from(3),
            Fp::from(4),
            Fp::from(5),
        ];

        let mut domain = vec![];
        let root_of_unity = Fp::root_of_unity();

        let subgroup_order = (coeffs.len() * 2).next_power_of_two();

        // Generator for the subgroup with order _subgroup_order_ in the field
        let generator = root_of_unity.pow(&[
            2u32.pow(32 - ((subgroup_order as f64).log2() as u32)) as u64,
            0,
            0,
            0,
        ]);

        for i in 0..subgroup_order {
            domain.push(generator.pow(&[i as u64, 0, 0, 0]));
        }

        let poly = UniPoly {
            coeffs: coeffs.clone(),
        };

        let mut evals = vec![];
        for val in &domain {
            evals.push(poly.eval(*val));
        }

        let interpolant = UniPoly::interpolate(&domain, &evals);

        assert!(interpolant.coeffs == poly.coeffs);
        assert!(interpolant.degree() == poly.degree());
    }
}
