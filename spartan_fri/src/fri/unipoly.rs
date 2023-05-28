use crate::fri::tree::CommittedMerkleTree;
use crate::transcript::Transcript;

use super::fft::{fft, ifft};
use crate::FieldExt;

#[derive(Clone)]
pub struct CommittedUniPoly<F>
where
    F: FieldExt,
{
    pub coeffs: Vec<F>,
    pub committed_merkle_tree: CommittedMerkleTree<F>,
}

impl<F> CommittedUniPoly<F>
where
    F: FieldExt,
{
    pub fn to_even_and_odd_coeffs(&self) -> (Vec<F>, Vec<F>) {
        let n = self.coeffs.len() / 2;
        let mut even_coeffs = Vec::with_capacity(n);
        let mut odd_coeffs = Vec::with_capacity(n);

        for (i, coeff) in self.coeffs.iter().enumerate() {
            if i % 2 == 0 {
                even_coeffs.push(*coeff);
            } else {
                odd_coeffs.push(*coeff);
            }
        }

        (even_coeffs, odd_coeffs)
    }

    pub fn eval(&self, x: F) -> F {
        let mut result = F::ZERO;
        for (i, coeff) in self.coeffs.iter().enumerate() {
            result += *coeff * x.pow(&[i as u64, 0, 0, 0]);
        }

        result
    }

    pub fn codeword(&self) -> Vec<F> {
        self.committed_merkle_tree.layers[0].clone()
    }
}

#[derive(Clone)]
pub struct UniPoly<F>
where
    F: FieldExt,
{
    pub coeffs: Vec<F>,
}

impl<F> UniPoly<F>
where
    F: FieldExt,
{
    pub fn new(coeffs: Vec<F>) -> Self {
        Self { coeffs } // [x^0, x^1, x^2, x^3...]
    }

    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    pub fn eval_fft(&self, domain: &[F]) -> Vec<F> {
        // Pad the coefficients vector with zeros to match the domain size
        let mut padded_coeffs = self.coeffs.clone();
        let mut pad = Vec::with_capacity(domain.len());
        pad.resize(domain.len() - self.coeffs.len(), F::ZERO);
        padded_coeffs.extend_from_slice(&pad);

        // Evaluate with FFT
        fft(&padded_coeffs, domain)
    }

    pub fn eval(&self, x: F) -> F {
        let mut result = F::ZERO;
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
            if coeffs[i] != F::ZERO {
                degree = i;
            }
        }

        Self {
            coeffs: coeffs[..(degree + 1)].to_vec(),
        }
    }

    pub fn merkle_commit(
        &self,
        domain: &[F],
        transcript: &mut Transcript<F>,
    ) -> CommittedUniPoly<F> {
        let evals = self.eval_fft(domain);
        let tree = CommittedMerkleTree::from_leaves(evals);
        transcript.append_fe(&tree.root());

        CommittedUniPoly {
            coeffs: self.coeffs.clone(),
            committed_merkle_tree: tree,
        }
    }
}

/*
impl<F: PrimeField<Repr = [> ; 32]>> Div for UniPoly<F {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let evals = fft(domain, self.coeffs.clone());
    }
}
 */

#[cfg(test)]
mod tests {
    use super::*;
    use ff::Field;
    use pasta_curves::group::ff::PrimeField;
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
        let root_of_unity = Fp::ROOT_OF_UNITY;

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
