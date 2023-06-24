use ark_std::{end_timer, iterable::Iterable, start_timer};

use crate::fri::tree::CommittedMerkleTree;
use crate::transcript::Transcript;

use super::fft::{fft, ifft};
use crate::FieldExt;

#[derive(Clone)]
pub struct CommittedSparseUniPoly<F>
where
    F: FieldExt,
{
    pub coeffs: Vec<(usize, F)>,
    pub degree: usize,
    pub committed_merkle_tree: CommittedMerkleTree<F>,
}

impl<F> CommittedSparseUniPoly<F>
where
    F: FieldExt,
{
    pub fn to_even_and_odd_coeffs(&self) -> (Vec<(usize, F)>, Vec<(usize, F)>) {
        let n = self.coeffs.len() / 2;
        let mut even_coeffs = Vec::with_capacity(n);
        let mut odd_coeffs = Vec::with_capacity(n);

        for coeff in &self.coeffs {
            if coeff.0 & 1 == 0 {
                even_coeffs.push((coeff.0 / 2, coeff.1));
            } else {
                odd_coeffs.push(((coeff.0 - 1) / 2, coeff.1));
            }
        }

        (even_coeffs, odd_coeffs)
    }

    pub fn eval(&self, x: F) -> F {
        let mut result = F::ZERO;
        for coeff in &self.coeffs {
            result += coeff.1 * x.pow(&[coeff.0 as u64, 0, 0, 0]);
        }

        result
    }

    pub fn eval_with_table(&self, x_powers_table: &[F]) -> F {
        debug_assert!(x_powers_table.len() >= self.coeffs.len());
        let mut result = F::ZERO;
        for coeff in &self.coeffs {
            result += coeff.1 * x_powers_table[coeff.0];
        }
        result
    }

    pub fn codeword(&self) -> Vec<F> {
        self.committed_merkle_tree.leaves.clone()
    }

    pub fn to_deep_poly(&self, z: F) -> SparseUniPoly<F> {
        // Divide this polynomial by (x - z) using [synthetic division](https://en.wikipedia.org/wiki/Synthetic_division)
        let mut sparse_quotient_coeffs = vec![F::ZERO; self.degree];
        sparse_quotient_coeffs[self.coeffs[0].0] = self.coeffs[0].1;

        let mut tmp = vec![F::ZERO; self.degree];

        tmp[self.coeffs[0].0 + 1] = self.coeffs[0].1 * z;

        for coeff in self.coeffs.iter().rev() {
            sparse_quotient_coeffs[coeff.0] = sparse_quotient_coeffs[coeff.0 - 1] * z + coeff.1;
        }

        let dense_quotient_coeffs = sparse_quotient_coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| (i, *c))
            .filter(|(_, c)| *c != F::ZERO)
            .collect::<Vec<(usize, F)>>();

        SparseUniPoly::new(dense_quotient_coeffs, self.degree - 1)
    }
}

#[derive(Clone)]
pub struct SparseUniPoly<F>
where
    F: FieldExt,
{
    pub coeffs: Vec<(usize, F)>,
    pub degree: usize,
}

impl<F> SparseUniPoly<F>
where
    F: FieldExt,
{
    pub fn new(coeffs: Vec<(usize, F)>, degree: usize) -> Self {
        Self { coeffs, degree } // [x^0, x^1, x^2, x^3...]
    }

    pub fn to_even_and_odd_coeffs(&self) -> (Vec<(usize, F)>, Vec<(usize, F)>) {
        let n = self.coeffs.len() / 2;
        let mut even_coeffs = Vec::with_capacity(n);
        let mut odd_coeffs = Vec::with_capacity(n);

        for coeff in &self.coeffs {
            if coeff.0 & 1 == 0 {
                even_coeffs.push((coeff.0 / 2, coeff.1));
            } else {
                odd_coeffs.push(((coeff.0 - 1) / 2, coeff.1));
            }
        }

        (even_coeffs, odd_coeffs)
    }

    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    pub fn eval(&self, x: F) -> F {
        let mut result = F::ZERO;
        for coeff in self.coeffs.iter() {
            result += coeff.1 * x.pow(&[coeff.0 as u64, 0, 0, 0]);
        }

        result
    }

    pub fn merkle_commit(
        &self,
        domain_generator: F,
        expansion_factor: usize,
        transcript: &mut Transcript<F>,
    ) -> CommittedSparseUniPoly<F> {
        // TODO: Use sparse FFT
        let evals = (0..((self.degree + 1) * expansion_factor))
            .map(|i| {
                let x = domain_generator.pow(&[i as u64, 0, 0, 0]);
                self.eval(x)
            })
            .collect::<Vec<F>>();
        let tree = CommittedMerkleTree::from_leaves(evals);
        transcript.append_bytes(&tree.root());

        CommittedSparseUniPoly {
            degree: self.degree,
            coeffs: self.coeffs.clone(),
            committed_merkle_tree: tree,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fri::unipoly::UniPoly;
    use ff::Field;
    use pasta_curves::Fp;

    type F = Fp;

    #[test]
    fn test_sparse_unipoly_eval() {
        let degree = 10;
        let sparse_coeffs = vec![
            (0, F::from(1u64)),
            (2, F::from(2u64)),
            (degree, F::from(3u64)),
        ];

        let mut coeffs = vec![F::ZERO; degree + 1];
        for sparse_coeff in &sparse_coeffs {
            coeffs[sparse_coeff.0] = sparse_coeff.1;
        }

        let dense_poly = UniPoly::new(coeffs);
        let sparse_poly = SparseUniPoly::new(sparse_coeffs, degree);

        for i in 0..100 {
            let x = F::from(i as u64);
            let dense_eval = dense_poly.eval(x);
            let sparse_eval = sparse_poly.eval(x);
            assert_eq!(dense_eval, sparse_eval, "x = {}", i);
        }
    }
}
