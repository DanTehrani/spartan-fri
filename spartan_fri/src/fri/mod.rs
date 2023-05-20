mod fft;
mod fri_prover;
mod fri_verifier;
mod tree;
mod unipoly;
mod utils;

use pasta_curves::arithmetic::FieldExt;
use tree::MerkleProof;

use crate::transcript::Transcript;
pub use fri_prover::FRIPolyCommitProver;
pub use fri_verifier::FRIPolyCommitVerifier;
pub use unipoly::UniPoly;

use crate::PolyCommitment;

pub struct BatchedPolyEvalProof<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub openings: Vec<Vec<(MerkleProof<F>, MerkleProof<F>, MerkleProof<F>)>>,
    pub reduced_codeword: Vec<F>,
    pub poly_openings: Vec<Vec<(MerkleProof<F>, MerkleProof<F>)>>,
    pub poly_degrees: Vec<usize>,
    pub z: Vec<F>,
    pub y: Vec<F>,
}
use std::marker::PhantomData;

#[derive(Clone)]
pub struct FRIMultilinearPCS<F> {
    _marker: PhantomData<F>,
}

impl<F: FieldExt> PolyCommitment<F> for FRIMultilinearPCS<F> {
    type Commitment = F; // tmp
    type Opening = F; // tmp
    fn new() -> Self {
        todo!()
    }

    fn commit(&self, evals: &[F]) -> Self::Commitment {
        todo!()
    }

    fn open(&self, point: &[F]) -> Self::Opening {
        todo!()
    }
}

#[derive(Clone)]
pub struct FRIConfig<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub R: usize,
    pub folding_factor: usize,
    pub num_queries: usize,
    pub L: Vec<Vec<F>>,
}

impl<F> FRIConfig<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub fn new(
        poly_degree: usize,
        expansion_factor: usize,
        folding_factor: usize,
        num_queries: usize,
        final_codeword_size: usize,
    ) -> Self {
        let root_of_unity = F::root_of_unity();
        let mut domain_order = (poly_degree * expansion_factor).next_power_of_two();
        let mut L = vec![];

        let mut codeword_size = poly_degree;
        while codeword_size > final_codeword_size {
            // Generator for the subgroup with order _subgroup_order_ in the field
            let domain_generator = root_of_unity.pow(&[
                2u32.pow(32 - ((domain_order as f64).log2() as u32)) as u64,
                0,
                0,
                0,
            ]);

            let L_i = (0..domain_order)
                .map(|i| domain_generator.pow(&[i as u64, 0, 0, 0]))
                .collect();
            codeword_size /= folding_factor;
            domain_order /= folding_factor;
            L.push(L_i);
        }

        Self {
            R: expansion_factor,
            folding_factor,
            num_queries,
            L,
        }
    }
}

impl<F> FRIConfig<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub fn final_codeword_size(&self) -> usize {
        self.L[self.L.len() - 1].len()
    }

    pub fn num_rounds(&self) -> usize {
        let mut num_rounds = 0;
        let mut codeword_size = self.L[0].len();
        while codeword_size > self.final_codeword_size() {
            num_rounds += 1;
            codeword_size /= self.folding_factor;
        }
        num_rounds
    }
}

#[cfg(test)]
mod tests {
    use super::unipoly::UniPoly;
    use super::*;

    type F = pasta_curves::Fp;

    #[test]
    fn test_fri_batch_poly_commit() {
        let num_polys = 5;
        let poly_degrees = (0..num_polys)
            .map(|i| 2usize.pow((i + 1) as u32))
            .collect::<Vec<usize>>();

        let num_queries = 10;
        let expansion_factor = 2;
        let folding_factor = 2;
        let final_codeword_size = 1;

        let fri_config = FRIConfig::<F>::new(
            (poly_degrees[poly_degrees.len() - 1] + 1).next_power_of_two(),
            expansion_factor,
            folding_factor,
            num_queries,
            final_codeword_size,
        );

        let polys = (0..num_polys)
            .map(|i| {
                UniPoly::new(
                    (0..(poly_degrees[i] + 1))
                        .map(|j| F::from((j + i) as u64))
                        .collect::<Vec<F>>(),
                )
            })
            .collect::<Vec<UniPoly<F>>>();

        let fri_prover = FRIPolyCommitProver::new(fri_config.clone());
        let mut transcript = Transcript::<F>::new(b"test_fri_batch_poly_commit");
        let proof = fri_prover.batch_prove(&polys, &mut transcript);

        let fri_verifier = FRIPolyCommitVerifier::new(fri_config);
        let mut transcript = Transcript::<F>::new(b"test_fri_batch_poly_commit");
        fri_verifier.batch_verify(&proof, &mut transcript);
    }
}
