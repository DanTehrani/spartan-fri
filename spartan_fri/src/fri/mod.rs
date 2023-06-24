mod fft;
mod fri;
mod sparse_fri;
mod sparse_unipoly;
pub mod tree;
mod unipoly;
mod utils;

use std::marker::PhantomData;

use crate::{FieldExt, MultilinearPCS, MultilinearPCSOpening};
use serde::{Deserialize, Serialize};
use tree::MerkleProof;

pub use fri::FRIMLPolyCommit;
pub use sparse_fri::SparseFRIMLPolyCommit;
pub use unipoly::UniPoly;

#[derive(Clone, Serialize, Deserialize)]
pub struct FirstRoundBoundedPolyOpenings<F>
where
    F: FieldExt,
{
    pub opening_at_s0: MerkleProof<F>,
    pub opening_at_s1: MerkleProof<F>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct QueryFirstRound<F>
where
    F: FieldExt,
{
    // Openings of f_0(\s_0), f_1(\s_0), ... f_{u-1}(s_0)
    // and f_0(\s_1), f_1(\s_1), ... f_{u-1}(s_1)
    pub bounded_poly_openings: Vec<FirstRoundBoundedPolyOpenings<F>>,
    // Openings of g_1(t)
    pub g_1_opening_at_t: MerkleProof<F>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct RoundQueries<F>
where
    F: FieldExt,
{
    pub queries: Vec<Query<F>>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Query<F>
where
    F: FieldExt,
{
    pub opening_at_s0: MerkleProof<F>,
    pub opening_at_s1: MerkleProof<F>,
    pub opening_at_t: MerkleProof<F>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct MLPolyEvalProof<F>
where
    F: FieldExt,
{
    pub queries_first_round: Vec<QueryFirstRound<F>>,
    pub queries: Vec<RoundQueries<F>>,
    pub reduced_codeword: Vec<F>,
    pub f_i_evals_beta: Vec<(F, F)>,
    pub f_i_evals_beta_squared: Vec<F>,
    pub y: F,
    pub x: Vec<F>,
}

impl<F: FieldExt> MultilinearPCSOpening<F> for MLPolyEvalProof<F> {
    fn x(&self) -> Vec<F> {
        self.x.clone()
    }

    fn y(&self) -> F {
        self.y
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct FRIConfig<F>
where
    F: FieldExt,
{
    pub expansion_factor: usize,
    pub folding_factor: usize,
    pub num_queries: usize,
    pub L: Vec<Vec<F>>,
}

impl<F> FRIConfig<F>
where
    F: FieldExt,
{
    pub fn new(
        poly_degree: usize,
        expansion_factor: usize,
        folding_factor: usize,
        num_queries: usize,
        final_codeword_size: usize,
    ) -> Self {
        debug_assert!(folding_factor.is_power_of_two());

        let root_of_unity = F::ROOT_OF_UNITY;
        let mut domain_order = (poly_degree * expansion_factor).next_power_of_two();
        let mut L = vec![];

        let mut codeword_size = poly_degree * expansion_factor;
        println!("final_codeword_size: {}", final_codeword_size);
        while codeword_size > final_codeword_size {
            println!("codeword_size: {}", codeword_size);
            if L.len() == 0 {
                // Generator for the subgroup with order _subgroup_order_ in the field
                let domain_generator = root_of_unity.pow(&[
                    2u32.pow(32 - ((domain_order as f64).log2() as u32)) as u64,
                    0,
                    0,
                    0,
                ]);

                let mut L_i = Vec::with_capacity(domain_order);
                L_i.push(F::ONE);
                for i in 0..(domain_order - 1) {
                    L_i.push(L_i[i] * domain_generator);
                }

                codeword_size /= folding_factor;
                domain_order /= folding_factor;
                L.push(L_i);
            } else {
                let L_i = (0..domain_order)
                    .map(|i| L[L.len() - 1][i * folding_factor])
                    .collect::<Vec<F>>();
                codeword_size /= folding_factor;
                domain_order /= folding_factor;
                L.push(L_i);
            }
        }

        Self {
            expansion_factor,
            folding_factor,
            num_queries,
            L,
        }
    }
}

impl<F> FRIConfig<F>
where
    F: FieldExt,
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

#[derive(Serialize, Deserialize, Clone)]
pub struct SparseFRIConfig<F>
where
    F: FieldExt,
{
    pub expansion_factor: usize,
    pub folding_factor: usize,
    pub num_queries: usize,
    pub poly_degree: usize,
    pub final_codeword_size: usize,
    pub _marker: PhantomData<F>,
}

impl<F> SparseFRIConfig<F>
where
    F: FieldExt,
{
    pub fn new(
        poly_degree: usize,
        expansion_factor: usize,
        folding_factor: usize,
        num_queries: usize,
        final_codeword_size: usize,
    ) -> Self {
        debug_assert!(folding_factor.is_power_of_two());

        Self {
            expansion_factor,
            folding_factor,
            num_queries,
            poly_degree,
            final_codeword_size,
            _marker: PhantomData,
        }
    }
}

impl<F> SparseFRIConfig<F>
where
    F: FieldExt,
{
    pub fn domain_order_in_round(&self, round: usize) -> usize {
        println!("round: {}", round);
        println!("poly_Degree {}", self.poly_degree);
        println!(
            "poly_Degree.next_power_of_two {}",
            self.poly_degree.next_power_of_two()
        );
        ((self.poly_degree + 1).next_power_of_two() * self.expansion_factor)
            / self.folding_factor.pow(round as u32)
    }

    pub fn num_rounds(&self) -> usize {
        let mut num_rounds = 0;
        let mut codeword_size = self.domain_order_in_round(0);
        while codeword_size > self.final_codeword_size {
            num_rounds += 1;
            codeword_size /= self.folding_factor;
        }
        num_rounds
    }
    pub fn domain_generator(&self, round: usize) -> F {
        let domain_order = self.domain_order_in_round(round);
        let root_of_unity = F::ROOT_OF_UNITY;
        let domain_generator = root_of_unity.pow(&[
            2u32.pow(32 - ((domain_order as f64).log2() as u32)) as u64,
            0,
            0,
            0,
        ]);
        domain_generator
    }
}

#[cfg(test)]
mod tests {
    use ark_std::{end_timer, start_timer};

    use super::*;
    use crate::spartan::polynomial::ml_poly::MlPoly;
    use crate::spartan::polynomial::sparse_ml_poly::SparseMLPoly;
    use crate::transcript::Transcript;
    use crate::MultilinearPCS;
    type F = pasta_curves::Fp;

    #[test]
    fn test_ml_poly_prove() {
        let m = 4;
        let n = 2usize.pow(m);

        let poly_degree = n;
        let num_queries = 31;
        let expansion_factor = 2;
        let folding_factor = 2;
        let final_codeword_size = 4;

        let evals = (0..n).map(|i| F::from(i as u64)).collect::<Vec<F>>();
        let mut ml_poly = MlPoly::new(evals);

        let compute_coeffs_timer = start_timer!(|| "compute coeffs");
        ml_poly.compute_coeffs();
        end_timer!(compute_coeffs_timer);

        let fri_config = FRIConfig::<F>::new(
            poly_degree,
            expansion_factor,
            folding_factor,
            num_queries,
            final_codeword_size,
        );

        let eval_at = (0..m).map(|i| F::from(i as u64)).collect::<Vec<F>>();

        let mut prover_transcript = Transcript::<F>::new(b"test");
        let fri = FRIMLPolyCommit::new(fri_config.clone());
        let prove_timer = start_timer!(|| "prove");
        let comm = fri.commit(&ml_poly);
        let proof = fri.prove_eval(&ml_poly, &eval_at, &mut prover_transcript);
        end_timer!(prove_timer);

        let mut verifier_transcript = Transcript::<F>::new(b"test");
        // fri.verify(&proof, &comm, &mut verifier_transcript);
    }
}
