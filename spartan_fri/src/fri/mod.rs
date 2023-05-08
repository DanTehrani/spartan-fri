mod fft;
mod fri_prover;
mod fri_verifier;
mod tree;
mod unipoly;
mod utils;

use pasta_curves::arithmetic::FieldExt;
use tree::MerkleProof;

pub use fri_prover::FriProver;
pub use fri_verifier::FriVerifier;
pub use merlin::Transcript;
pub use unipoly::UniPoly;

use crate::PolyCommitment;

#[derive(Clone, Debug)]
pub struct LayerProof<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub q_openings: Vec<(MerkleProof<F>, MerkleProof<F>, MerkleProof<F>)>,
    pub f_openings: Vec<(MerkleProof<F>, MerkleProof<F>, MerkleProof<F>)>,
}

#[derive(Clone)]
pub struct FriProof<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub reduced_codeword: Vec<F>,
    pub queries: Vec<LayerProof<F>>,
}

pub struct FriEvalProof<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub z: F,
    pub y: F,
    pub q_ld_proof: FriProof<F>,
    pub f_comm: F,
}

use std::marker::PhantomData;

pub struct FRIMultilinearPCS<F> {
    _marker: PhantomData<F>,
}

impl<F: FieldExt> PolyCommitment<F> for FRIMultilinearPCS<F> {
    type Commitment = F; // tmp
    type Opening = F; // tmp
    fn new() -> Self {
        todo!()
    }

    fn commit(&self) -> Self::Commitment {
        todo!()
    }

    fn open(&self, point: &[F]) -> Self::Opening {
        todo!()
    }
}

#[cfg(test)]
mod tests {

    use super::unipoly::UniPoly;
    use super::*;
    use crate::transcript::Transcript;
    use pasta_curves::Fp;

    type F = Fp;

    /*
    #[test]
    fn test_fri() {
        let poly_degree = 2u32.pow(13);

        let mut coeffs = vec![];
        for i in 0..(poly_degree + 1) {
            coeffs.push(F::from(i as u64));
        }

        let poly = UniPoly::new(coeffs);
        let prover = FriProver::<F>::new(poly.degree());

        let mut transcript = Transcript::new(b"test_fri");
        let mut proof = prover.prove_degree(&poly, &mut transcript);

        // C of the first round is the polynomial we're committing to.
        let poly_commitment = proof.queries[0].openings[0].2.root;

        let num_rounds = proof.queries.len();
        let verifier = FriVerifier::<F>::new(poly.degree(), num_rounds);

        verifier.verify_deg(proof, poly_commitment);
    }
     */

    #[test]
    fn test_fri_eval() {
        let poly_degree = 2u32.pow(5u32);

        let mut coeffs = vec![];
        for i in 0..(poly_degree + 1) {
            coeffs.push(F::from(i as u64));
        }

        let poly = UniPoly::new(coeffs);
        let z = F::from(33u64);

        let prover = FriProver::<Fp>::new(poly.degree());

        let mut prover_transcript = Transcript::new(b"test_fri_eval");
        let fri_eval_proof = prover.prove_eval(&poly, z, &mut prover_transcript);

        let num_rounds = fri_eval_proof.q_ld_proof.queries.len();
        let verifier = FriVerifier::<Fp>::new(poly.degree(), num_rounds);

        let mut verifier_transcript = Transcript::new(b"test_fri_eval");
        verifier.verify_eval(&fri_eval_proof, &mut verifier_transcript);
    }
}
