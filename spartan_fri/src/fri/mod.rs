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

#[derive(Debug)]
pub struct LayerProof<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub openings: Vec<(MerkleProof<F>, MerkleProof<F>, MerkleProof<F>)>,
}

pub struct FriProof<F>
where
    F: FieldExt<Repr = [u8; 32]>,
{
    pub reduced_codeword: Vec<F>,
    pub queries: Vec<LayerProof<F>>,
}

#[cfg(test)]
mod tests {
    use super::unipoly::UniPoly;
    use super::*;
    use crate::transcript::Transcript;
    use pasta_curves::Fp;

    #[test]
    fn test_prove() {
        let poly_degree = 2u32.pow(4u32);

        let mut coeffs = vec![];
        for i in 0..(poly_degree + 1) {
            coeffs.push(Fp::from(i as u64));
        }

        let poly = UniPoly::new(coeffs);
        let prover = FriProver::<Fp>::new(poly.degree());

        let mut transcript = Transcript::new(b"test_fri");
        let proof = prover.prove_degree(&poly, &mut transcript);

        // C of the first round is the polynomial we're committing to.
        let poly_commitment = proof.queries[0].openings[0].2.root;

        let verifier = FriVerifier::<Fp>::new(poly.degree());

        verifier.verify(proof, poly_commitment);
    }
}
