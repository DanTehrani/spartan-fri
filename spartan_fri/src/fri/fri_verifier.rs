use std::marker::PhantomData;

use super::unipoly::UniPoly;
use super::utils::sample_indices;
use super::FriProof;
use crate::transcript::Transcript;
use pasta_curves::arithmetic::FieldExt;
use pasta_curves::group::ff::PrimeField;

pub struct FriVerifier<F: PrimeField<Repr = [u8; 32]> + FieldExt> {
    domain: Vec<F>,
    expansion_factor: usize, // (i.e. expansion factor) (info bits) / (total bits)
    num_colinearity_checks: usize,
}

impl<F: FieldExt<Repr = [u8; 32]>> FriVerifier<F> {
    pub fn new(max_degree: usize) -> Self {
        // TODO: Allow arbitrary degree
        assert!(max_degree.is_power_of_two());

        // Are these params OK?
        let expansion_factor = 2;
        let num_colinearity_checks = 2;

        let root_of_unity = F::root_of_unity();
        let root_of_unity_log_2 = F::S;

        let domain_order = (max_degree * expansion_factor).next_power_of_two();

        // Generator for the subgroup with order _subgroup_order_ in the field
        let domain_generator = root_of_unity.pow(&[
            2u32.pow(32 - ((domain_order as f64).log2() as u32)) as u64,
            0,
            0,
            0,
        ]);

        let domain = (0..domain_order)
            .map(|i| domain_generator.pow(&[i as u64, 0, 0, 0]))
            .collect();

        Self {
            domain,
            expansion_factor,
            num_colinearity_checks,
        }
    }

    pub fn verify(&self, proof: FriProof<F>, com: F) {
        let mut transcript =
            Transcript::<F>::new(b"Fast Reed-Solomon Interactive Oracle Proof of Proximity");

        let final_codeword = proof.reduced_codeword;

        // TODO: Precompute this
        let mut domain_reduced = self.domain.clone();
        for _ in 0..proof.queries.len() {
            let mut domain_unique = vec![];
            domain_reduced.iter().map(|x| x.square()).for_each(|x| {
                if !domain_unique.contains(&x) {
                    domain_unique.push(x);
                }
            });
            domain_reduced = domain_unique;
        }

        let interpolant = UniPoly::interpolate(&domain_reduced, &final_codeword);

        let degree = if final_codeword.len() == 1 {
            0
        } else {
            final_codeword.len() / self.expansion_factor
        };

        assert!(interpolant.degree() == degree);

        let domain_length = self.domain.len();

        let mut indices = sample_indices(
            self.num_colinearity_checks,
            domain_length,
            domain_length / 2usize.pow(proof.queries.len() as u32),
            &mut transcript,
        );

        let mut hash_count = 0;
        for (i, layer) in proof.queries.iter().enumerate() {
            assert!(
                layer.openings.len() == self.num_colinearity_checks,
                "Invalid number of colinearity checks"
            );

            // Halve the range of the indices
            indices = indices
                .iter()
                .map(|index| index % (domain_length / 2 >> (i + 1)))
                .collect::<Vec<usize>>();

            let a_indices = indices.clone();
            let b_indices = indices
                .iter()
                .map(|index| (domain_length / 2 >> (i + 1)) + index)
                .collect::<Vec<usize>>();
            let c_indices = indices.clone();

            // Colinearity checks
            for (j, (a, b, c)) in layer.openings.iter().enumerate() {
                let a_y = a.leaf;
                let b_y = b.leaf;
                let c_y = c.leaf;

                let a_x = self.domain[a_indices[j]];
                let b_x = self.domain[b_indices[j]];
                let c_x = self.domain[c_indices[j]];

                // Check that (a_x, a_y), (b_x, b_y) and (c_x, c_y) are on a straight line.
                let coeff = (a_y - b_y) * (a_x - b_x).invert().unwrap();
                let intercept = a_y - coeff * a_x;
                assert!(c_y == coeff * c_x + intercept);

                // Check Merkle proofs

                a.verify();
                b.verify();
                c.verify();
                hash_count += a.siblings.len() + b.siblings.len() + c.siblings.len();

                // Check that the root is correct
                assert_eq!(
                    a.root, b.root,
                    "Roots of the two Merkle proofs are not equal"
                );

                if i == 0 && j == 0 {
                    assert_eq!(c.root, com, "c.root != com");
                } else if j > 0 {
                    assert_eq!(
                        layer.openings[j - 1].0.root,
                        c.root,
                        "a_prev.root != c.root"
                    );

                    assert_eq!(
                        layer.openings[j - 1].1.root,
                        c.root,
                        "b_prev.root != c.root"
                    );
                }
            }
        }
    }
}
