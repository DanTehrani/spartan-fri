use super::unipoly::UniPoly;
use super::utils::sample_indices;
use super::FriEvalProof;
use crate::transcript::Transcript;
use pasta_curves::arithmetic::FieldExt;
use pasta_curves::group::ff::PrimeField;

pub struct FriVerifier<F: PrimeField<Repr = [u8; 32]> + FieldExt> {
    domain: Vec<F>,
    domain_reduced: Vec<F>,
    expansion_factor: usize, // (i.e. expansion factor) (info bits) / (total bits)
    num_colinearity_checks: usize,
}

impl<F: FieldExt<Repr = [u8; 32]>> FriVerifier<F> {
    pub fn new(max_degree: usize, num_rounds: usize) -> Self {
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
            .collect::<Vec<F>>();

        let final_domain_size = domain.len() / (2usize.pow(num_rounds as u32));
        let final_domain_generator = root_of_unity.pow(&[
            2u32.pow(32 - ((final_domain_size as f64).log2() as u32)) as u64,
            0,
            0,
            0,
        ]);
        let domain_reduced = (0..final_domain_size)
            .map(|i| final_domain_generator.pow(&[i as u64, 0, 0, 0]))
            .collect::<Vec<F>>();

        Self {
            domain,
            domain_reduced,
            expansion_factor,
            num_colinearity_checks,
        }
    }

    /*
    pub fn verify_deg(&self, proof: FriProof<F>, com: F) {
        let mut transcript =
            Transcript::<F>::new(b"Fast Reed-Solomon Interactive Oracle Proof of Proximity");

        let final_codeword = proof.reduced_codeword;
        let interpolant = UniPoly::interpolate(&self.domain_reduced, &final_codeword);

        let degree = if final_codeword.len() == 1 {
            0
        } else {
            (final_codeword.len() / self.expansion_factor)
        };

        assert_eq!(interpolant.degree(), degree);

        let domain_length = self.domain.len();

        let mut indices = sample_indices(
            self.num_colinearity_checks,
            domain_length,
            domain_length / 2usize.pow(proof.queries.len() as u32),
            &mut transcript,
        );

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
                let a_x = self.domain[a_indices[j]];
                let b_x = self.domain[b_indices[j]];
                let c_x = self.domain[c_indices[j]];

                let a_y = a.leaf * a_x;
                let b_y = b.leaf * b_x;
                let c_y = c.leaf * c_x;

                // Check that (a_x, a_y), (b_x, b_y) and (c_x, c_y) are on a straight line.
                let coeff = (a_y - b_y) * (a_x - b_x).invert().unwrap();
                let intercept = a_y - coeff * a_x;
                assert!(c_y == coeff * c_x + intercept);

                // Check Merkle proofs

                a.verify();
                b.verify();
                c.verify();

                // Check that the root is correct
                assert_eq!(
                    a.root, b.root,
                    "Roots of the two Merkle proofs are not equal"
                );

                if i == 0 && j == 0 {
                    // TODO: Uncomment this
                    // assert_eq!(c.root, com, "c.root != com");
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
     */
    pub fn verify_eval(&self, proof: &FriEvalProof<F>, transcript: &mut Transcript<F>) {
        let root = &proof.q_ld_proof.queries[0].q_openings[0].0.root;

        let final_codeword = &proof.q_ld_proof.reduced_codeword;
        let interpolant = UniPoly::interpolate(&self.domain_reduced, &final_codeword);

        // TODO Append all intermidiate Merkle roots to the transcript
        for i in 0..proof.q_ld_proof.queries.len() {
            transcript.append_fe(&proof.q_ld_proof.queries[i].q_openings[0].0.root);
            transcript.challenge_fe();
        }

        let degree = if final_codeword.len() == 1 {
            0
        } else {
            (final_codeword.len() / self.expansion_factor) - 1
        };

        assert_eq!(interpolant.degree(), degree);

        let domain_length = self.domain.len();

        let mut indices = sample_indices(
            self.num_colinearity_checks,
            domain_length,
            domain_length / 2usize.pow((proof.q_ld_proof.queries.len() - 1) as u32),
            transcript,
        );

        for (i, layer) in proof.q_ld_proof.queries.iter().enumerate() {
            assert!(
                layer.q_openings.len() == self.num_colinearity_checks,
                "Invalid number of colinearity checks"
            );

            // Halve the range of the indices
            indices = indices
                .iter()
                .map(|index| index % (domain_length / 2 >> i))
                .collect::<Vec<usize>>();

            let a_indices = indices.clone();
            let b_indices = indices
                .iter()
                .map(|index| (domain_length / 2 >> i) + index)
                .collect::<Vec<usize>>();
            let c_indices = indices.clone();

            let y = proof.y;
            let z = proof.z;
            // Colinearity checks
            for j in 0..layer.q_openings.len() {
                let a_x = self.domain[a_indices[j]];
                let b_x = self.domain[b_indices[j]];
                let c_x = self.domain[c_indices[j]];

                let q_openings = &layer.q_openings[j];
                let f_openings = &layer.f_openings[j];

                // a_y = q(a_x) = (f(a_x) - y) (a_x - z)
                //                let a_y = a.leaf;
                let a_y = (f_openings.0.leaf - y) * (a_x - z).invert().unwrap();
                let b_y = (f_openings.1.leaf - y) * (b_x - z).invert().unwrap();
                let c_y = (f_openings.2.leaf - y) * (c_x - z).invert().unwrap();

                println!(
                    "v q(x) {:?}",
                    (f_openings.0.leaf - y) * (a_x - z).invert().unwrap()
                );

                let c_x = self.domain[c_indices[j]];

                let a = &q_openings.0;
                let b = &q_openings.1;
                let c = &q_openings.2;

                assert_eq!(a_y, a.leaf);
                assert_eq!(b_y, b.leaf);
                assert_eq!(c_y, c.leaf);

                // Check that (a_x, a_y), (b_x, b_y) and (c_x, c_y) are on a straight line.
                let coeff = (a_y - b_y) * (a_x - b_x).invert().unwrap();
                let intercept = a_y - coeff * a_x;
                assert_eq!(c_y, coeff * c_x + intercept);

                // Check Merkle proofs

                f_openings.0.verify();
                f_openings.1.verify();
                f_openings.2.verify();

                assert_eq!(f_openings.0.root, proof.f_comm);
                assert_eq!(f_openings.1.root, proof.f_comm);
                assert_eq!(f_openings.2.root, proof.f_comm);

                a.verify();
                b.verify();
                c.verify();

                // Check that the root is correct
                assert_eq!(
                    a.root, b.root,
                    "Roots of the two Merkle proofs are not equal"
                );

                if i == 0 && j == 0 {
                    assert_eq!(c.root, *root, "c.root != root");
                } else if j > 0 {
                    assert_eq!(
                        layer.q_openings[j - 1].0.root,
                        c.root,
                        "a_prev.root != c.root"
                    );

                    assert_eq!(
                        layer.q_openings[j - 1].1.root,
                        c.root,
                        "b_prev.root != c.root"
                    );
                }
            }
        }
    }
}
