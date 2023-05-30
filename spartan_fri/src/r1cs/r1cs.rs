use crate::FieldExt;

#[derive(Clone)]
pub struct Matrix<F>(Vec<(usize, usize, F)>)
where
    F: FieldExt;

impl<F> Matrix<F>
where
    F: FieldExt,
{
    pub fn mul_vector(&self, num_rows: usize, vec: &[F]) -> Vec<F> {
        let mut result = vec![F::ZERO; num_rows];
        for i in 0..self.0.len() {
            let row = self.0[i].0;
            let col = self.0[i].1;
            let val = self.0[i].2;
            result[row] += val * vec[col];
        }
        result
    }
}

#[derive(Clone)]
pub struct R1CS<F>
where
    F: FieldExt,
{
    pub A: Matrix<F>,
    pub B: Matrix<F>,
    pub C: Matrix<F>,
    pub public_input: Vec<F>,
    pub num_cons: usize,
    pub num_vars: usize,
    pub num_input: usize,
}

impl<F> R1CS<F>
where
    F: FieldExt,
{
    pub fn hadamard_prod(a: &[F], b: &[F]) -> Vec<F> {
        assert_eq!(a.len(), b.len());
        let mut result = vec![F::ZERO; a.len()];
        for i in 0..a.len() {
            result[i] = a[i] * b[i];
        }
        result
    }

    pub fn produce_synthetic_r1cs(
        num_cons: usize,
        num_vars: usize,
        num_input: usize,
    ) -> (Self, Vec<F>) {
        //        assert_eq!(num_cons, num_vars);
        let mut public_input = Vec::with_capacity(num_input);
        let mut witness = Vec::with_capacity(num_vars);

        for i in 0..num_input {
            public_input.push(F::from((i + 1) as u64));
        }

        for i in 0..num_vars {
            witness.push(F::from((i + 1) as u64));
        }

        let z: Vec<F> = vec![public_input.clone(), witness.clone()].concat();

        let mut A: Vec<(usize, usize, F)> = vec![];
        let mut B: Vec<(usize, usize, F)> = vec![];
        let mut C: Vec<(usize, usize, F)> = vec![];

        for i in 0..num_cons {
            let A_col = i % num_vars;
            let B_col = (i + 1) % num_vars;
            let C_col = (i + 2) % num_vars;

            // For the i'th constraint,
            // add the value 1 at the (i % num_vars)th column of A, B.
            // Compute the corresponding C_column value so that A_i * B_i = C_i
            // we apply multiplication since the Hadamard product is computed for Az ・ Bz,

            // We only _enable_ a single variable in each constraint.
            A.push((i, A_col, F::ONE));
            B.push((i, B_col, F::ONE));
            C.push((i, C_col, (z[A_col] * z[B_col]) * z[C_col].invert().unwrap()));
        }

        (
            Self {
                A: Matrix(A),
                B: Matrix(B),
                C: Matrix(C),
                public_input,
                num_cons,
                num_vars,
                num_input,
            },
            witness,
        )
    }

    pub fn is_sat(&self, witness: &Vec<F>, public_input: &Vec<F>) -> bool {
        let mut z = Vec::with_capacity(witness.len() + public_input.len() + 1);
        z.extend(public_input);
        z.extend(witness);

        let Az = self.A.mul_vector(self.num_cons, &z);
        let Bz = self.B.mul_vector(self.num_cons, &z);
        let Cz = self.C.mul_vector(self.num_cons, &z);

        Self::hadamard_prod(&Az, &Bz) == Cz
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pasta_curves::Fp;

    #[test]
    fn test_r1cs() {
        let num_cons = 20;
        let num_vars = 10;
        let num_input = 5;

        let (r1cs, mut witness) = R1CS::<Fp>::produce_synthetic_r1cs(num_cons, num_vars, num_input);

        assert_eq!(witness.len(), num_vars);
        assert_eq!(r1cs.public_input.len(), num_input);

        assert!(r1cs.is_sat(&witness, &r1cs.public_input));

        // Should assert if the witness is invalid
        witness[0] = witness[0] + Fp::one();
        assert!(r1cs.is_sat(&r1cs.public_input, &witness) == false);
        witness[0] = witness[0] - Fp::one();

        // Should assert if the public input is invalid
        let mut public_input = r1cs.public_input.clone();
        public_input[0] = public_input[0] + Fp::one();
        assert!(r1cs.is_sat(&witness, &public_input) == false);
    }
}
