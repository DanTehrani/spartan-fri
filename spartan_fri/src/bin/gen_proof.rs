use std::env;
use std::fs::File;
use std::io::Write;

use bincode::serialize;
use spartan_fri::fri::{FRIConfig, FRIMLPolyCommitProver};
use spartan_fri::spartan::polynomial::ml_poly::MlPoly;
use spartan_fri::transcript::Transcript;
use spartan_fri::MultilinearPCS;

type F = pasta_curves::Fp;

fn gen_config(num_vars: usize) -> FRIConfig<F> {
    let n = 2usize.pow(num_vars as u32);

    let poly_degree = n;
    let num_queries = 30;
    let expansion_factor = 2;
    let folding_factor = 2;
    let final_codeword_size = 1;

    let fri_config = FRIConfig::<F>::new(
        poly_degree,
        expansion_factor,
        folding_factor,
        num_queries,
        final_codeword_size,
    );
    fri_config
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let num_vars = usize::from_str_radix(&args[1], 10).unwrap();
    println!("num_vars: {}", num_vars);
    let fri_config: FRIConfig<F> = gen_config(num_vars);

    let prover: FRIMLPolyCommitProver<pasta_curves::Fp> =
        FRIMLPolyCommitProver::new(fri_config.clone());
    let mut transcript = Transcript::new(b"test");

    let n = 2usize.pow(num_vars as u32);
    let evals = (0..n).map(|i| F::from(i as u64)).collect::<Vec<F>>();
    let mut ml_poly = MlPoly::new(evals);
    ml_poly.compute_coeffs();

    let eval_at = (0..num_vars).map(|i| F::from(i as u64)).collect::<Vec<F>>();
    let proof = prover.prove_eval(&ml_poly, &eval_at, &mut transcript);
    let proof_ser = serialize(&proof).unwrap();

    let mut proof_file = File::create("proof.bin").unwrap();
    proof_file.write_all(&proof_ser).unwrap();
}
