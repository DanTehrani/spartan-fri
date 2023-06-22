use criterion::{black_box, criterion_group, criterion_main, Criterion};
use spartan_fri::fri::{FRIConfig, FRIMLPolyCommit};
use spartan_fri::spartan::prover::SpartanProver;
use spartan_fri::transcript::Transcript;
use spartan_fri::{MultilinearPCS, R1CS};

type F = pasta_curves::Fp;

fn criterion_benchmark(c: &mut Criterion) {
    let num_cons = 2usize.pow(10);
    let num_vars = num_cons;
    let num_input = 0;

    let (r1cs, witness) = R1CS::<F>::produce_synthetic_r1cs(num_cons, num_vars, num_input);

    let witness_poly_degree = num_cons;
    let num_queries = 30;
    let expansion_factor = 2;
    let folding_factor = 2;
    let final_codeword_size = 1;

    let indexer_fri_config = FRIConfig::<F>::new(
        num_cons * num_vars,
        expansion_factor,
        folding_factor,
        num_queries,
        final_codeword_size,
    );

    let fri_config = FRIConfig::<F>::new(
        witness_poly_degree,
        expansion_factor,
        folding_factor,
        num_queries,
        final_codeword_size,
    );
    let pcs_witness = FRIMLPolyCommit::<F>::new(fri_config);
    let pcs_indexer = FRIMLPolyCommit::<F>::new(indexer_fri_config);

    c.bench_function("prove", |b| {
        b.iter(|| {
            let prover = SpartanProver::<F, FRIMLPolyCommit<F>>::new(
                r1cs.clone(),
                pcs_witness.clone(),
                pcs_indexer.clone(),
            );
            let mut transcript = Transcript::new(b"bench");
            prover.prove(black_box(&witness), black_box(&mut transcript));
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
