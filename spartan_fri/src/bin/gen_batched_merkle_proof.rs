use std::env;
use std::fs::File;
use std::io::Write;

use ff::PrimeField;
use rand;
use rand::Rng;
use spartan_fri::fri::tree::{BatchedMerkleProof, CommittedMerkleTree};

use pasta_curves::group::ff::Field;

type F = pasta_curves::Fp;

pub fn save_file(name: &str, data: &[u8]) {
    let mut file = File::create(name).unwrap();
    file.write_all(data).unwrap();
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let depth = usize::from_str_radix(&args[1], 10).unwrap();
    let batch_size = usize::from_str_radix(&args[2], 10).unwrap();
    let num_leaves = 2usize.pow(depth as u32);

    let mut rng = rand::thread_rng();
    let leaves = (0..num_leaves)
        .map(|_| F::random(&mut rng))
        .collect::<Vec<F>>();

    let tree = CommittedMerkleTree::from_leaves(leaves.clone());

    let mut proofs = Vec::with_capacity(num_leaves);

    let mut opened_leaves = Vec::with_capacity(batch_size);
    for _ in 0..batch_size {
        let leaf = leaves[rng.gen_range(0..leaves.len())];
        let proof = tree.open(leaf);
        // assert!(proof.verify(), "Proof failed to verify");
        proofs.push(proof);
        opened_leaves.push(leaf);
    }

    let batched_proof: BatchedMerkleProof<pasta_curves::Fp> =
        BatchedMerkleProof::from_proofs(proofs);

    batched_proof.verify(opened_leaves.clone(), tree.root());

    let hashes_ser = batched_proof
        .hashes
        .iter()
        .map(|x| *x)
        .flatten()
        .collect::<Vec<_>>();

    let indices_ser = batched_proof
        .indices
        .iter()
        .map(|x| {
            let mut bytes = x.to_le_bytes().to_vec();
            bytes.resize(32, 0);
            bytes
        })
        .flatten()
        .collect::<Vec<_>>();

    let leaf_indices = batched_proof
        .leaf_indices
        .iter()
        .map(|x| {
            let mut bytes = x.to_le_bytes().to_vec();
            bytes.resize(32, 0);
            bytes
        })
        .flatten()
        .collect::<Vec<_>>();

    let leaves_ser = opened_leaves
        .iter()
        .map(|x| x.to_repr())
        .flatten()
        .collect::<Vec<u8>>();

    save_file("hashes.bin", &hashes_ser);
    save_file("indices.bin", &indices_ser);
    save_file("leaf_indices.bin", &leaf_indices);
    save_file("leaves.bin", &leaves_ser);
    save_file("root.bin", &tree.root());
}
