use std::marker::PhantomData;

use super::utils::hash_two;
use crate::FieldExt;
use ethers::types::U256;
use serde::{Deserialize, Serialize};

#[derive(Clone)]
pub struct CommittedMerkleTree<F> {
    pub layers: Vec<Vec<[u8; 32]>>,
    pub leaves: Vec<F>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct MerkleProof<F: FieldExt> {
    pub hashes: Vec<[u8; 32]>,
    pub indices: Vec<usize>,
    pub leaf_index: usize,
    pub leaf: [u8; 32],
    pub root: [u8; 32],
    _marker: PhantomData<F>,
}

impl<F: FieldExt> MerkleProof<F> {
    pub fn verify(&self) -> bool {
        let mut current_hash = self.leaf;
        let mut leaf_index = self.leaf_index;
        for index in &self.indices {
            current_hash = if leaf_index & 1 == 0 {
                hash_two(&[current_hash, self.hashes[*index]])
            } else {
                hash_two(&[self.hashes[*index], current_hash])
            };

            leaf_index >>= 1;
        }

        current_hash == self.root
    }

    pub fn leaf_fe(&self) -> F {
        F::from_repr(self.leaf).unwrap()
    }
}

impl<F: FieldExt> CommittedMerkleTree<F> {
    pub fn from_leaves(leaves: Vec<F>) -> Self {
        let n = leaves.len();
        assert!(n.is_power_of_two());

        let num_layers = (n as f64).log2() as usize + 1;
        let mut layers = Vec::with_capacity(num_layers);

        // Add a dummy leaf if the number of leaves is odd.
        let mut first_layer = leaves
            .iter()
            .map(|x| x.to_repr())
            .collect::<Vec<[u8; 32]>>();
        if n % 2 == 1 {
            first_layer.push(F::ZERO.to_repr());
        }
        layers.push(first_layer);

        while layers.len() < num_layers {
            let prev_layer = &layers[layers.len() - 1];
            let mut layer = vec![];
            for i in (0..layers[layers.len() - 1].len()).step_by(2) {
                let left = prev_layer[i];
                let right = prev_layer[i + 1];
                let parent = hash_two(&[left, right]);
                layer.push(parent);
            }
            layers.push(layer);
        }

        assert_eq!(layers.len(), num_layers);
        assert_eq!(layers[num_layers - 1].len(), 1);

        Self { layers, leaves }
    }

    pub fn root(&self) -> [u8; 32] {
        self.layers[self.layers.len() - 1][0]
    }

    pub fn leaves(&self) -> Vec<F> {
        self.leaves.clone()
    }

    pub fn depth(&self) -> usize {
        self.layers.len()
    }

    pub fn open_index(&self, leaf_index: usize) -> MerkleProof<F> {
        let m = self.layers.len();
        let mut sibling_indices = Vec::with_capacity(m);

        // Find the indices of the siblings through all layers
        for i in 0..(self.layers.len() - 1) {
            if i == 0 {
                // For the first layer, the sibling is the next/previous leaf
                sibling_indices.push(if leaf_index & 1 == 0 {
                    // leaf_index is even, so take the next leaf
                    leaf_index + 1
                } else {
                    // leaf_index is odd, so take the previous leaf
                    leaf_index - 1
                });
            } else {
                let current_index = sibling_indices[i - 1] / 2;
                sibling_indices.push(if current_index & 1 == 0 {
                    if current_index == self.layers[i].len() - 1 {
                        current_index
                    } else {
                        current_index + 1
                    }
                } else {
                    current_index - 1
                });
            }
        }

        let mut siblings = Vec::with_capacity(m);
        // Get the sibling values from the indices
        for (i, index) in sibling_indices.iter().enumerate() {
            siblings.push(self.layers[i][*index]);
        }

        MerkleProof {
            indices: (0..siblings.len()).map(|i| i).collect(),
            hashes: siblings,
            leaf_index,
            leaf: self.layers[0][leaf_index],
            root: self.root(),
            _marker: PhantomData,
        }
    }

    pub fn open(&self, leaf: F) -> MerkleProof<F> {
        // Find the index of the leaf
        let leaf_index = self.layers[0]
            .iter()
            .position(|&x| x == leaf.to_repr())
            .unwrap();
        self.open_index(leaf_index)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct BatchedMerkleProof<F: FieldExt> {
    pub depth: usize,
    pub hashes: Vec<[u8; 32]>,
    pub indices: Vec<usize>,
    pub leaf_indices: Vec<usize>,
    _marker: PhantomData<F>,
}

impl<F: FieldExt> BatchedMerkleProof<F> {
    pub fn from_proofs(proofs: Vec<MerkleProof<F>>) -> Self {
        let num_proofs = proofs.len();
        let depth = proofs[0].hashes.len();

        let hashes: Vec<[u8; 32]> = proofs.iter().map(|p| p.hashes.clone()).flatten().collect();
        let mut unique_hashes = vec![];
        for hash in hashes {
            if !unique_hashes.contains(&hash) {
                unique_hashes.push(hash);
            }
        }

        let mut indices = Vec::with_capacity(num_proofs);
        let leaf_indices = proofs.iter().map(|p| p.leaf_index).collect();
        for proof in proofs {
            let mut proof_indices = Vec::with_capacity(depth);
            for index in proof.indices {
                let sibling = proof.hashes[index];
                let sibling_index = unique_hashes.iter().position(|&x| x == sibling).unwrap();
                proof_indices.push(sibling_index);
            }
            indices.push(proof_indices);
        }

        Self {
            depth,
            hashes: unique_hashes,
            indices: indices.into_iter().flatten().collect(),
            leaf_indices,
            _marker: PhantomData,
        }
    }

    pub fn verify(&self, leaves: Vec<F>, root: [u8; 32]) {
        let num_proofs = self.indices.len() / self.depth;
        assert_eq!(leaves.len(), num_proofs);

        for i in 0..num_proofs {
            let indices = self.indices[(i * self.depth)..(i * self.depth) + self.depth].to_vec();
            assert_eq!(indices.len(), self.depth);
            let mut current_hash = leaves[i].to_repr();
            let mut leaf_index = self.leaf_indices[i];
            for index in indices {
                // println!("lhs: {:?}", U256::from_little_endian(&current_hash));
                current_hash = if leaf_index & 1 == 0 {
                    hash_two(&[current_hash, self.hashes[index]])
                } else {
                    hash_two(&[self.hashes[index], current_hash])
                };
                /*
                println!("leaf_index: {:?}", leaf_index);
                println!("sibling_index: {:?}", index);
                println!(
                    "sibling: {:?}",
                    U256::from_little_endian(&self.hashes[index])
                );
                println!("hash: {:?}", U256::from_little_endian(&current_hash));
                 */
                leaf_index >>= 1;
            }
            assert_eq!(current_hash, root);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff::Field;
    use pasta_curves::Fp;

    type F = Fp;

    #[test]
    fn test_tree() {
        let num_leaves = 2usize.pow(5);

        let leaves = (0..num_leaves)
            .map(|x| F::from(x as u64))
            .collect::<Vec<F>>();

        let tree = CommittedMerkleTree::<F>::from_leaves(leaves.clone());
        let root = tree.root();
        let open_leaf = leaves[3];

        for i in 0..leaves.len() {
            let proof = tree.open(leaves[i]);
            assert!(proof.verify(), "Proof failed to verify");
        }

        // Should assert invalid opening proof
        let mut proof = tree.open(open_leaf);
        proof.indices[0] = proof.indices[0] + 1;
        assert!(!proof.verify(), "Verifycation should fail");
    }
}
