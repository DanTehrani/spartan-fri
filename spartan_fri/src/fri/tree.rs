use super::utils::hash_two;
use crate::FieldExt;
use serde::{Deserialize, Serialize};

#[derive(Clone)]
pub struct CommittedMerkleTree<F: FieldExt> {
    pub layers: Vec<Vec<F>>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct MerkleProof<F: FieldExt> {
    pub root: F,
    pub leaf: F,
    pub siblings: Vec<F>,
    pub leaf_index: usize,
}

impl<F: FieldExt> MerkleProof<F> {
    pub fn verify(&self) -> bool {
        let mut current_hash = self.leaf;
        let mut index = self.leaf_index;
        for sibling in &self.siblings {
            current_hash = if index & 1 == 0 {
                hash_two(&[current_hash, *sibling])
            } else {
                hash_two(&[*sibling, current_hash])
            };

            index >>= 1;
        }

        current_hash == self.root
    }
}

impl<F: FieldExt> CommittedMerkleTree<F> {
    pub fn from_leaves(leaves: Vec<F>) -> Self {
        let n = leaves.len();
        assert!(n.is_power_of_two());

        let num_layers = (n as f64).log2() as usize + 1;
        let mut layers = Vec::with_capacity(num_layers);

        // Add a dummy leaf if the number of leaves is odd.
        let mut leaves = leaves.to_vec();
        if n % 2 == 1 {
            leaves.push(F::ZERO);
        }

        // The first layer is the leaves.
        layers.push(leaves.clone());

        while leaves.len() != 1 {
            let mut layer = vec![];
            for i in (0..leaves.len()).step_by(2) {
                let left = leaves[i];
                let right = leaves[i + 1];
                let parent = hash_two(&[left, right]);
                layer.push(parent);
            }
            layers.push(layer.clone());
            leaves = layer;
        }

        assert_eq!(layers.len(), num_layers);
        assert_eq!(layers[num_layers - 1].len(), 1);

        Self { layers }
    }

    pub fn root(&self) -> F {
        self.layers[self.layers.len() - 1][0]
    }

    pub fn leaves(&self) -> Vec<F> {
        self.layers[0].clone()
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
            root: self.layers.last().unwrap()[0],
            leaf: self.layers[0][leaf_index],
            siblings,
            leaf_index,
        }
    }

    pub fn open(&self, leaf: F) -> MerkleProof<F> {
        // Find the index of the leaf
        let leaf_index = self.layers[0].iter().position(|&x| x == leaf).unwrap();
        self.open_index(leaf_index)
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

        let tree = CommittedMerkleTree::from_leaves(leaves.clone());

        for i in 0..leaves.len() {
            let proof = tree.open(leaves[i]);
            assert!(proof.verify(), "Proof failed to verify");
        }

        // Should assert invalid opening proof
        let mut proof = tree.open(Fp::from(5));
        proof.siblings[0] = proof.siblings[0] + F::ONE;
        assert!(!proof.verify(), "Verifycation should fail");
    }
}
