use super::utils::hash_two;
use ff::Field;
use pasta_curves::arithmetic::FieldExt;

pub struct MerkleTree<F: FieldExt<Repr = [u8; 32]>> {
    pub layers: Vec<Vec<F>>, // to root
}

#[derive(Clone, Debug)]
pub struct MerkleProof<F: FieldExt<Repr = [u8; 32]>> {
    pub root: F,
    pub leaf: F,
    pub siblings: Vec<F>,
    pub leaf_index: usize,
    pub tree_depth: usize,
}

impl<F: FieldExt<Repr = [u8; 32]>> MerkleProof<F> {
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

impl<F: FieldExt<Repr = [u8; 32]>> MerkleTree<F> {
    pub fn new() -> Self {
        Self { layers: vec![] }
    }

    fn compute_layers(&mut self, leaves: Vec<F>) -> Vec<F> {
        self.layers.push(leaves.clone());

        if leaves.len() == 1 {
            return leaves;
        }

        let mut parent_nodes = vec![];
        for i in (0..leaves.len()).step_by(2) {
            if leaves.len() - 1 == i {
                // Sibling.
                let parent = hash_two(&[leaves[i], leaves[i]]);
                parent_nodes.push(parent);
            } else {
                let left = leaves[i];
                let right = leaves[i + 1];
                let parent = hash_two(&[left, right]);
                parent_nodes.push(parent);
            }
        }
        self.compute_layers(parent_nodes)
    }

    pub fn commit(&mut self, leaves: &[F]) -> F {
        let n = leaves.len();
        assert!(n.is_power_of_two());

        // Add a dummy leaf if the number of leaves is odd.
        let mut leaves = leaves.to_vec();
        if n % 2 == 1 {
            leaves.push(F::zero());
        }

        // The first layer is the leaves.
        self.layers.push(leaves.clone());

        while leaves.len() != 1 {
            let mut layer = vec![];
            for i in (0..leaves.len()).step_by(2) {
                let left = leaves[i];
                let right = leaves[i + 1];
                let parent = hash_two(&[left, right]);
                layer.push(parent);
            }
            self.layers.push(layer.clone());
            leaves = layer;
        }

        leaves[0]
    }

    pub fn open(&self, leaf: F) -> MerkleProof<F> {
        let m = self.layers.len();
        // Find the index of the leaf
        let leaf_index = self.layers[0].iter().position(|&x| x == leaf).unwrap();

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
            leaf,
            siblings,
            leaf_index,
            tree_depth: self.layers.len() - 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pasta_curves::Fp;

    type F = Fp;

    #[test]
    fn test_tree() {
        let num_leaves = 2usize.pow(5);
        let mut tree = MerkleTree::new();

        let leaves = (0..num_leaves)
            .map(|x| F::from(x as u64))
            .collect::<Vec<F>>();

        tree.commit(&leaves);

        for i in 0..leaves.len() {
            let proof = tree.open(leaves[i]);
            assert!(proof.verify(), "Proof failed to verify");
        }

        // Should assert invalid opening proof
        let mut proof = tree.open(Fp::from(5));
        proof.siblings[0] = proof.siblings[0] + F::one();
        assert!(!proof.verify(), "Verifycation should fail");
    }
}
