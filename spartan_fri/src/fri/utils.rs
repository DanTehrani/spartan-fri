use crate::transcript::Transcript;
use crate::FieldExt;
use ethers::abi::{encode_packed, Token};
use ethers::core::utils::keccak256;
use ethers::types::U256;

pub fn hash_two(values: &[[u8; 32]; 2]) -> [u8; 32] {
    let a = U256::from_little_endian(&values[0]);
    let b = U256::from_little_endian(&values[1]);
    let mut hashed = keccak256(&encode_packed(&[Token::Uint(a), Token::Uint(b)]).unwrap());
    hashed.reverse();
    hashed
}

fn sample_index(random_bytes: [u8; 64], size: usize) -> usize {
    let mut acc: u64 = 0;
    for b in random_bytes {
        acc = acc << 8 ^ (b as u64);
    }

    (acc % (size as u64)) as usize
}

pub fn sample_indices<F: FieldExt>(
    num_indices: usize,
    max_index: usize,
    reduced_max_index: usize,
    transcript: &mut Transcript<F>,
) -> Vec<usize> {
    println!("num_indices: {}", num_indices);
    println!("max_index: {}", max_index);
    let mut indices = Vec::with_capacity(num_indices);
    let mut counter: u32 = 0;

    while indices.len() < num_indices {
        let mut random_bytes = [0u8; 64];

        transcript.append_bytes(&counter.to_le_bytes());
        transcript.challenge_bytes(&mut random_bytes);

        let index = sample_index(random_bytes, max_index);
        if !indices.contains(&index) {
            indices.push(index);
        }
        counter += 1;
    }

    indices
}

// Have the range of the indices
pub fn reduce_indices(indices: &mut Vec<usize>, max_index: usize) {
    indices.iter_mut().for_each(|index| {
        if max_index == 1 {
            *index = 0;
        } else {
            *index = *index % max_index;
        }
    });
}
