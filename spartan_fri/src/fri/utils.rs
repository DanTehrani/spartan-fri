use crate::transcript::Transcript;
use pasta_curves::arithmetic::FieldExt;
use sha3::Digest;
use sha3::Keccak256;

pub fn hash_two<F>(values: &[F; 2]) -> F
where
    F: FieldExt<Repr = [u8; 32]>,
{
    /*
    let mut bytes = vec![];
    bytes.extend_from_slice(&values[0].to_repr());
    bytes.resize(32, 0);
    bytes.extend_from_slice(&values[1].to_repr());
    bytes.resize(32, 0);
     */

    let mut hasher = Keccak256::new();
    hasher.update(values[0].to_repr());
    hasher.update(values[1].to_repr());
    let result = hasher.finalize();
    let bytes_wide = vec![result.to_vec(), vec![0; 32]].concat();

    let val = F::from_bytes_wide(&bytes_wide.try_into().unwrap());

    /*
    let mut bytes_8 = vec![];
    for i in 0..(bytes.len() / 8) {
        let mut acc: u64 = 0;
        for j in 0..8 {
            acc += (bytes[8 * i + j] as u64) << j;
        }
        bytes_8.push(acc);
    }

    bytes_8.resize(25, 33);

    keccak::f1600(&mut bytes_8.as_slice().try_into().unwrap());

    let mut bytes = bytes_8[0..16]
        .iter()
        .flat_map(|x| x.to_le_bytes().to_vec())
        .collect::<Vec<u8>>();
    bytes.reverse();
    println!("bytes: {:?}", bytes);

    let val = F::from_bytes_wide(&bytes[0..64].try_into().unwrap());
     */

    val
}

fn sample_index(random_bytes: [u8; 64], size: usize) -> usize {
    let mut acc: u64 = 0;
    for b in random_bytes {
        acc = acc << 8 ^ (b as u64);
    }

    (acc % (size as u64)) as usize
}

pub fn sample_indices<F: FieldExt<Repr = [u8; 32]>>(
    num_indices: usize,
    max_index: usize,
    reduced_max_index: usize,
    transcript: &mut Transcript<F>,
) -> Vec<usize> {
    assert!(num_indices <= 2 * reduced_max_index, "not enough entropy!");
    assert!(num_indices <= reduced_max_index);

    let mut indices = vec![];

    let mut reduced_indices = vec![];

    let mut counter: u32 = 0;

    while indices.len() < num_indices {
        let mut random_bytes = [0u8; 64];

        transcript.append_bytes(&counter.to_le_bytes());
        transcript.challenge_bytes(&mut random_bytes);

        let index = sample_index(random_bytes, max_index);
        let reduced_index = index % reduced_max_index;

        counter += 1;
        if !reduced_indices.contains(&reduced_index) {
            reduced_indices.push(reduced_index);
            indices.push(index);
        };
    }

    indices
}
