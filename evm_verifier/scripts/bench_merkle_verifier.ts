import { ethers } from "hardhat";
import { exec } from "child_process";
import { promisify } from "util";
import fs from "fs";

const execAsync = promisify(exec);

type PrefixedHex = `0x${string}`;

const bigIntToPrefixHex = (num: BigInt): PrefixedHex => {
  const hex = num.toString(16).split("").reverse().join("");
  return `0x${hex}`;
};

const deserialize = (bytes: Buffer): bigint[] => {
  const deserialized: bigint[] = [];
  for (let i = 0; i < bytes.byteLength / 32; i++) {
    let hash = BigInt(0);
    for (let j = 0; j < 32; j++) {
      hash += BigInt(bytes[i * 32 + j]) << BigInt(8 * j);
    }
    deserialized.push(hash);
  }
  return deserialized;
};

async function main() {
  const MerkleVerifier = await ethers.getContractFactory("MerkleVerifier");

  const depth = 9;
  const batchSize = 2;
  // Generate proof with syscall
  const merkleVerifier = await (
    await MerkleVerifier.deploy(depth, batchSize)
  ).deployed();

  const result = await execAsync(
    `cd ../spartan_fri && cargo run --release --bin gen_batched_merkle_proof ${depth} ${batchSize}`
  );

  console.log(result.stdout);

  const hashes_bin = fs.readFileSync("../spartan_fri/hashes.bin");
  const indices_bin = fs.readFileSync("../spartan_fri/indices.bin");
  const leafIndicesBin = fs.readFileSync("../spartan_fri/leaf_indices.bin");
  const leaves_bin = fs.readFileSync("../spartan_fri/leaves.bin");
  const root_bin = fs.readFileSync("../spartan_fri/root.bin");

  const proofSize =
    hashes_bin.byteLength +
    indices_bin.byteLength +
    leafIndicesBin.byteLength +
    leaves_bin.byteLength;

  console.log("Num hashes", hashes_bin.byteLength / 32);
  console.log("Num indices", indices_bin.byteLength / 32);
  console.log("depth", depth);
  console.log("batch size", batchSize);
  console.log("proof size", proofSize, "bytes");

  const hashes = deserialize(hashes_bin);
  const indices = deserialize(indices_bin);
  const leafIndices = deserialize(leafIndicesBin);
  const leaves = deserialize(leaves_bin);
  const root = deserialize(root_bin)[0];

  const tx = await merkleVerifier.verifyBatch(
    hashes.map(hash => Buffer.from(hash.toString(16).padStart(64, "0"), "hex")),
    indices,
    leaves.map(hash => Buffer.from(hash.toString(16).padStart(64, "0"), "hex")),
    leafIndices,
    Buffer.from(root.toString(16).padStart(64, "0"), "hex")
  );
  const receipt = await tx.wait();
  console.log("gas used", receipt.gasUsed.toString());
}

main().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
