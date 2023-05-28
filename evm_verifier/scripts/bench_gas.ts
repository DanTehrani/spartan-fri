import { ethers } from "hardhat";
import { exec } from "child_process";
import { promisify } from "util";
import fs from "fs";

const execAsync = promisify(exec);

async function main() {
  const SpartanSCVerifier = await ethers.getContractFactory("MultiLinearFRIVerifier");

  const numVars = 13;
  const numRound = 5;
  // Generate proof with syscall
  const lock = await SpartanSCVerifier.deploy(numVars, numRound);
  const result = await execAsync(`cd ../spartan_fri && cargo run --release --bin gen_proof ${numVars}`);

  const proof_bin = fs.readFileSync("../spartan_fri/proof.bin");
  console.log("proof_bin", proof_bin);
  console.log("proof_bin", proof_bin.byteLength);

  await lock.deployed();
}

// We recommend this pattern to be able to use async/await everywhere
// and properly handle errors.
main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
