// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.18;

import "hardhat/console.sol";

contract MerkleVerifier {
  uint256 public depth;
  uint256 public batchSize;

  constructor(uint256 _depth, uint256 _batchSize) {
    depth = _depth;
    batchSize = _batchSize;
  }

  function hashTwo(bytes32 a, bytes32 b) internal pure returns (bytes32) {
    bytes32 hashed = keccak256(abi.encodePacked(a, b));
    return hashed;
  }

  function verify(
    bytes32[] memory hashes,
    uint16[] memory indices,
    uint16 leafIndex,
    bytes32 leaf,
    bytes32 root
  ) internal view returns (bool) {
    for (uint256 i = 0; i < depth; ++i) {
      uint16 index = indices[i];
      if (leafIndex & 1 == 0) {
        leaf = hashTwo(leaf, hashes[index]);
      } else {
        leaf = hashTwo(hashes[index], leaf);
      }

      leafIndex >>= 1;
    }

    return leaf == root;
  }

  function verifyBatch(
    bytes32[] calldata hashes,
    uint16[] calldata indices, // index per byte
    bytes32[] calldata leaves,
    uint16[] calldata leafIndices, // index per byte
    bytes32 root
  ) external {
    for (uint256 i = 0; i < batchSize; ++i) {
      uint16[] memory indices_i = indices[(i * depth):(i * depth + depth)];
      bool result = verify(hashes, indices_i, leafIndices[i], leaves[i], root);
      assert(result);
    }
  }
}
