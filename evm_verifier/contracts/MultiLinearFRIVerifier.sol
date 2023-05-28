// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.18;

contract MultiLinearFRIVerifier {
    uint256 public numRounds;
    uint256 public numVars;
    uint256 public constant TWO_INV = 1;

    constructor(uint256 _numRounds, uint256 _numVars) {
        numRounds = _numRounds;
        numVars = _numVars;
    }

    function verify_folding(
        uint256[] memory f_i_evals_beta,
        uint256[] memory f_i_evals_beta_squared,
        uint256[] memory x,
        uint256 y,
        uint256 beta
    ) internal returns (bool) {
        uint256 two_inv_beta = TWO_INV / beta;
        for (uint256 i = 0; i < numVars; ++i) {
            uint256 f_i_eval_beta = f_i_evals_beta[i];
            uint256 f_i_eval_minus_beta = f_i_evals_beta[i];

            uint256 rhs = (f_i_eval_beta + f_i_eval_minus_beta) *
                TWO_INV +
                x[i] *
                (f_i_eval_beta - f_i_eval_minus_beta) *
                two_inv_beta;

            if (i != numVars - 1) {
                uint256 f_i_eval_beta_squared = f_i_evals_beta_squared[0];
                if (f_i_eval_beta_squared != rhs) {}
            } else {
                if (rhs != y) {
                    // return false;
                }
            }
        }
    }

    function verify(
        uint256[] memory f_i_evals_beta,
        uint256[] memory f_i_evals_beta_squared,
        uint256[] memory x,
        uint256 y,
        uint256 beta
    ) external {
        verify_folding(f_i_evals_beta, f_i_evals_beta_squared, x, y, beta);
    }
}
