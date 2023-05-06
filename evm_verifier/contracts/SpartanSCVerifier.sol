// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.19;

contract SpartanSumCheckVerifier {
    // Just store the claims so we can simulate state updates
    uint256[] public sumClaims;
    uint256[] public finalPolyEvalClaims;

    // Modulus of the Pallas curve
    uint256 constant MODULUS =
        0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001;
    uint256 POLY_NUM_VARS;

    constructor(uint256 _POLY_NUM_VARS) {
        POLY_NUM_VARS = _POLY_NUM_VARS;
    }

    function evalRoundPoly(
        uint256[] memory poly,
        uint256 x
    ) internal pure returns (uint256) {
        if (poly.length == 4) {
            uint256 x_sq = mulmod(x, x, MODULUS);
            uint256 x_cub = mulmod(x_sq, x, MODULUS);

            uint256 evalHigh = addmod(
                mulmod(poly[0], x_cub, MODULUS),
                mulmod(poly[1], x_sq, MODULUS),
                MODULUS
            );

            uint256 evalLow = addmod(
                mulmod(poly[2], x, MODULUS),
                poly[3],
                MODULUS
            );
            uint256 eval = addmod(evalHigh, evalLow, MODULUS);
            return eval;
        } else if (poly.length == 3) {
            uint256 x_sq = mulmod(x, x, MODULUS);
            uint256 evalHigh = addmod(
                mulmod(poly[0], x_sq, MODULUS),
                mulmod(poly[1], x, MODULUS),
                MODULUS
            );

            return addmod(evalHigh, poly[2], MODULUS);
        } else {
            revert("Invalid poly length");
        }
    }

    event log(string);

    function _verify(
        uint256[] memory roundPolys,
        uint256[] memory challenge,
        uint256 polyDegree
    ) internal returns (uint256) {
        uint256[] memory firstRoundPoly = new uint256[](polyDegree);

        for (uint256 i = 0; i < polyDegree; i++) {
            firstRoundPoly[i] = roundPolys[i];
        }

        uint256 firstPolyEval0 = evalRoundPoly(firstRoundPoly, 0);
        uint256 firstPolyEval1 = evalRoundPoly(firstRoundPoly, 1);
        uint256 sumClaim = addmod(firstPolyEval0, firstPolyEval1, MODULUS);

        for (uint256 i = 0; i < POLY_NUM_VARS - 1; i++) {
            uint256 roundPolyIndex = polyDegree * i;
            uint256 nextRondPolyIndex = polyDegree * (i + 1);
            uint256 r = challenge[i];

            uint256[] memory roundPoly = new uint256[](polyDegree);
            uint256[] memory nextRoundPoly = new uint256[](polyDegree);

            for (uint256 j = 0; j < polyDegree; j++) {
                roundPoly[j] = roundPolys[roundPolyIndex + j];
                nextRoundPoly[j] = roundPolys[nextRondPolyIndex + j];
            }

            uint256 eval0 = evalRoundPoly(nextRoundPoly, 0);
            uint256 eval1 = evalRoundPoly(nextRoundPoly, 1);
            uint256 lhs = addmod(eval0, eval1, MODULUS);
            uint256 rhs = evalRoundPoly(roundPoly, r);

            // require(lhs == rhs, "Sum check failed");
        }

        return sumClaim;
    }

    function challenge_vec(uint256 n) internal view returns (uint256[] memory) {
        uint256[] memory challenge = new uint256[](n);
        for (uint256 i = 0; i < n; i++) {
            challenge[i] = uint256(
                keccak256(abi.encodePacked(block.timestamp, i))
            );
        }

        return challenge;
    }

    function verify(
        uint256[] calldata _phase1RoundPolys,
        uint256[] calldata _phase2RoundPolys
    ) external {
        uint256[] memory phase1RoundPolys = _phase1RoundPolys;
        uint256[] memory phase2RoundPolys = _phase2RoundPolys;

        uint256[] memory rx = challenge_vec(POLY_NUM_VARS);

        uint256 sumClaimPhase1 = _verify(phase1RoundPolys, rx, 4);

        // uint256[] memory r = challenge_vec(3);

        uint256[] memory ry = challenge_vec(POLY_NUM_VARS);

        uint256 sumClaimPhase2 = _verify(phase2RoundPolys, ry, 3);

        sumClaims.push(sumClaimPhase1);
        sumClaims.push(sumClaimPhase2);
    }
}
