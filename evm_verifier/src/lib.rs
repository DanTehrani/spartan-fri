use ethers::{contract::abigen, types::U256};
use pasta_curves::arithmetic::FieldExt;
use spartan_fri::SpartanFRIProof;

// Generate the type-safe contract bindings by providing the ABI
// definition in human readable format
abigen!(
    SpartanSumCheckVerifier,
    r#"[
        function verify(uint256[] calldata _phase1RoundPolys, uint256[] calldata _phase2RoundPolys) external
    ]"#,
    event_derives(serde::Deserialize, serde::Serialize)
);

pub trait ToCallData<F: FieldExt<Repr = [u8; 32]>> {
    type Calldata;
    fn to_calldata(&self) -> Self::Calldata;
}

pub struct SpartanFRIProofCalldata {
    pub phase1_round_polys: Vec<U256>,
    pub phase2_round_polys: Vec<U256>,
}

impl<F: FieldExt<Repr = [u8; 32]>> ToCallData<F> for SpartanFRIProof<F> {
    type Calldata = SpartanFRIProofCalldata;
    fn to_calldata(&self) -> Self::Calldata {
        let phase1_round_polys = self
            .sc_proof_1
            .round_polys
            .iter()
            .map(|round_poly| {
                round_poly
                    .coeffs
                    .iter()
                    .map(|coeff| U256::from_big_endian(&coeff.to_repr()))
                    .collect::<Vec<U256>>()
            })
            .flatten()
            .collect::<Vec<U256>>();

        let phase2_round_polys = self
            .sc_proof_2
            .round_polys
            .iter()
            .map(|round_poly| {
                round_poly
                    .coeffs
                    .iter()
                    .map(|coeff| U256::from_big_endian(&coeff.to_repr()))
                    .collect::<Vec<U256>>()
            })
            .flatten()
            .collect::<Vec<U256>>();

        SpartanFRIProofCalldata {
            phase1_round_polys,
            phase2_round_polys,
        }
    }
}
