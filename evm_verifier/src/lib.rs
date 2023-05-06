use ark_std::{end_timer, start_timer};
use ethers::{
    contract::{abigen, ContractFactory},
    core::utils::Anvil,
    middleware::{contract::ContractError, SignerMiddleware},
    prelude::{
        k256::elliptic_curve::{consts::U25, Field},
        multicall_contract::Call3,
        ContractInstance,
    },
    providers::{Http, JsonRpcError, Provider},
    signers::{LocalWallet, Signer},
    solc::{artifacts::deserialize_bytes, Artifact, Project, ProjectPathsConfig},
    types::{TransactionReceipt, H160, U256},
};
use eyre::Result;
use pasta_curves::arithmetic::FieldExt;
use pasta_curves::group::ff::PrimeField;
use serde_json::de::Deserializer;
use spartan_fri::{SpartanFRIProof, SpartanFRIProver, Transcript, R1CS};
use std::{os::unix::prelude::FileTypeExt, path::PathBuf, sync::Arc, time::Duration};

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

async fn main() -> Result<()> {
    type F = pasta_curves::Fp;

    // 2. instantiate our wallet & anvil
    let anvil = Anvil::new().spawn();
    let wallet: LocalWallet = anvil.keys()[0].clone().into();

    //    3. connect to the network
    let provider =
        Provider::<Http>::try_from(anvil.endpoint())?.interval(Duration::from_millis(10u64));

    // 4. instantiate the client with the wallet
    let client = SignerMiddleware::new(provider, wallet.with_chain_id(anvil.chain_id()));
    let client = Arc::new(client);

    Ok(())
}

#[tokio::test]
async fn test_main() {
    let result = main().await;

    println!("result: {:?}", result);
}
