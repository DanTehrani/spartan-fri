use ethers::{
    contract::ContractFactory,
    core::utils::Anvil,
    middleware::SignerMiddleware,
    providers::{Http, Provider},
    signers::{LocalWallet, Signer},
    solc::{Artifact, Project, ProjectPathsConfig},
    types::{TransactionReceipt, H160, U256},
};
use evm_verifier::{SpartanSumCheckVerifier, ToCallData};
use eyre::Result;
use pasta_curves::arithmetic::FieldExt;
use spartan_fri::{SpartanFRIProof, SpartanFRIProver, SpartanPP, R1CS};
use std::{path::PathBuf, sync::Arc, time::Duration};
use tokio::runtime::Runtime;

async fn deploy_sc_verifier(
    client: &Arc<SignerMiddleware<Provider<Http>, LocalWallet>>,
    poly_num_vars: usize,
) -> Result<H160> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("contracts");
    let paths = ProjectPathsConfig::builder()
        .root(&root)
        .sources(&root)
        .build()
        .unwrap();

    // get the solc project instance using the paths above
    let project = Project::builder()
        .paths(paths)
        .ephemeral()
        .no_artifacts()
        .build()
        .unwrap();

    // compile the project and get the artifacts
    let output = project.compile().unwrap();
    let contract = output
        .find_first("SpartanSumCheckVerifier")
        .expect("could not find contract")
        .clone();

    let (abi, bytecode, _) = contract.into_parts();

    // create a factory which will be used to deploy instances of the contract
    let factory = ContractFactory::new(abi.unwrap(), bytecode.unwrap(), client.clone());

    // deploy it with the constructor arguments
    let contract = factory.deploy(U256::from(poly_num_vars))?.send().await?;

    //  get the contract's address
    let addr = contract.address();

    Ok(addr)
}

fn gen_proof<F: FieldExt<Repr = [u8; 32]>>(num_vars: usize) -> SpartanFRIProof<F> {
    let num_cons = num_vars;
    let num_vars = num_cons;
    let num_input = 5;

    let (r1cs, witness) = R1CS::produce_synthetic_r1cs(num_cons, num_vars, num_input);

    let pp = SpartanPP::new(r1cs, b"test_evm_verify");
    let prover: SpartanFRIProver<F> = SpartanFRIProver::<F>::new(pp);
    prover.prove(&witness)
}

async fn bench_gas(client: Arc<SignerMiddleware<Provider<Http>, LocalWallet>>) -> Result<()> {
    type F = pasta_curves::Fp;

    for num_var_exp in 10..14 {
        let num_vars = 2usize.pow(num_var_exp);
        let contract_addr = deploy_sc_verifier(&client, num_var_exp as usize).await?;
        let proof = gen_proof::<F>(num_vars);

        let proof_calldata = proof.to_calldata();

        let verifier = SpartanSumCheckVerifier::new(contract_addr, client.clone());
        let result: Option<TransactionReceipt> = verifier
            .verify(
                proof_calldata.phase1_round_polys,
                proof_calldata.phase2_round_polys,
            )
            .send()
            .await?
            .await?;

        println!(
            "Gas used for num_vars {}: {:?}",
            num_vars,
            result.unwrap().gas_used
        );
    }

    Ok(())
}

fn main() {
    let rt = Runtime::new().unwrap();

    // 2. instantiate our wallet & anvil
    let anvil = Anvil::new().spawn();
    let wallet: LocalWallet = anvil.keys()[0].clone().into();

    //    3. connect to the network
    let provider = Provider::<Http>::try_from(anvil.endpoint())
        .unwrap()
        .interval(Duration::from_millis(10u64));

    // 4. instantiate the client with the wallet
    let client = SignerMiddleware::new(provider, wallet.with_chain_id(anvil.chain_id()));
    let client = Arc::new(client);

    rt.block_on(async move {
        bench_gas(client).await.unwrap();
    });
}
