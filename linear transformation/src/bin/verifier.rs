use ark_bls12_381::{Bls12_381, Fr};
use ark_groth16::{
    prepare_verifying_key,
    Groth16,
    Proof,
    VerifyingKey,
};
use ark_serialize::CanonicalDeserialize;
use ark_snark::SNARK;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

fn exercise_path(file_name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(file_name)
}

fn f(value: u64) -> Fr {
    Fr::from(value)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. 開啟 vk.bin
    let vk_file = File::open(exercise_path("vk.bin"))?;
    let mut vk_reader = BufReader::new(vk_file);

    // 2. 反序列化 verifying key
    let vk =
        VerifyingKey::<Bls12_381>::deserialize_compressed(
            &mut vk_reader,
        )?;

    // 3. 預先計算 verifying key
    let prepared_vk = prepare_verifying_key(&vk);

    // 4. 開啟 compressed proof
    let proof_file = File::open(exercise_path("proof_compressed.data"))?;
    let mut proof_reader = BufReader::new(proof_file);

    // 5. 反序列化 proof
    let proof =
        Proof::<Bls12_381>::deserialize_compressed(
            &mut proof_reader,
        )?;

    // 6. 宣告 public inputs
    //
    // 配置順序必須和 circuit 中的 y_vars 相同：
    // y[0] = 15
    // y[1] = 24
    let public_inputs = [
        f(15),
        f(24),
    ];

    // 7. 驗證正確的 public inputs
    let is_valid =
        Groth16::<Bls12_381>::verify_with_processed_vk(
            &prepared_vk,
            &public_inputs,
            &proof,
        )?;

    println!("正確 public inputs 的驗證結果：{is_valid}");

    // 8. 額外測試錯誤的 public inputs
    let wrong_public_inputs = [
        f(15),
        f(25),
    ];

    let wrong_is_valid =
        Groth16::<Bls12_381>::verify_with_processed_vk(
            &prepared_vk,
            &wrong_public_inputs,
            &proof,
        )?;

    println!(
        "錯誤 public inputs 的驗證結果：{wrong_is_valid}"
    );

    Ok(())
}
