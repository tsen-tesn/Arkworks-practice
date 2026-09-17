use ark_bls12_381::Bls12_381;
use ark_groth16::Groth16;
use ark_serialize::CanonicalSerialize;
use ark_std::test_rng;
use linear_transformation::MatrixLinearCircuit;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

fn exercise_path(file_name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(file_name)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. 建立只有電路結構、沒有實際 assignments 的 circuit
    let circuit = MatrixLinearCircuit::default();

    // 2. 建立 setup 使用的隨機數產生器
    let mut rng = test_rng();

    // 3. 執行 Groth16 trusted setup
    let pk = Groth16::<Bls12_381>::generate_random_parameters_with_reduction(circuit, &mut rng)?;

    // 4. 將 proving key 寫入 pk.bin
    let pk_file = File::create(exercise_path("pk.bin"))?;
    let mut pk_writer = BufWriter::new(pk_file);

    pk.serialize_compressed(&mut pk_writer)?;
    pk_writer.flush()?;

    // 5. 將 verifying key 寫入 vk.bin
    let vk_file = File::create(exercise_path("vk.bin"))?;
    let mut vk_writer = BufWriter::new(vk_file);

    pk.vk.serialize_compressed(&mut vk_writer)?;
    vk_writer.flush()?;

    println!("Groth16 setup 完成");
    println!("Proving key 已寫入 pk.bin");
    println!("Verifying key 已寫入 vk.bin");

    Ok(())
}
