use ark_bls12_381::{Bls12_381, Fr};
use ark_groth16::{Groth16, ProvingKey};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_snark::SNARK;
use ark_std::rand::{rngs::StdRng, SeedableRng};
use linear_transformation::MatrixLinearCircuit;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

fn exercise_path(file_name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(file_name)
}

fn f(value: u64) -> Fr {
    Fr::from(value)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. 開啟 pk.bin
    let pk_file = File::open(exercise_path("pk.bin"))?;
    let mut pk_reader = BufReader::new(pk_file);

    // 2. 反序列化 proving key
    let pk = ProvingKey::<Bls12_381>::deserialize_compressed(&mut pk_reader)?;

    // 3. 建立帶有實際 assignments 的 circuit
    let circuit = MatrixLinearCircuit {
        w: [[Some(f(2)), Some(f(1))], [Some(f(3)), Some(f(4))]],
        x: [Some(f(5)), Some(f(2))],
        b: [Some(f(3)), Some(f(1))],
        y: [Some(f(15)), Some(f(24))],
    };

    // 4. 建立 prover 使用的隨機數產生器
    let mut rng = StdRng::seed_from_u64(42);

    // 5. 產生 Groth16 proof
    let proof = Groth16::<Bls12_381>::prove(&pk, circuit, &mut rng)?;

    // 6. 儲存 compressed proof
    let compressed_file = File::create(exercise_path("proof_compressed.data"))?;

    let mut compressed_writer = BufWriter::new(compressed_file);

    proof.serialize_compressed(&mut compressed_writer)?;
    compressed_writer.flush()?;

    // 7. 儲存 uncompressed proof
    let uncompressed_file = File::create(exercise_path("proof_uncompressed.data"))?;

    let mut uncompressed_writer = BufWriter::new(uncompressed_file);

    proof.serialize_uncompressed(&mut uncompressed_writer)?;
    uncompressed_writer.flush()?;

    println!("Groth16 proof 產生完成");
    println!("Compressed proof 已寫入 proof_compressed.data");
    println!("Uncompressed proof 已寫入 proof_uncompressed.data");

    Ok(())
}
