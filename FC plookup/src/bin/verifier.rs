use ark_bls12_381::{Bls12_381, Fr};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_poly_commit::kzg10::Commitment;
use ark_serialize::CanonicalDeserialize;
use fc_plookup::{
    arithmetic::{self, ArithmeticProof},
    binding::{self, OpeningsProof},
    proof_io, MAX_DEGREE,
};
use merlin::Transcript;
use plookup::{
    kzg10,
    lookup::{
        proof::{compress_column, TABLE_SIZE},
        table::relu::{encode_signed, ReLUTable},
    },
    multiset::MultiSet,
};
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

fn exercise_path(file_name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(file_name)
}

/// 必須跟 prover 使用同一個壓縮挑戰。
const THETA: u64 = 7;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 手動宣告的公開資訊：只有最終輸出 y。z 不再需要手動宣告或事先知道，
    // 它只活在 arithmetic proof 與 Plookup proof 共用的 commitment C_z 裡。
    let y_public: [i64; 3] = [0, 0, 1];
    let y_fr: Vec<Fr> = y_public.iter().map(|v| encode_signed(*v)).collect();

    let theta = Fr::from(THETA);
    let (powers, vk) = kzg10::trusted_setup(MAX_DEGREE);
    let domain = Radix2EvaluationDomain::<Fr>::new(TABLE_SIZE).unwrap();

    let arithmetic_proof = ArithmeticProof::deserialize_compressed(BufReader::new(File::open(
        exercise_path("arithmetic_proof.data"),
    )?))?;
    let plookup_proof =
        proof_io::read_equality_proof(BufReader::new(File::open(exercise_path(
            "plookup_proof.data",
        ))?))?;
    let y_commit = Commitment::<Bls12_381>::deserialize_compressed(BufReader::new(File::open(
        exercise_path("y_commit.data"),
    )?))?;
    let y_openings = OpeningsProof::deserialize_compressed(BufReader::new(File::open(
        exercise_path("y_openings.data"),
    )?))?;

    // 必須跟 prover 使用同一個 transcript label / 附加順序，才能重現一樣的挑戰值
    let mut transcript = Transcript::new(b"fc-plookup-relu-layer");

    // ---------- 1. 驗證 Arithmetic proof：z(X) 滿足 z = W·X + b ----------
    let arithmetic_ok = arithmetic::verify(&arithmetic_proof, &domain, &vk, &mut transcript);
    println!("Arithmetic proof（z = W·X + b，對 C_z 的 quotient 檢查）驗證結果：{arithmetic_ok}");

    // ---------- 2. 驗證 Plookup：(z_i, y_i) 都在公開 ReLU table 內 ----------
    let table = ReLUTable::new();
    let compressed_table = compress_column(
        &table.encode_input_column(),
        &table.encode_output_column(),
        theta,
    );
    let t = MultiSet::from_slice(&compressed_table);
    let t_poly = t.to_polynomial(&domain);
    let t_commit = kzg10::commit(&powers, &t_poly);

    let plookup_ok = plookup_proof.verify(TABLE_SIZE, &vk, t_commit, &mut transcript);
    println!("Plookup (z_i, y_i) ⊆ ReLU table 驗證結果：{plookup_ok}");

    // ---------- 3. 同源性檢查：Plookup 用的 witness commitment 真的是同一個 C_z ----------
    // f(X) = z(X) + theta * y(X)（IFFT 是線性運算），所以 KZG commitment 也滿足
    // C_f = C_z + theta * C_y。這是一個精確的群元素等式，不是機率性的取樣檢查。
    let expected_c_f =
        kzg10::aggregate_commitments(vec![&arithmetic_proof.commitments.z, &y_commit], theta);
    let binding_ok = expected_c_f == plookup_proof.commitments.f;
    println!("C_f == C_z + theta * C_y（z 在兩段證明中是同一個 commitment）驗證結果：{binding_ok}");

    // ---------- 4. 驗證 C_y 真的對應到宣稱的公開輸出 y ----------
    let y_opening_ok = binding::verify_openings(&y_commit, &y_openings, &domain, &y_fr, &vk);
    println!("y 的 commitment 對應到宣稱的公開輸出驗證結果：{y_opening_ok}");

    let overall_ok = arithmetic_ok && plookup_ok && binding_ok && y_opening_ok;
    println!("整體證明 y = ReLU(W·X + b) 驗證結果：{overall_ok}");

    Ok(())
}
