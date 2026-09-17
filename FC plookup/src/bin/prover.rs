use ark_bls12_381::Fr;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::CanonicalSerialize;
use fc_plookup::{arithmetic, binding, columns::poly_from_rows, proof_io, MAX_DEGREE};
use merlin::Transcript;
use plookup::{
    kzg10,
    lookup::{
        lookup::LookUp,
        proof::{compress_column, encode_and_pad_witness, TABLE_SIZE},
        table::relu::{encode_signed, ReLUTable},
    },
    multiset::{proof::EqualityProof, MultiSet},
};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

fn exercise_path(file_name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(file_name)
}

/// Plookup 壓縮 (input, output) 兩欄用的隨機挑戰，prover 與 verifier 必須用同一個值。
const THETA: u64 = 7;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ---------- 0. 明文輸入：W, X, b（3x3 全連接層） ----------
    let w_i64: [[i64; 3]; 3] = [[1, 2, -1], [0, 1, 1], [-1, 0, 2]];
    let x_i64: [i64; 3] = [3, -2, 1];
    let b_i64: [i64; 3] = [1, -1, 2];

    let mut z_i64 = [0i64; 3];
    for i in 0..3 {
        let mut acc = b_i64[i];
        for j in 0..3 {
            acc += w_i64[i][j] * x_i64[j];
        }
        z_i64[i] = acc;
    }

    // z 必須落在公開 ReLU table 的範圍 [-32, 31] 內才能查表
    let mut lookup = LookUp::new(ReLUTable::new());
    for z in z_i64 {
        assert!(lookup.read(z), "z = {z} 不在 ReLU table 範圍內");
    }
    let y_i64: Vec<i64> = lookup.output_wires().to_vec();

    println!("z = {z_i64:?}");
    println!("y = {y_i64:?}");

    let theta = Fr::from(THETA);
    let (powers, _vk_kzg) = kzg10::trusted_setup(MAX_DEGREE);
    let domain = Radix2EvaluationDomain::<Fr>::new(TABLE_SIZE).unwrap();

    let mut prover_transcript = Transcript::new(b"fc-plookup-relu-layer");

    // ---------- 1. Arithmetic proof：直接在 C_z 上證明 z = W·X + b ----------
    let arithmetic_input = arithmetic::ProverInput {
        w: [
            [
                encode_signed(w_i64[0][0]),
                encode_signed(w_i64[0][1]),
                encode_signed(w_i64[0][2]),
            ],
            [
                encode_signed(w_i64[1][0]),
                encode_signed(w_i64[1][1]),
                encode_signed(w_i64[1][2]),
            ],
            [
                encode_signed(w_i64[2][0]),
                encode_signed(w_i64[2][1]),
                encode_signed(w_i64[2][2]),
            ],
        ],
        x: [
            encode_signed(x_i64[0]),
            encode_signed(x_i64[1]),
            encode_signed(x_i64[2]),
        ],
        b: [
            encode_signed(b_i64[0]),
            encode_signed(b_i64[1]),
            encode_signed(b_i64[2]),
        ],
        z: [
            encode_signed(z_i64[0]),
            encode_signed(z_i64[1]),
            encode_signed(z_i64[2]),
        ],
    };

    let (arithmetic_proof, z_poly) =
        arithmetic::prove(&arithmetic_input, &domain, &powers, &mut prover_transcript);

    // ---------- 2. y(X) 承諾：Plookup witness column 的另一半欄位 ----------
    let y_fr: Vec<Fr> = y_i64.iter().map(|v| encode_signed(*v)).collect();
    let y_poly = poly_from_rows(&domain, &y_fr);
    let y_commit = kzg10::commit(&powers, &y_poly);

    // ---------- 3. Plookup：witness column 直接由同一組 z_i64 / y_i64 構造，
    //              其 commitment 在數學上必然等於 C_z + theta * C_y（見 verifier 的同源性檢查） ----------
    let table = ReLUTable::new();
    let compressed_table = compress_column(
        &table.encode_input_column(),
        &table.encode_output_column(),
        theta,
    );
    let t = MultiSet::from_slice(&compressed_table);

    let (encoded_z, encoded_y) =
        encode_and_pad_witness(lookup.input_wires(), lookup.output_wires());
    let compressed_witness = compress_column(&encoded_z, &encoded_y, theta);
    let witness_ms = MultiSet::from_slice(&compressed_witness);

    let plookup_proof = EqualityProof::prove(witness_ms, t, &powers, &mut prover_transcript);

    // ---------- 4. y 的公開性：把 y(X) 在前 3 個 domain 點的值開給 verifier ----------
    let y_openings = binding::prove_openings(&y_poly, &domain, &y_fr, &powers);

    // ---------- 5. 落地存檔 ----------
    let mut arithmetic_writer =
        BufWriter::new(File::create(exercise_path("arithmetic_proof.data"))?);
    arithmetic_proof.serialize_compressed(&mut arithmetic_writer)?;
    arithmetic_writer.flush()?;

    let mut plookup_writer = BufWriter::new(File::create(exercise_path("plookup_proof.data"))?);
    proof_io::write_equality_proof(&plookup_proof, &mut plookup_writer)?;
    plookup_writer.flush()?;

    let mut y_commit_writer = BufWriter::new(File::create(exercise_path("y_commit.data"))?);
    y_commit.serialize_compressed(&mut y_commit_writer)?;
    y_commit_writer.flush()?;

    let mut y_openings_writer = BufWriter::new(File::create(exercise_path("y_openings.data"))?);
    y_openings.serialize_compressed(&mut y_openings_writer)?;
    y_openings_writer.flush()?;

    // z_poly 不需要落地：verifier 只需要 arithmetic_proof 裡的 C_z 就能做同源性檢查。
    let _ = z_poly;

    println!("Arithmetic proof（z = W·X + b）已寫入 arithmetic_proof.data");
    println!("Plookup ReLU proof 已寫入 plookup_proof.data");
    println!("y 的 commitment 與 openings 已寫入 y_commit.data / y_openings.data");

    Ok(())
}
