use ark_bls12_381::Fr;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use merlin::Transcript;
use plookup::{
    kzg10,
    lookup::{
        lookup::LookUp,
        proof::{compress_column, encode_and_pad_witness, TABLE_SIZE},
        table::relu::ReLUTable,
    },
    multiset::{proof::EqualityProof, MultiSet},
};

/// 壓縮 table 並承諾 t(X)，是 prover 和 verifier 共用的前置步驟。
fn setup_table(
    powers: &ark_poly_commit::kzg10::Powers<ark_bls12_381::Bls12_381>,
    theta: Fr,
) -> (MultiSet, ark_poly_commit::kzg10::Commitment<ark_bls12_381::Bls12_381>) {
    let table = ReLUTable::new();
    let compressed = compress_column(
        &table.encode_input_column(),
        &table.encode_output_column(),
        theta,
    );
    let t = MultiSet::from_slice(&compressed);
    let domain: Radix2EvaluationDomain<Fr> =
        EvaluationDomain::new(t.len()).unwrap();
    let t_poly = t.to_polynomial(&domain);
    let t_commit = kzg10::commit(powers, &t_poly);
    (t, t_commit)
}

/// 合法 witness：整數 ReLU 案例，proof 應該驗證成功。
#[test]
fn valid_relu_witness_proof_verifies() {
    // Trusted Setup：max_degree 要大於 quotient poly 的 degree（約 3n = 189）
    let (powers, vk) = kzg10::trusted_setup(256);

    let theta = Fr::from(7u64);
    let (t, _) = setup_table(&powers, theta);

    // 建立 witness
    let mut lookup = LookUp::new(ReLUTable::new());
    for input in [2, 4, -1, 0, 15, 6, -24, -2, 15] {
        assert!(lookup.read(input), "read 失敗：{input} 不在 table 中");
    }
    let (enc_in, enc_out) =
        encode_and_pad_witness(lookup.input_wires(), lookup.output_wires());
    let compressed_witness = compress_column(&enc_in, &enc_out, theta);
    let f = MultiSet::from_slice(&compressed_witness);

    // Prove
    let mut prover_transcript = Transcript::new(b"plookup-relu");
    let proof = EqualityProof::prove(f, t, &powers, &mut prover_transcript);

    // Verify（transcript 重置，走同一條路得到相同挑戰）
    let (_, t_commit_v) = setup_table(&powers, theta);
    let mut verifier_transcript = Transcript::new(b"plookup-relu");
    let ok = proof.verify(TABLE_SIZE, &vk, t_commit_v, &mut verifier_transcript);
    assert!(ok, "合法 witness 應該驗證成功");
}