use crate::gate::layout::{build_y_column, build_z_copy_column, compute_wires, MATRIX_DIM};
use crate::gate::proof::{self as gate_proof, GateProof};
use ark_bls12_381::{Bls12_381, Fr};
use ark_poly::{univariate::DensePolynomial as Polynomial, DenseUVPolynomial, EvaluationDomain};
use ark_poly_commit::kzg10::{Commitment, Powers, VerifierKey};
use merlin::Transcript;
use plookup::kzg10;
use plookup::lookup::proof::{compress_column, TABLE_SIZE, WITNESS_SIZE};
use plookup::lookup::table::relu::ReLUTable;
use plookup::multiset::proof::EqualityProof;
use plookup::multiset::MultiSet;
use plookup::transcript::TranscriptProtocol;

/// 整個 FC_plookup 的最終證明：「我知道私密的 w,x,b,y，使得 z=Wx+b 算術/wiring 正確，
/// 而且每個 (z_i, y_i) 都落在公開的 ReLU table 裡」，全程 w,x,b,z,y 都不會被揭露。
///
/// M1+M2+M3：gate_proof 證明 z=Wx+b 且 wiring 正確，並額外 commit 了一條跟 z 同態綁定的
/// D 欄位（gate_proof.commitments.d）。
/// M4：y_commit 是 y 的獨立 commitment；binding 靠 verifier 自己用同態性檢查
/// `commit(D) + theta * commit(Y) == equality_proof 內部的 f_commit`，
/// 而不是要求 prover「自稱」f 是怎麼來的。
/// M5：因為 D 欄位本身就跟 lookup witness 活在同一個 64-domain、同一組 row 位置上，
/// 不需要額外的跨 domain 開值——這是跟 design.md 原始方案最大的不同之處。
/// M6：gate_proof 跟 equality_proof 共用同一條 Fiat-Shamir transcript。
pub struct FCProof {
    pub gate_proof: GateProof,
    pub y_commit: Commitment<Bls12_381>,
    pub equality_proof: EqualityProof,
}

fn compressed_table(theta: Fr) -> MultiSet {
    let table = ReLUTable::new();
    let compressed = compress_column(
        &table.encode_input_column(),
        &table.encode_output_column(),
        theta,
    );
    MultiSet::from_slice(&compressed)
}

pub fn prove(
    w: &[[Fr; MATRIX_DIM]; MATRIX_DIM],
    x: &[Fr; MATRIX_DIM],
    b: &[Fr; MATRIX_DIM],
    y: &[Fr; MATRIX_DIM],
    powers: &Powers<Bls12_381>,
) -> FCProof {
    let mut transcript = Transcript::new(b"fc-plookup");

    let (gate_proof, _z_values) = gate_proof::prove(w, x, b, powers, &mut transcript);

    let domain = gate_proof::domain();
    let y_evals = build_y_column(y);
    let y_poly = Polynomial::from_coefficients_vec(domain.ifft(&y_evals));
    let y_commit = kzg10::commit(powers, &y_poly);
    transcript.append_commitment(b"y_commit", &y_commit);

    let theta = transcript.challenge_scalar(b"theta");

    // D、Y 欄位長度是 NUM_GATES(=64)，跟 lookup 需要的 WITNESS_SIZE(=63) 只差最後一格
    // padding（已驗證恆為 0，見 gate::layout 的 z_copy_column_matches_and_pads_with_zero 測試），
    // 所以直接取前 WITNESS_SIZE 個元素跟 encode_and_pad_witness 產生的形狀完全一致。
    let layout = crate::gate::layout::build_gate_layout();
    let (_, _, c_evals) = compute_wires(&layout, w, x, b);
    let z_copy_full = build_z_copy_column(&c_evals);
    let f_inputs = &z_copy_full[..WITNESS_SIZE];
    let f_outputs = &y_evals[..WITNESS_SIZE];
    let f_compressed = compress_column(f_inputs, f_outputs, theta);
    let f = MultiSet::from_slice(&f_compressed);

    let t = compressed_table(theta);

    let equality_proof = EqualityProof::prove(f, t, powers, &mut transcript);

    FCProof {
        gate_proof,
        y_commit,
        equality_proof,
    }
}

/// `powers` 在這個練習版本裡也被 verifier 用來重新 commit 公開的 ReLU table
/// （table 的壓縮方式依賴 theta，theta 是證明產生當下才決定的挑戰，沒辦法在
/// trusted setup 階段就固定 commit 起來）。真正部署時通常會希望 verifier
/// 只需要 VerifierKey，這裡為了跟現有 plookup 測試（見 plookup/tests/lookup.rs
/// 的 setup_table）保持一致的做法而簡化，沒有解決這個問題。
pub fn verify(
    proof: &FCProof,
    powers: &Powers<Bls12_381>,
    vk: &VerifierKey<Bls12_381>,
) -> bool {
    let mut transcript = Transcript::new(b"fc-plookup");

    if !proof.gate_proof.verify(vk, &mut transcript) {
        return false;
    }

    transcript.append_commitment(b"y_commit", &proof.y_commit);
    let theta = transcript.challenge_scalar(b"theta");

    let t = compressed_table(theta);
    let domain = gate_proof::domain();
    let t_commit = kzg10::commit(powers, &t.to_polynomial(&domain));

    if !proof
        .equality_proof
        .verify(TABLE_SIZE, vk, t_commit, &mut transcript)
    {
        return false;
    }

    // M4：verifier 自己用同態性把 D、Y 的 commitment 組合起來，
    // 檢查跟 lookup 論證內部宣稱的 f_commit 是不是同一個群元素——
    // 不需要開值、不需要相信 prover，z 全程沒有以明文形式出現過。
    let expected_f_commit = kzg10::aggregate_commitments(
        vec![&proof.gate_proof.commitments.d, &proof.y_commit],
        theta,
    );

    expected_f_commit == proof.equality_proof.commitments.f
}

#[cfg(test)]
mod tests {
    use super::*;

    fn f(v: u64) -> Fr {
        Fr::from(v)
    }

    fn sample_witness() -> (
        [[Fr; MATRIX_DIM]; MATRIX_DIM],
        [Fr; MATRIX_DIM],
        [Fr; MATRIX_DIM],
        [Fr; MATRIX_DIM],
    ) {
        // z0 = 1+2+3+1 = 7,  z1 = 4+5+6+1 = 16,  z2 = 7+8+9+1 = 25
        // ReLU(z) = z（全部是正數），所以合法的 y 就是 z 本身
        let w = [
            [f(1), f(2), f(3)],
            [f(4), f(5), f(6)],
            [f(7), f(8), f(9)],
        ];
        let x = [f(1), f(1), f(1)];
        let b = [f(1), f(1), f(1)];
        let y = [f(7), f(16), f(25)];
        (w, x, b, y)
    }

    #[test]
    fn valid_fc_proof_verifies() {
        let (powers, vk) = kzg10::trusted_setup(1024);
        let (w, x, b, y) = sample_witness();

        let proof = prove(&w, &x, &b, &y, &powers);
        assert!(verify(&proof, &powers, &vk));
    }

    #[test]
    fn swapped_z_commitment_is_rejected_by_binding_check() {
        // 這是整個綁定機制要防的攻擊，也是最初 design.md 提出的情境：
        // 換一個完全合法（gate identity + wiring 都對）但對應到不同 w,x,b 的
        // gate_proof 塞進同一個 FCProof，binding check（commitment 相等）必須抓到，
        // 光靠「gate_proof 自己 verify 會過」跟「lookup 自己 verify 會過」都不夠。
        let (powers, vk) = kzg10::trusted_setup(1024);
        let (w, x, b, y) = sample_witness();

        let mut proof = prove(&w, &x, &b, &y, &powers);

        // 換成另一組「z,y 也都合法」的 gate_proof（對應到不同的 w,x,b），
        // 藉此拿到一個 D commitment 跟原本的 y_commit 對不上的組合。
        let other_w = [
            [f(2), f(2), f(2)],
            [f(2), f(2), f(2)],
            [f(2), f(2), f(2)],
        ];
        let other_x = [f(1), f(1), f(1)];
        let other_b = [f(0), f(0), f(0)];
        let mut swap_transcript = Transcript::new(b"fc-plookup");
        let (other_gate_proof, _) =
            gate_proof::prove(&other_w, &other_x, &other_b, &powers, &mut swap_transcript);
        proof.gate_proof = other_gate_proof;

        assert!(!verify(&proof, &powers, &vk));
    }
}