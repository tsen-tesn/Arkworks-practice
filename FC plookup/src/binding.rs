use ark_bls12_381::{Bls12_381, Fr};
use ark_poly::{
    polynomial::univariate::DensePolynomial as UniPoly, EvaluationDomain, Polynomial as _,
    Radix2EvaluationDomain,
};
use ark_poly_commit::kzg10::{Commitment, Powers, VerifierKey};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use plookup::kzg10;

/// z 是否為同一組數字現在靠 arithmetic proof 與 Plookup commitment 的同源性直接檢查
/// （見 verifier.rs 裡的 `C_f == C_z + theta * C_y`），不需要另外開洞。
/// 這裡只需要把「最終公開輸出」y(X) 在前 n 個 domain 點的值開給 verifier，
/// 確認 y(X) 真的對應到一組公開已知的數值。
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct OpeningsProof {
    pub witness_commits: Vec<Commitment<Bls12_381>>,
}

pub fn prove_openings(
    poly: &UniPoly<Fr>,
    domain: &Radix2EvaluationDomain<Fr>,
    values: &[Fr],
    powers: &Powers<Bls12_381>,
) -> OpeningsProof {
    let mut witness_commits = Vec::with_capacity(values.len());

    for (i, value) in values.iter().enumerate() {
        let point = domain.element(i);
        debug_assert_eq!(
            poly.evaluate(&point),
            *value,
            "多項式在 domain 點 {i} 的值跟宣稱的公開值不符"
        );

        let witness_poly = kzg10::compute_witness(poly, point);
        witness_commits.push(kzg10::commit(powers, &witness_poly));
    }

    OpeningsProof { witness_commits }
}

pub fn verify_openings(
    commitment: &Commitment<Bls12_381>,
    proof: &OpeningsProof,
    domain: &Radix2EvaluationDomain<Fr>,
    values: &[Fr],
    vk: &VerifierKey<Bls12_381>,
) -> bool {
    if proof.witness_commits.len() != values.len() {
        return false;
    }

    for (i, value) in values.iter().enumerate() {
        let point = domain.element(i);

        if !kzg10::verify(vk, commitment, &proof.witness_commits[i], point, *value) {
            return false;
        }
    }

    true
}
