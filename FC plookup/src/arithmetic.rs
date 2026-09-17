use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::Zero;
use ark_poly::{
    polynomial::univariate::DensePolynomial, EvaluationDomain, Polynomial, Radix2EvaluationDomain,
};
use ark_poly_commit::kzg10::{Commitment, Powers, VerifierKey};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use plookup::{kzg10, transcript::TranscriptProtocol};

use crate::columns::poly_from_rows;

/// 這個 crate 只處理 3x3 全連接層： z = W·X + b。
pub const N_ROWS: usize = 3;

/// W 的每一欄 w_j、X 的每個分量 x_j 都各自是一條跨 3 個 row 的多項式，
/// b(X)、z(X) 也是同樣的 row-indexed 多項式。
/// z(X) 就是直接拿去餵給 Plookup witness column 的同一顆多項式，不會另外重算。
pub struct ProverInput {
    /// row-major：w[i][j] 是 W 的第 i 列第 j 行
    pub w: [[Fr; N_ROWS]; N_ROWS],
    pub x: [Fr; N_ROWS],
    pub b: [Fr; N_ROWS],
    pub z: [Fr; N_ROWS],
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Commitments {
    pub w: [Commitment<Bls12_381>; N_ROWS],
    pub x: [Commitment<Bls12_381>; N_ROWS],
    pub b: Commitment<Bls12_381>,
    pub z: Commitment<Bls12_381>,
    pub q: Commitment<Bls12_381>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Evaluations {
    pub w: [Fr; N_ROWS],
    pub x: [Fr; N_ROWS],
    pub b: Fr,
    pub z: Fr,
    pub q: Fr,
}

/// 直接對 z(X)（以及 W, X, b 的委諾）證明 z = W·X + b，
/// 手法跟 plookup 的 quotient_poly.rs / multiset_equality.rs 一致：
/// 把等式差建成一個在整個 domain 上都要為零的多項式，除以 vanishing polynomial 檢查整除性，
/// 再用一個隨機挑戰點做批次 KZG 開啟證明。
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct ArithmeticProof {
    pub commitments: Commitments,
    pub evaluations: Evaluations,
    pub aggregate_witness_comm: Commitment<Bls12_381>,
}

pub fn prove(
    input: &ProverInput,
    domain: &Radix2EvaluationDomain<Fr>,
    powers: &Powers<Bls12_381>,
    transcript: &mut dyn TranscriptProtocol,
) -> (ArithmeticProof, DensePolynomial<Fr>) {
    let w_col_polys: [DensePolynomial<Fr>; N_ROWS] = std::array::from_fn(|j| {
        let column = [input.w[0][j], input.w[1][j], input.w[2][j]];
        poly_from_rows(domain, &column)
    });
    let x_polys: [DensePolynomial<Fr>; N_ROWS] =
        std::array::from_fn(|j| poly_from_rows(domain, &[input.x[j]; N_ROWS]));
    let b_poly = poly_from_rows(domain, &input.b);
    let z_poly = poly_from_rows(domain, &input.z);

    // A(X) = z(X) - [ w_0(X)x_0(X) + w_1(X)x_1(X) + w_2(X)x_2(X) + b(X) ]
    let mut rhs_poly = &w_col_polys[0] * &x_polys[0];
    for j in 1..N_ROWS {
        rhs_poly = &rhs_poly + &(&w_col_polys[j] * &x_polys[j]);
    }
    rhs_poly = &rhs_poly + &b_poly;
    let a_poly = &z_poly - &rhs_poly;

    let (q_poly, remainder) = a_poly.divide_by_vanishing_poly(*domain);
    assert!(
        remainder.is_zero(),
        "z = W·X + b 不成立：constraint polynomial 沒辦法被 Z_H(X) 整除"
    );

    let w_commits: [Commitment<Bls12_381>; N_ROWS] =
        std::array::from_fn(|j| kzg10::commit(powers, &w_col_polys[j]));
    let x_commits: [Commitment<Bls12_381>; N_ROWS] =
        std::array::from_fn(|j| kzg10::commit(powers, &x_polys[j]));
    let b_commit = kzg10::commit(powers, &b_poly);
    let z_commit = kzg10::commit(powers, &z_poly);
    let q_commit = kzg10::commit(powers, &q_poly);

    for commit in w_commits.iter() {
        transcript.append_commitment(b"arith_w", commit);
    }
    for commit in x_commits.iter() {
        transcript.append_commitment(b"arith_x", commit);
    }
    transcript.append_commitment(b"arith_b", &b_commit);
    transcript.append_commitment(b"arith_z", &z_commit);
    transcript.append_commitment(b"arith_q", &q_commit);

    let evaluation_challenge = transcript.challenge_scalar(b"arith_evaluation_challenge");

    let w_evals: [Fr; N_ROWS] =
        std::array::from_fn(|j| w_col_polys[j].evaluate(&evaluation_challenge));
    let x_evals: [Fr; N_ROWS] = std::array::from_fn(|j| x_polys[j].evaluate(&evaluation_challenge));
    let b_eval = b_poly.evaluate(&evaluation_challenge);
    let z_eval = z_poly.evaluate(&evaluation_challenge);
    let q_eval = q_poly.evaluate(&evaluation_challenge);

    for eval in w_evals.iter() {
        transcript.append_scalar(b"arith_w_eval", eval);
    }
    for eval in x_evals.iter() {
        transcript.append_scalar(b"arith_x_eval", eval);
    }
    transcript.append_scalar(b"arith_b_eval", &b_eval);
    transcript.append_scalar(b"arith_z_eval", &z_eval);
    transcript.append_scalar(b"arith_q_eval", &q_eval);

    let aggregation_challenge = transcript.challenge_scalar(b"arith_aggregation");

    let mut ordered_polys: Vec<&DensePolynomial<Fr>> = Vec::with_capacity(2 * N_ROWS + 3);
    ordered_polys.extend(w_col_polys.iter());
    ordered_polys.extend(x_polys.iter());
    ordered_polys.push(&b_poly);
    ordered_polys.push(&z_poly);
    ordered_polys.push(&q_poly);

    let agg_witness =
        kzg10::compute_aggregate_witness(ordered_polys, evaluation_challenge, aggregation_challenge);
    let aggregate_witness_comm = kzg10::commit(powers, &agg_witness);

    let proof = ArithmeticProof {
        commitments: Commitments {
            w: w_commits,
            x: x_commits,
            b: b_commit,
            z: z_commit,
            q: q_commit,
        },
        evaluations: Evaluations {
            w: w_evals,
            x: x_evals,
            b: b_eval,
            z: z_eval,
            q: q_eval,
        },
        aggregate_witness_comm,
    };

    (proof, z_poly)
}

pub fn verify(
    proof: &ArithmeticProof,
    domain: &Radix2EvaluationDomain<Fr>,
    vk: &VerifierKey<Bls12_381>,
    transcript: &mut dyn TranscriptProtocol,
) -> bool {
    for commit in proof.commitments.w.iter() {
        transcript.append_commitment(b"arith_w", commit);
    }
    for commit in proof.commitments.x.iter() {
        transcript.append_commitment(b"arith_x", commit);
    }
    transcript.append_commitment(b"arith_b", &proof.commitments.b);
    transcript.append_commitment(b"arith_z", &proof.commitments.z);
    transcript.append_commitment(b"arith_q", &proof.commitments.q);

    let evaluation_challenge = transcript.challenge_scalar(b"arith_evaluation_challenge");

    for eval in proof.evaluations.w.iter() {
        transcript.append_scalar(b"arith_w_eval", eval);
    }
    for eval in proof.evaluations.x.iter() {
        transcript.append_scalar(b"arith_x_eval", eval);
    }
    transcript.append_scalar(b"arith_b_eval", &proof.evaluations.b);
    transcript.append_scalar(b"arith_z_eval", &proof.evaluations.z);
    transcript.append_scalar(b"arith_q_eval", &proof.evaluations.q);

    let aggregation_challenge = transcript.challenge_scalar(b"arith_aggregation");

    // z(ζ) - [Σ w_j(ζ)x_j(ζ) + b(ζ)] 必須等於 Q(ζ) * Z_H(ζ)
    let mut rhs_eval = Fr::zero();
    for j in 0..N_ROWS {
        rhs_eval += proof.evaluations.w[j] * proof.evaluations.x[j];
    }
    rhs_eval += proof.evaluations.b;
    let a_eval = proof.evaluations.z - rhs_eval;
    let vanishing_eval = domain.evaluate_vanishing_polynomial(evaluation_challenge);

    if a_eval != proof.evaluations.q * vanishing_eval {
        return false;
    }

    let mut ordered_commitments: Vec<&Commitment<Bls12_381>> = Vec::with_capacity(2 * N_ROWS + 3);
    ordered_commitments.extend(proof.commitments.w.iter());
    ordered_commitments.extend(proof.commitments.x.iter());
    ordered_commitments.push(&proof.commitments.b);
    ordered_commitments.push(&proof.commitments.z);
    ordered_commitments.push(&proof.commitments.q);

    let mut ordered_values: Vec<&Fr> = Vec::with_capacity(2 * N_ROWS + 3);
    ordered_values.extend(proof.evaluations.w.iter());
    ordered_values.extend(proof.evaluations.x.iter());
    ordered_values.push(&proof.evaluations.b);
    ordered_values.push(&proof.evaluations.z);
    ordered_values.push(&proof.evaluations.q);

    let agg_commitment = kzg10::aggregate_commitments(ordered_commitments, aggregation_challenge);
    let agg_value = kzg10::aggregate_values(ordered_values, aggregation_challenge);

    kzg10::verify(
        vk,
        &agg_commitment,
        &proof.aggregate_witness_comm,
        evaluation_challenge,
        agg_value,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use merlin::Transcript;
    use plookup::lookup::proof::TABLE_SIZE;

    fn f(value: i64) -> Fr {
        if value >= 0 {
            Fr::from(value as u64)
        } else {
            -Fr::from(value.unsigned_abs())
        }
    }

    // W = [[1,2,-1],[0,1,1],[-1,0,2]], X = [3,-2,1], b = [1,-1,2]
    // z1 = 3 - 4 - 1 + 1 = -1, z2 = 0 - 2 + 1 - 1 = -2, z3 = -3 + 0 + 2 + 2 = 1
    fn valid_input() -> ProverInput {
        ProverInput {
            w: [
                [f(1), f(2), f(-1)],
                [f(0), f(1), f(1)],
                [f(-1), f(0), f(2)],
            ],
            x: [f(3), f(-2), f(1)],
            b: [f(1), f(-1), f(2)],
            z: [f(-1), f(-2), f(1)],
        }
    }

    #[test]
    fn accepts_correct_z() {
        let (powers, vk) = kzg10::trusted_setup(256);
        let domain: Radix2EvaluationDomain<Fr> = EvaluationDomain::new(TABLE_SIZE).unwrap();

        let mut prover_transcript = Transcript::new(b"arithmetic-test");
        let (proof, _z_poly) = prove(&valid_input(), &domain, &powers, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"arithmetic-test");
        assert!(verify(&proof, &domain, &vk, &mut verifier_transcript));
    }

    #[test]
    #[should_panic(expected = "沒辦法被 Z_H(X) 整除")]
    fn rejects_wrong_z_at_prove_time() {
        let (powers, _vk) = kzg10::trusted_setup(256);
        let domain: Radix2EvaluationDomain<Fr> = EvaluationDomain::new(TABLE_SIZE).unwrap();

        let mut input = valid_input();
        input.z[0] = f(0); // 故意填錯

        let mut prover_transcript = Transcript::new(b"arithmetic-test");
        prove(&input, &domain, &powers, &mut prover_transcript);
    }

    #[test]
    fn verify_rejects_tampered_evaluation() {
        let (powers, vk) = kzg10::trusted_setup(256);
        let domain: Radix2EvaluationDomain<Fr> = EvaluationDomain::new(TABLE_SIZE).unwrap();

        let mut prover_transcript = Transcript::new(b"arithmetic-test");
        let (mut proof, _z_poly) = prove(&valid_input(), &domain, &powers, &mut prover_transcript);

        // 竄改其中一個 evaluation，quotient 等式應該對不上
        proof.evaluations.b += Fr::from(1u64);

        let mut verifier_transcript = Transcript::new(b"arithmetic-test");
        assert!(!verify(&proof, &domain, &vk, &mut verifier_transcript));
    }
}
