use crate::gate::layout::{
    build_gate_layout, build_z_copy_column, compute_wires, extract_z_outputs, MATRIX_DIM,
    NUM_GATES,
};
use crate::gate::permutation;
use crate::gate::quotient;
use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::One;
use ark_poly::{
    univariate::DensePolynomial as Polynomial, DenseUVPolynomial, EvaluationDomain, Polynomial as _,
    Radix2EvaluationDomain,
};
use ark_poly_commit::kzg10::{Commitment, Powers, VerifierKey};
use plookup::kzg10;
use plookup::transcript::TranscriptProtocol;

/// 電路四條 wire (a,b,c,d) 跟 permutation accumulator (z) 在挑戰點的開值。
/// z_omega 是 z 在 zeta*omega 的開值，transition check 需要它。
pub struct GateEvaluations {
    pub a: Fr,
    pub b: Fr,
    pub c: Fr,
    pub d: Fr,
    pub z: Fr,
    pub z_omega: Fr,
}

pub struct GateCommitments {
    pub a: Commitment<Bls12_381>,
    pub b: Commitment<Bls12_381>,
    pub c: Commitment<Bls12_381>,
    pub d: Commitment<Bls12_381>,
    pub z: Commitment<Bls12_381>,
    pub q1: Commitment<Bls12_381>, // gate identity（M2）的 quotient commitment
    pub q2: Commitment<Bls12_381>, // permutation（M3）的 quotient commitment
}

/// 證明「我知道私密的 w, x, b，使得 z = Wx+b 算術正確、wiring 正確」，
/// 並且額外公開 commit 了一條 D 欄位（透過 permutation 強制等於 z），
/// 供 caller（fc_proof.rs）拿去跟 plookup 的 lookup 論證做同態綁定（M4）。
///
/// z 本身全程不會以明文形式出現在這個 proof 結構裡：commitments 只是群元素，
/// evaluations 是在一個 Fiat-Shamir 隨機挑戰點的開值，不是電路裡任何一個
/// 有意義位置（例如 z 所在的 gate 11/14/17）的值。
pub struct GateProof {
    pub commitments: GateCommitments,
    pub evaluations: GateEvaluations,
    pub aggregate_witness_comm: Commitment<Bls12_381>,
    pub shifted_aggregate_witness_comm: Commitment<Bls12_381>,
}

pub fn domain() -> Radix2EvaluationDomain<Fr> {
    EvaluationDomain::new(NUM_GATES).unwrap()
}

fn eval_public_poly(domain: &Radix2EvaluationDomain<Fr>, evals: &[Fr], point: Fr) -> Fr {
    Polynomial::from_coefficients_vec(domain.ifft(evals)).evaluate(&point)
}

pub fn prove(
    w: &[[Fr; MATRIX_DIM]; MATRIX_DIM],
    x: &[Fr; MATRIX_DIM],
    b: &[Fr; MATRIX_DIM],
    proving_key: &Powers<Bls12_381>,
    transcript: &mut dyn TranscriptProtocol,
) -> (GateProof, [Fr; MATRIX_DIM]) {
    let domain = domain();
    let layout = build_gate_layout();
    let (a_evals, b_evals, c_evals) = compute_wires(&layout, w, x, b);
    let d_evals = build_z_copy_column(&c_evals);
    let z_values = extract_z_outputs(&c_evals);

    let a_poly = Polynomial::from_coefficients_vec(domain.ifft(&a_evals));
    let b_poly = Polynomial::from_coefficients_vec(domain.ifft(&b_evals));
    let c_poly = Polynomial::from_coefficients_vec(domain.ifft(&c_evals));
    let d_poly = Polynomial::from_coefficients_vec(domain.ifft(&d_evals));

    let a_commit = kzg10::commit(proving_key, &a_poly);
    let b_commit = kzg10::commit(proving_key, &b_poly);
    let c_commit = kzg10::commit(proving_key, &c_poly);
    let d_commit = kzg10::commit(proving_key, &d_poly);

    transcript.append_commitment(b"gate_a", &a_commit);
    transcript.append_commitment(b"gate_b", &b_commit);
    transcript.append_commitment(b"gate_c", &c_commit);
    transcript.append_commitment(b"gate_d", &d_commit);

    let beta = transcript.challenge_scalar(b"perm_beta");
    let gamma = transcript.challenge_scalar(b"perm_gamma");

    let id_e = permutation::build_id_evaluations(&domain);
    let sigma_e = permutation::build_sigma_evaluations(&domain);
    let wires_evals: [&[Fr]; permutation::NUM_COLUMNS] = [&a_evals, &b_evals, &c_evals, &d_evals];
    let z_evals =
        permutation::compute_accumulator_values(&wires_evals, &id_e, &sigma_e, beta, gamma);
    let z_poly = Polynomial::from_coefficients_vec(domain.ifft(&z_evals));
    let z_commit = kzg10::commit(proving_key, &z_poly);
    transcript.append_commitment(b"perm_z", &z_commit);

    let alpha = transcript.challenge_scalar(b"alpha");

    let (q_m, q_l, q_r, q_o, q_c) = quotient::build_selector_evaluations(&layout);
    let (q1_poly, _) = quotient::compute(
        &domain, &a_evals, &b_evals, &c_evals, &q_m, &q_l, &q_r, &q_o, &q_c,
    );
    let (q2_poly, _) =
        permutation::compute(&domain, &a_evals, &b_evals, &c_evals, &d_evals, beta, gamma, alpha);

    let q1_commit = kzg10::commit(proving_key, &q1_poly);
    let q2_commit = kzg10::commit(proving_key, &q2_poly);
    transcript.append_commitment(b"gate_q1", &q1_commit);
    transcript.append_commitment(b"gate_q2", &q2_commit);

    let zeta = transcript.challenge_scalar(b"zeta");
    transcript.append_scalar(b"zeta", &zeta);
    let zeta_omega = zeta * domain.group_gen;

    let a_eval = a_poly.evaluate(&zeta);
    let b_eval = b_poly.evaluate(&zeta);
    let c_eval = c_poly.evaluate(&zeta);
    let d_eval = d_poly.evaluate(&zeta);
    let z_eval = z_poly.evaluate(&zeta);
    let q1_eval = q1_poly.evaluate(&zeta);
    let q2_eval = q2_poly.evaluate(&zeta);
    let z_omega_eval = z_poly.evaluate(&zeta_omega);

    transcript.append_scalar(b"a_eval", &a_eval);
    transcript.append_scalar(b"b_eval", &b_eval);
    transcript.append_scalar(b"c_eval", &c_eval);
    transcript.append_scalar(b"d_eval", &d_eval);
    transcript.append_scalar(b"z_eval", &z_eval);
    transcript.append_scalar(b"q1_eval", &q1_eval);
    transcript.append_scalar(b"q2_eval", &q2_eval);
    transcript.append_scalar(b"z_omega_eval", &z_omega_eval);

    let aggregation_challenge = transcript.challenge_scalar(b"gate_witness_aggregation");

    let agg_witness = kzg10::compute_aggregate_witness(
        vec![&a_poly, &b_poly, &c_poly, &d_poly, &z_poly, &q1_poly, &q2_poly],
        zeta,
        aggregation_challenge,
    );
    let agg_witness_comm = kzg10::commit(proving_key, &agg_witness);

    let shifted_agg_witness =
        kzg10::compute_aggregate_witness(vec![&z_poly], zeta_omega, aggregation_challenge);
    let shifted_agg_witness_comm = kzg10::commit(proving_key, &shifted_agg_witness);

    (
        GateProof {
            commitments: GateCommitments {
                a: a_commit,
                b: b_commit,
                c: c_commit,
                d: d_commit,
                z: z_commit,
                q1: q1_commit,
                q2: q2_commit,
            },
            evaluations: GateEvaluations {
                a: a_eval,
                b: b_eval,
                c: c_eval,
                d: d_eval,
                z: z_eval,
                z_omega: z_omega_eval,
            },
            aggregate_witness_comm: agg_witness_comm,
            shifted_aggregate_witness_comm: shifted_agg_witness_comm,
        },
        z_values,
    )
}

impl GateProof {
    pub fn verify(
        &self,
        verification_key: &VerifierKey<Bls12_381>,
        transcript: &mut dyn TranscriptProtocol,
    ) -> bool {
        let domain = domain();
        let layout = build_gate_layout();

        transcript.append_commitment(b"gate_a", &self.commitments.a);
        transcript.append_commitment(b"gate_b", &self.commitments.b);
        transcript.append_commitment(b"gate_c", &self.commitments.c);
        transcript.append_commitment(b"gate_d", &self.commitments.d);

        let beta = transcript.challenge_scalar(b"perm_beta");
        let gamma = transcript.challenge_scalar(b"perm_gamma");

        transcript.append_commitment(b"perm_z", &self.commitments.z);
        let alpha = transcript.challenge_scalar(b"alpha");

        transcript.append_commitment(b"gate_q1", &self.commitments.q1);
        transcript.append_commitment(b"gate_q2", &self.commitments.q2);

        let zeta = transcript.challenge_scalar(b"zeta");
        transcript.append_scalar(b"zeta", &zeta);
        let zeta_omega = zeta * domain.group_gen;

        let (q1_eval, q2_eval) = self.compute_quotient_evaluations(&layout, &domain, beta, gamma, alpha, zeta);

        transcript.append_scalar(b"a_eval", &self.evaluations.a);
        transcript.append_scalar(b"b_eval", &self.evaluations.b);
        transcript.append_scalar(b"c_eval", &self.evaluations.c);
        transcript.append_scalar(b"d_eval", &self.evaluations.d);
        transcript.append_scalar(b"z_eval", &self.evaluations.z);
        transcript.append_scalar(b"q1_eval", &q1_eval);
        transcript.append_scalar(b"q2_eval", &q2_eval);
        transcript.append_scalar(b"z_omega_eval", &self.evaluations.z_omega);

        let aggregation_challenge = transcript.challenge_scalar(b"gate_witness_aggregation");

        let agg_commitment = kzg10::aggregate_commitments(
            vec![
                &self.commitments.a,
                &self.commitments.b,
                &self.commitments.c,
                &self.commitments.d,
                &self.commitments.z,
                &self.commitments.q1,
                &self.commitments.q2,
            ],
            aggregation_challenge,
        );
        let agg_value = kzg10::aggregate_values(
            vec![
                &self.evaluations.a,
                &self.evaluations.b,
                &self.evaluations.c,
                &self.evaluations.d,
                &self.evaluations.z,
                &q1_eval,
                &q2_eval,
            ],
            aggregation_challenge,
        );

        let shifted_agg_commitment =
            kzg10::aggregate_commitments(vec![&self.commitments.z], aggregation_challenge);
        let shifted_agg_value =
            kzg10::aggregate_values(vec![&self.evaluations.z_omega], aggregation_challenge);

        kzg10::batch_verify(
            verification_key,
            vec![agg_commitment, shifted_agg_commitment],
            vec![self.aggregate_witness_comm, self.shifted_aggregate_witness_comm],
            vec![zeta, zeta_omega],
            vec![agg_value, shifted_agg_value],
        )
    }

    /// Verifier 端從已開值的 a,b,c,d,z,z_omega 加上「公開、可自行重算」的
    /// selector / id / sigma / Lagrange / vanishing 多項式在 zeta 的值，
    /// 反推出 q1(zeta)（gate identity）跟 q2(zeta)（permutation）應該是多少。
    /// 這跟 plookup::multiset::proof::EqualityProof::compute_quotient_evaluation
    /// 是同一種手法：verifier 不需要「相信」prover 講的 q_eval，是自己從其他
    /// 已經被 KZG 開值證明約束住的數字算出來的。
    fn compute_quotient_evaluations(
        &self,
        layout: &[crate::gate::layout::GateLayout; NUM_GATES],
        domain: &Radix2EvaluationDomain<Fr>,
        beta: Fr,
        gamma: Fr,
        alpha: Fr,
        zeta: Fr,
    ) -> (Fr, Fr) {
        let ev = &self.evaluations;
        let vanishing_at_zeta = domain.evaluate_vanishing_polynomial(zeta);

        // --- q1: gate identity ---
        let (q_m, q_l, q_r, q_o, q_c) = quotient::build_selector_evaluations(layout);
        let q_m_z = eval_public_poly(domain, &q_m, zeta);
        let q_l_z = eval_public_poly(domain, &q_l, zeta);
        let q_r_z = eval_public_poly(domain, &q_r, zeta);
        let q_o_z = eval_public_poly(domain, &q_o, zeta);
        let q_c_z = eval_public_poly(domain, &q_c, zeta);

        let gate_check_at_zeta =
            q_m_z * ev.a * ev.b + q_l_z * ev.a + q_r_z * ev.b + q_o_z * ev.c + q_c_z;
        let q1_eval = gate_check_at_zeta / vanishing_at_zeta;

        // --- q2: permutation ---
        let id_e = permutation::build_id_evaluations(domain);
        let sigma_e = permutation::build_sigma_evaluations(domain);
        let wires_at_zeta = [ev.a, ev.b, ev.c, ev.d];

        let mut num = Fr::one();
        let mut den = Fr::one();
        for col in 0..permutation::NUM_COLUMNS {
            let id_col_z = eval_public_poly(domain, &id_e[col], zeta);
            let sigma_col_z = eval_public_poly(domain, &sigma_e[col], zeta);
            num *= wires_at_zeta[col] + beta * id_col_z + gamma;
            den *= wires_at_zeta[col] + beta * sigma_col_z + gamma;
        }
        let transition_at_zeta = ev.z * num - ev.z_omega * den;

        let lagrange_evaluations = domain.evaluate_all_lagrange_coefficients(zeta);
        let l0_z = lagrange_evaluations[0];
        let boundary_at_zeta = l0_z * (ev.z - Fr::one());

        let q2_eval = (transition_at_zeta + alpha * boundary_at_zeta) / vanishing_at_zeta;

        (q1_eval, q2_eval)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use merlin::Transcript;

    fn f(v: u64) -> Fr {
        Fr::from(v)
    }

    fn sample_wxb() -> ([[Fr; MATRIX_DIM]; MATRIX_DIM], [Fr; MATRIX_DIM], [Fr; MATRIX_DIM]) {
        (
            [
                [f(1), f(2), f(3)],
                [f(4), f(5), f(6)],
                [f(7), f(8), f(9)],
            ],
            [f(1), f(1), f(1)],
            [f(1), f(1), f(1)],
        )
    }

    #[test]
    fn valid_proof_verifies_and_returns_correct_z() {
        let (powers, vk) = kzg10::trusted_setup(512);
        let (w, x, b) = sample_wxb();

        let mut prover_transcript = Transcript::new(b"fc-plookup-gate");
        let (proof, z_values) = prove(&w, &x, &b, &powers, &mut prover_transcript);

        assert_eq!(z_values, [f(1 + 2 + 3 + 1), f(4 + 5 + 6 + 1), f(7 + 8 + 9 + 1)]);

        let mut verifier_transcript = Transcript::new(b"fc-plookup-gate");
        assert!(proof.verify(&vk, &mut verifier_transcript));
    }

    #[test]
    fn tampered_evaluation_is_rejected() {
        let (powers, vk) = kzg10::trusted_setup(512);
        let (w, x, b) = sample_wxb();

        let mut prover_transcript = Transcript::new(b"fc-plookup-gate");
        let (mut proof, _z_values) = prove(&w, &x, &b, &powers, &mut prover_transcript);

        proof.evaluations.a += f(1);

        let mut verifier_transcript = Transcript::new(b"fc-plookup-gate");
        assert!(!proof.verify(&vk, &mut verifier_transcript));
    }
}

