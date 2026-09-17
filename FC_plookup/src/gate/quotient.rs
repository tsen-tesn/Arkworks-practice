use crate::gate::layout::{GateLayout, GateType, NUM_GATES};
use ark_bls12_381::Fr;
use ark_ff::{One, Zero};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
};

/// 依照 layout 產生五條 selector 的 evaluation 向量（長度 32，跟 w,x,b 完全無關）。
///
/// 乘法 gate： q_M=1, q_O=-1，其餘為 0   -> 恆等式 a*b - c = 0
/// 加法 gate： q_L=1, q_R=1, q_O=-1     -> 恆等式 a+b - c = 0
pub fn build_selector_evaluations(
    layout: &[GateLayout; NUM_GATES],
) -> (Vec<Fr>, Vec<Fr>, Vec<Fr>, Vec<Fr>, Vec<Fr>) {
    let mut q_m = Vec::with_capacity(NUM_GATES);
    let mut q_l = Vec::with_capacity(NUM_GATES);
    let mut q_r = Vec::with_capacity(NUM_GATES);
    let mut q_o = Vec::with_capacity(NUM_GATES);
    let mut q_c = Vec::with_capacity(NUM_GATES);

    for gate in layout.iter() {
        match gate.gate_type {
            GateType::Mul => {
                q_m.push(Fr::one());
                q_l.push(Fr::zero());
                q_r.push(Fr::zero());
                q_o.push(-Fr::one());
                q_c.push(Fr::zero());
            }
            GateType::Add => {
                q_m.push(Fr::zero());
                q_l.push(Fr::one());
                q_r.push(Fr::one());
                q_o.push(-Fr::one());
                q_c.push(Fr::zero());
            }
        }
    }

    (q_m, q_l, q_r, q_o, q_c)
}

fn to_poly(domain: &Radix2EvaluationDomain<Fr>, evals: &[Fr]) -> DensePolynomial<Fr> {
    DensePolynomial::from_coefficients_vec(domain.ifft(evals))
}

/// GateCheck(X) = q_M(X)a(X)b(X) + q_L(X)a(X) + q_R(X)b(X) + q_O(X)c(X) + q_C(X)
///
/// qM*a*b 這一項最高到三個 degree (N-1) 多項式相乘，係數長度可能到 3N-2，
/// 所以要在夠大的 domain（4N，比照 plookup/multiset/quotient_poly.rs 的作法）上
/// 逐點相乘再 IFFT 回係數，而不是直接對 DensePolynomial 做乘法。
fn compute_gate_check_poly(
    domain: &Radix2EvaluationDomain<Fr>,
    q_m: &DensePolynomial<Fr>,
    q_l: &DensePolynomial<Fr>,
    q_r: &DensePolynomial<Fr>,
    q_o: &DensePolynomial<Fr>,
    q_c: &DensePolynomial<Fr>,
    a: &DensePolynomial<Fr>,
    b: &DensePolynomial<Fr>,
    c: &DensePolynomial<Fr>,
) -> DensePolynomial<Fr> {
    let domain_4n: Radix2EvaluationDomain<Fr> = EvaluationDomain::new(4 * domain.size()).unwrap();

    let q_m_e = domain_4n.fft(&q_m.coeffs);
    let q_l_e = domain_4n.fft(&q_l.coeffs);
    let q_r_e = domain_4n.fft(&q_r.coeffs);
    let q_o_e = domain_4n.fft(&q_o.coeffs);
    let q_c_e = domain_4n.fft(&q_c.coeffs);
    let a_e = domain_4n.fft(&a.coeffs);
    let b_e = domain_4n.fft(&b.coeffs);
    let c_e = domain_4n.fft(&c.coeffs);

    let gate_check_evals: Vec<Fr> = (0..domain_4n.size())
        .map(|i| {
            q_m_e[i] * a_e[i] * b_e[i] + q_l_e[i] * a_e[i] + q_r_e[i] * b_e[i] + q_o_e[i] * c_e[i]
                + q_c_e[i]
        })
        .collect();

    DensePolynomial::from_coefficients_vec(domain_4n.ifft(&gate_check_evals))
}

/// 主函數：把 gate identity 合成一個多項式，除以 vanishing polynomial Z_H(X) = X^N - 1，
/// 回傳 (quotient Q(X), remainder)。若 witness 合法，remainder 必須為零多項式。
pub fn compute(
    domain: &Radix2EvaluationDomain<Fr>,
    a_evals: &[Fr],
    b_evals: &[Fr],
    c_evals: &[Fr],
    q_m_evals: &[Fr],
    q_l_evals: &[Fr],
    q_r_evals: &[Fr],
    q_o_evals: &[Fr],
    q_c_evals: &[Fr],
) -> (DensePolynomial<Fr>, DensePolynomial<Fr>) {
    let a = to_poly(domain, a_evals);
    let b = to_poly(domain, b_evals);
    let c = to_poly(domain, c_evals);
    let q_m = to_poly(domain, q_m_evals);
    let q_l = to_poly(domain, q_l_evals);
    let q_r = to_poly(domain, q_r_evals);
    let q_o = to_poly(domain, q_o_evals);
    let q_c = to_poly(domain, q_c_evals);

    let gate_check = compute_gate_check_poly(domain, &q_m, &q_l, &q_r, &q_o, &q_c, &a, &b, &c);
    gate_check.divide_by_vanishing_poly(*domain)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gate::layout::{build_gate_layout, compute_wires};
    use ark_poly::Radix2EvaluationDomain;

    fn f(v: u64) -> Fr {
        Fr::from(v)
    }

    fn sample_domain() -> Radix2EvaluationDomain<Fr> {
        EvaluationDomain::new(NUM_GATES).unwrap()
    }

    #[test]
    fn valid_witness_gives_zero_remainder() {
        let layout = build_gate_layout();
        let w = [
            [f(1), f(2), f(3)],
            [f(4), f(5), f(6)],
            [f(7), f(8), f(9)],
        ];
        let x = [f(1), f(1), f(1)];
        let b = [f(1), f(1), f(1)];

        let (a_evals, b_evals, c_evals) = compute_wires(&layout, &w, &x, &b);
        let (q_m, q_l, q_r, q_o, q_c) = build_selector_evaluations(&layout);
        let domain = sample_domain();

        let (_quotient, remainder) = compute(
            &domain, &a_evals, &b_evals, &c_evals, &q_m, &q_l, &q_r, &q_o, &q_c,
        );

        assert!(remainder.is_zero());
    }

    #[test]
    fn corrupted_witness_gives_nonzero_remainder() {
        let layout = build_gate_layout();
        let w = [
            [f(1), f(2), f(3)],
            [f(4), f(5), f(6)],
            [f(7), f(8), f(9)],
        ];
        let x = [f(1), f(1), f(1)];
        let b = [f(1), f(1), f(1)];

        let (a_evals, b_evals, mut c_evals) = compute_wires(&layout, &w, &x, &b);
        c_evals[5] += f(1); // 破壞其中一個 gate 的輸出，使該 gate 的恆等式不成立

        let (q_m, q_l, q_r, q_o, q_c) = build_selector_evaluations(&layout);
        let domain = sample_domain();

        let (_quotient, remainder) = compute(
            &domain, &a_evals, &b_evals, &c_evals, &q_m, &q_l, &q_r, &q_o, &q_c,
        );

        assert!(!remainder.is_zero());
    }
}
