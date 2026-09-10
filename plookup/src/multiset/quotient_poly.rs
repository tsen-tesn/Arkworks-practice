use ark_bls12_381::Fr;
use ark_ff::{One, Zero};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial,
    Radix2EvaluationDomain,
};

/// 主函數：將四條 constraints 合成一個多項式 C(X)，
/// 除以 vanishing polynomial Z_H(X) = X^n - 1，
/// 回傳 (quotient Q(X), remainder)。
/// 若 witness 合法，remainder 必須為零多項式。
pub fn compute(
    domain: &Radix2EvaluationDomain<Fr>,
    z_poly: &DensePolynomial<Fr>,
    f_poly: &DensePolynomial<Fr>,
    t_poly: &DensePolynomial<Fr>,
    h_1_poly: &DensePolynomial<Fr>,
    h_2_poly: &DensePolynomial<Fr>,
    beta: Fr,
    gamma: Fr,
) -> (DensePolynomial<Fr>, DensePolynomial<Fr>) {
    let point_check = compute_point_checks(z_poly, domain);
    let interval_check = compute_interval_check(h_1_poly, h_2_poly, domain);
    let term_check = compute_term_check(
        domain, z_poly, f_poly, t_poly, h_1_poly, h_2_poly, beta, gamma,
    );
    let sum = &(&interval_check + &point_check) + &term_check;
    sum.divide_by_vanishing_poly(*domain)
}

/// C1 + C4：合併起點與終點的 Z 邊界條件。
/// (L_0 + L_n) * (Z(X) - 1)
fn compute_point_checks(
    z_poly: &DensePolynomial<Fr>,
    domain: &Radix2EvaluationDomain<Fr>,
) -> DensePolynomial<Fr> {
    let l0_poly = compute_n_lagrange_poly(domain, 0);
    let ln_poly = compute_n_lagrange_poly(domain, domain.size() - 1);
    // Z(X) - 1
    let z_prime_poly = z_poly - &DensePolynomial::from_coefficients_vec(vec![Fr::one()]);
    // (L_0 + L_n)
    let l_poly = &l0_poly + &ln_poly;
    &z_prime_poly * &l_poly
}

/// C3：h₁/h₂ 銜接條件。
/// L_n(X) * (h₁(X) - h₂(ωX))
///
/// 用 2n domain 來計算「ωX 的 shift」：
/// 在 2n domain 中，原 n domain 的 generator ω 對應步進 2，
/// 所以 h₂(ωX) 在 index i 處等於 h_2_evals[i + 2]。
fn compute_interval_check(
    h_1_poly: &DensePolynomial<Fr>,
    h_2_poly: &DensePolynomial<Fr>,
    domain: &Radix2EvaluationDomain<Fr>,
) -> DensePolynomial<Fr> {
    let domain_2n: Radix2EvaluationDomain<Fr> =
        EvaluationDomain::new(2 * domain.size()).unwrap();

    // 把 L_n 的 evaluations（在 n domain）升到 2n domain 上
    let ln_evals = compute_n_lagrange_evaluations(domain.size(), domain.size() - 1);
    let ln_2n_evals = domain_2n.fft(&domain.ifft(&ln_evals));

    let h_1_evals = domain_2n.fft(&h_1_poly.coeffs);

    // h_2 在 2n domain 上的 evaluations；
    // 補 2 個 wrap-around 讓 i+2 的索引安全
    let mut h_2_evals = domain_2n.fft(&h_2_poly.coeffs);
    h_2_evals.push(h_2_evals[0]);
    h_2_evals.push(h_2_evals[1]);

    let i_evals: Vec<Fr> = (0..domain_2n.size())
        .map(|i| {
            let ln_i = ln_2n_evals[i];
            let h_1_i = h_1_evals[i];
            let h_2_i_next = h_2_evals[i + 2]; // h₂(ωX)
            ln_i * (h_1_i - h_2_i_next)
        })
        .collect();

    DensePolynomial::from_coefficients_vec(domain_2n.ifft(&i_evals))
}

/// C2：grand-product transition。
/// (X - ωⁿ) * [Z(X) * (1+β) * F(X) - Z(ωX) * G(X)] = 0
///
/// 分成 a（分子項）和 b（分母項）分別計算，再相減。
pub fn compute_term_check(
    domain: &Radix2EvaluationDomain<Fr>,
    z_poly: &DensePolynomial<Fr>,
    f_poly: &DensePolynomial<Fr>,
    t_poly: &DensePolynomial<Fr>,
    h_1_poly: &DensePolynomial<Fr>,
    h_2_poly: &DensePolynomial<Fr>,
    beta: Fr,
    gamma: Fr,
) -> DensePolynomial<Fr> {
    let part_a = compute_term_check_a(domain, z_poly, f_poly, t_poly, beta, gamma);
    let part_b = compute_term_check_b(domain, z_poly, h_1_poly, h_2_poly, beta, gamma);
    &part_a - &part_b
}

/// 分子側：(X - ωⁿ) * Z(X) * (1+β) * (γ + f(X)) * (γ(1+β) + t(X) + β·t(ωX))
///
/// 用 4n domain 計算，讓 t(ωX) 的步進 +4 對應原 n domain 的 ω 步進。
fn compute_term_check_a(
    domain: &Radix2EvaluationDomain<Fr>,
    z_poly: &DensePolynomial<Fr>,
    f_poly: &DensePolynomial<Fr>,
    t_poly: &DensePolynomial<Fr>,
    beta: Fr,
    gamma: Fr,
) -> DensePolynomial<Fr> {
    let domain_4n: Radix2EvaluationDomain<Fr> =
        EvaluationDomain::new(4 * domain.size()).unwrap();

    let z_evals = domain_4n.fft(&z_poly.coeffs);
    let f_evals = domain_4n.fft(&f_poly.coeffs);

    // t 補 4 個 wrap-around 讓 i+4 的索引安全
    let mut t_evals = domain_4n.fft(&t_poly.coeffs);
    t_evals.push(t_evals[0]);
    t_evals.push(t_evals[1]);
    t_evals.push(t_evals[2]);
    t_evals.push(t_evals[3]);

    let beta_one = Fr::one() + beta;
    let g_n = domain.elements().last().unwrap(); // ωⁿ（最後一個 domain 元素）

    let i_evals: Vec<Fr> = (0..domain_4n.size())
        .zip(domain_4n.elements())
        .map(|(i, root_i)| {
            let z_i = z_evals[i];
            let f_i = f_evals[i];
            let t_i = t_evals[i];
            let t_i_next = t_evals[i + 4]; // t(ωX)
            let a = root_i - g_n;          // (X - ωⁿ)：在最後一點歸零
            let b = z_i * beta_one;
            let c = gamma + f_i;
            let d = (gamma * beta_one) + t_i + (beta * t_i_next);
            a * b * c * d
        })
        .collect();

    let i_poly = DensePolynomial::from_coefficients_vec(domain_4n.ifft(&i_evals));
    // 在 ωⁿ 處必須為 0（(X - ωⁿ) 因子保證）
    assert_eq!(i_poly.evaluate(&g_n), Fr::zero());
    i_poly
}

/// 分母側：(X - ωⁿ) * Z(ωX) * (γ(1+β) + h₁(X) + β·h₁(ωX)) * (γ(1+β) + h₂(X) + β·h₂(ωX))
fn compute_term_check_b(
    domain: &Radix2EvaluationDomain<Fr>,
    z_poly: &DensePolynomial<Fr>,
    h_1_poly: &DensePolynomial<Fr>,
    h_2_poly: &DensePolynomial<Fr>,
    beta: Fr,
    gamma: Fr,
) -> DensePolynomial<Fr> {
    let domain_4n: Radix2EvaluationDomain<Fr> =
        EvaluationDomain::new(4 * domain.size()).unwrap();

    // z, h_1, h_2 各補 4 個 wrap-around 讓 i+4 安全
    let mut z_evals = domain_4n.fft(&z_poly.coeffs);
    z_evals.push(z_evals[0]);
    z_evals.push(z_evals[1]);
    z_evals.push(z_evals[2]);
    z_evals.push(z_evals[3]);

    let mut h_1_evals = domain_4n.fft(&h_1_poly.coeffs);
    h_1_evals.push(h_1_evals[0]);
    h_1_evals.push(h_1_evals[1]);
    h_1_evals.push(h_1_evals[2]);
    h_1_evals.push(h_1_evals[3]);

    let mut h_2_evals = domain_4n.fft(&h_2_poly.coeffs);
    h_2_evals.push(h_2_evals[0]);
    h_2_evals.push(h_2_evals[1]);
    h_2_evals.push(h_2_evals[2]);
    h_2_evals.push(h_2_evals[3]);

    let beta_one = Fr::one() + beta;
    let g_n = domain.elements().last().unwrap();

    let i_evals: Vec<Fr> = (0..domain_4n.size())
        .zip(domain_4n.elements())
        .map(|(i, root_i)| {
            let z_i_next = z_evals[i + 4];   // Z(ωX)
            let h_1_i = h_1_evals[i];
            let h_1_i_next = h_1_evals[i + 4]; // h₁(ωX)
            let h_2_i = h_2_evals[i];
            let h_2_i_next = h_2_evals[i + 4]; // h₂(ωX)
            let a = (root_i - g_n) * z_i_next;
            let b = (gamma * beta_one) + h_1_i + (beta * h_1_i_next);
            let c = (gamma * beta_one) + h_2_i + (beta * h_2_i_next);
            a * b * c
        })
        .collect();

    let i_poly = DensePolynomial::from_coefficients_vec(domain_4n.ifft(&i_evals));
    assert_eq!(i_poly.evaluate(&g_n), Fr::zero());
    i_poly
}

/// Lagrange basis polynomial L_n(X)：在 domain[n] 處等於 1，其他點等於 0。
pub fn compute_n_lagrange_poly(
    domain: &Radix2EvaluationDomain<Fr>,
    n: usize,
) -> DensePolynomial<Fr> {
    assert!(n <= domain.size() - 1);
    let mut evaluations = compute_n_lagrange_evaluations(domain.size(), n);
    domain.ifft_in_place(&mut evaluations);
    DensePolynomial::from_coefficients_vec(evaluations)
}

fn compute_n_lagrange_evaluations(domain_size: usize, n: usize) -> Vec<Fr> {
    let mut evals = vec![Fr::zero(); domain_size];
    evals[n] = Fr::one();
    evals
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_poly::EvaluationDomain;
    use crate::multiset::multiset::MultiSet;
    use crate::multiset::multiset_equality::{compute_accumulator_values, compute_h1_h2};
    use crate::lookup::lookup::LookUp;
    use crate::lookup::proof::{compress_column, encode_and_pad_witness};
    use crate::lookup::table::relu::ReLUTable;

    /// 小型驗煙測試（對齊參考 repo 的 test_quotient_poly）：
    /// f = [2,3,4]，t = [2,3,4,5]，合法 witness，remainder 必須為零。
    #[test]
    fn quotient_remainder_is_zero_for_valid_witness() {
        let f = MultiSet::from_slice(&[2u64, 3, 4].map(Fr::from));
        let t = MultiSet::from_slice(&[2u64, 3, 4, 5].map(Fr::from));

        let domain: Radix2EvaluationDomain<Fr> = EvaluationDomain::new(f.len()).unwrap();
        let beta = Fr::from(10u64);
        let gamma = Fr::from(11u64);

        let (h_1, h_2) = compute_h1_h2(&f, &t);
        let f_poly = f.to_polynomial(&domain);             
        let t_poly = t.to_polynomial(&domain);
        let h_1_poly = h_1.to_polynomial(&domain);
        let h_2_poly = h_2.to_polynomial(&domain);

        let z_evaluations = compute_accumulator_values(&f, &t, &h_1, &h_2, beta, gamma);
        let z_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&z_evaluations));

        let (_, remainder) = compute(
            &domain, &z_poly, &f_poly, &t_poly, &h_1_poly, &h_2_poly, beta, gamma,
        );
        assert!(remainder.is_zero());
    }

    /// ReLU 整數案例：合法 witness，四條 constraints 全部為零，remainder 為零。
    #[test]
    fn quotient_remainder_is_zero_for_relu_witness() {
        let table = ReLUTable::new();
        let mut lookup = LookUp::new(ReLUTable::new());
        for input in [2, 4, -1, 0, 15, 6, -24, -2, 15] {
            assert!(lookup.read(input));
        }
        let (encoded_inputs, encoded_outputs) =
            encode_and_pad_witness(lookup.input_wires(), lookup.output_wires());
        let theta = Fr::from(7u64);
        let compressed_witness =
            compress_column(&encoded_inputs, &encoded_outputs, theta);
        let compressed_table = compress_column(
            &table.encode_input_column(),
            &table.encode_output_column(),
            theta,
        );

        let f = MultiSet::from_slice(&compressed_witness);
        let t = MultiSet::from_slice(&compressed_table);

        let domain: Radix2EvaluationDomain<Fr> = EvaluationDomain::new(t.len()).unwrap();
        let beta = Fr::from(5u64);
        let gamma = Fr::from(6u64);

        let (h_1, h_2) = compute_h1_h2(&f, &t);
        let f_poly = f.to_polynomial(&domain);
        let t_poly = t.to_polynomial(&domain);
        let h_1_poly = h_1.to_polynomial(&domain);
        let h_2_poly = h_2.to_polynomial(&domain);

        let z_evaluations = compute_accumulator_values(&f, &t, &h_1, &h_2, beta, gamma);
        let z_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&z_evaluations));

        let (_, remainder) = compute(
            &domain, &z_poly, &f_poly, &t_poly, &h_1_poly, &h_2_poly, beta, gamma,
        );
        assert!(remainder.is_zero());
    }
}
