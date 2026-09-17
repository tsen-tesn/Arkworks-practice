use ark_bls12_381::Fr;
use ark_ff::Zero;
use ark_poly::{
    polynomial::univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain,
    Radix2EvaluationDomain,
};

/// 把前 `rows.len()` 個 domain 點設成 `rows`，其餘補 0，再用 IFFT 內插成多項式。
/// Padding 慣例跟 `plookup::multiset::MultiSet::to_polynomial` 一致，
/// 這樣同一組數字不管走哪條路徑內插，得到的多項式都完全相同。
pub fn poly_from_rows(domain: &Radix2EvaluationDomain<Fr>, rows: &[Fr]) -> DensePolynomial<Fr> {
    let mut evals = vec![Fr::zero(); domain.size()];
    evals[..rows.len()].copy_from_slice(rows);
    DensePolynomial::from_coefficients_vec(domain.ifft(&evals))
}
