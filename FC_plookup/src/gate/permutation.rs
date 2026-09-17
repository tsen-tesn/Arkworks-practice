use crate::gate::layout::NUM_GATES;
use ark_bls12_381::Fr;
use ark_ff::{One, Zero};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
};

/// 四條 wire column 的代號。A,B,C 是 M1/M2 電路的乘法/加法 wire；
/// D 是專門用來把 z 綁進 plookup 論證的「z 副本」欄位（見 build_copy_constraint_groups）。
/// D 完全不參與 gate identity（M2），只透過這裡的 permutation argument 被强制等於 C 欄位
/// 在 Z_GATE_INDICES 位置的值。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Column {
    A,
    B,
    C,
    D,
}

pub const NUM_COLUMNS: usize = 4;

/// 用來讓四條欄位的「身分標籤」落在彼此不相交的陪集（coset）上，
/// 這是標準 PLONK permutation argument 的作法：
/// id(col, row) = coset_shift(col) * ω^row，只要 1, COSET_B, COSET_C, COSET_D
/// 兩兩不落在同一個 domain 的陪集裡，id 值就不會跨欄位碰撞。
pub const COSET_B: u64 = 5;
pub const COSET_C: u64 = 7;
pub const COSET_D: u64 = 11;

fn coset_shift(col: Column) -> Fr {
    match col {
        Column::A => Fr::one(),
        Column::B => Fr::from(COSET_B),
        Column::C => Fr::from(COSET_C),
        Column::D => Fr::from(COSET_D),
    }
}

fn column_index(col: Column) -> usize {
    match col {
        Column::A => 0,
        Column::B => 1,
        Column::C => 2,
        Column::D => 3,
    }
}

fn column_from_index(idx: usize) -> Column {
    match idx {
        0 => Column::A,
        1 => Column::B,
        2 => Column::C,
        _ => Column::D,
    }
}

/// 所有 copy constraint 群組：同一群組內的 cell 全部視為必須相等。
/// 對應電路裡三種接線需求：
/// 1. 乘法閘輸出 -> 加法閘輸入（M1 的 wiring）
/// 2. x_j 在三個列的乘法閘裡被重複使用
/// 3. z 的 gate 輸出（C 欄位）-> z 的專用 copy 欄位（D 欄位），把 z 綁進 lookup 論證
pub fn build_copy_constraint_groups() -> Vec<Vec<(Column, usize)>> {
    use Column::*;
    vec![
        // 乘法 -> 加法 wiring
        vec![(C, 0), (A, 9)],
        vec![(C, 1), (B, 9)],
        vec![(C, 9), (A, 10)],
        vec![(C, 2), (B, 10)],
        vec![(C, 10), (A, 11)],
        vec![(C, 3), (A, 12)],
        vec![(C, 4), (B, 12)],
        vec![(C, 12), (A, 13)],
        vec![(C, 5), (B, 13)],
        vec![(C, 13), (A, 14)],
        vec![(C, 6), (A, 15)],
        vec![(C, 7), (B, 15)],
        vec![(C, 15), (A, 16)],
        vec![(C, 8), (B, 16)],
        vec![(C, 16), (A, 17)],
        // x_j 重複使用
        vec![(B, 0), (B, 3), (B, 6)],
        vec![(B, 1), (B, 4), (B, 7)],
        vec![(B, 2), (B, 5), (B, 8)],
        // z 的 gate 輸出 -> z 的專用 copy 欄位（D column，row 0/1/2）
        vec![(C, 11), (D, 0)],
        vec![(C, 14), (D, 1)],
        vec![(C, 17), (D, 2)],
    ]
}

fn flat_index(col: Column, row: usize) -> usize {
    column_index(col) * NUM_GATES + row
}

/// 把群組展開成一個長度 NUM_COLUMNS*NUM_GATES 的置換陣列（flat index = col*N + row）。
/// 群組內的 cell 依序做 cyclic rotation；沒被任何群組提到的 cell 是 fixed point（sigma[i] = i）。
pub fn build_sigma() -> [usize; NUM_COLUMNS * NUM_GATES] {
    let mut sigma: [usize; NUM_COLUMNS * NUM_GATES] = std::array::from_fn(|i| i);

    for group in build_copy_constraint_groups() {
        let n = group.len();
        for i in 0..n {
            let (col_i, row_i) = group[i];
            let (col_next, row_next) = group[(i + 1) % n];
            sigma[flat_index(col_i, row_i)] = flat_index(col_next, row_next);
        }
    }

    sigma
}

/// 四條欄位各自的 id(col, row) evaluation 向量：id_col[row] = coset_shift(col) * ω^row
pub fn build_id_evaluations(domain: &Radix2EvaluationDomain<Fr>) -> [Vec<Fr>; NUM_COLUMNS] {
    let elements: Vec<Fr> = domain.elements().collect();
    let id_for = |col: Column| -> Vec<Fr> {
        elements.iter().map(|e| coset_shift(col) * e).collect()
    };
    [
        id_for(Column::A),
        id_for(Column::B),
        id_for(Column::C),
        id_for(Column::D),
    ]
}

/// 四條欄位各自的 sigma(col, row) evaluation 向量：sigma_col[row] = id 值 of sigma[flat(col,row)]
pub fn build_sigma_evaluations(domain: &Radix2EvaluationDomain<Fr>) -> [Vec<Fr>; NUM_COLUMNS] {
    let sigma = build_sigma();
    let elements: Vec<Fr> = domain.elements().collect();

    let value_at_flat = |flat: usize| -> Fr {
        let col = column_from_index(flat / NUM_GATES);
        let row = flat % NUM_GATES;
        coset_shift(col) * elements[row]
    };

    let sigma_for = |col: Column| -> Vec<Fr> {
        (0..NUM_GATES)
            .map(|row| value_at_flat(sigma[flat_index(col, row)]))
            .collect()
    };

    [
        sigma_for(Column::A),
        sigma_for(Column::B),
        sigma_for(Column::C),
        sigma_for(Column::D),
    ]
}

/// Grand-product accumulator：
/// Z(ω^0) = 1
/// Z(ω^{i+1}) = Z(ω^i) * ∏_col (col(ω^i) + β·id_col(ω^i) + γ) / ∏_col (col(ω^i) + β·σ_col(ω^i) + γ)
///
/// 回傳長度 NUM_GATES 的 evaluation 向量，剛好對應 domain 的 NUM_GATES 個點。
pub fn compute_accumulator_values(
    wires: &[&[Fr]; NUM_COLUMNS],
    id: &[Vec<Fr>; NUM_COLUMNS],
    sigma: &[Vec<Fr>; NUM_COLUMNS],
    beta: Fr,
    gamma: Fr,
) -> Vec<Fr> {
    let n = wires[0].len();
    let mut z = Vec::with_capacity(n);
    z.push(Fr::one());

    for i in 0..n - 1 {
        let mut num = Fr::one();
        let mut den = Fr::one();
        for col in 0..NUM_COLUMNS {
            num *= wires[col][i] + beta * id[col][i] + gamma;
            den *= wires[col][i] + beta * sigma[col][i] + gamma;
        }
        let prev = *z.last().unwrap();
        z.push(prev * num / den);
    }

    z
}

fn to_poly(domain: &Radix2EvaluationDomain<Fr>, evals: &[Fr]) -> DensePolynomial<Fr> {
    DensePolynomial::from_coefficients_vec(domain.ifft(evals))
}

fn lagrange_poly(domain: &Radix2EvaluationDomain<Fr>, n: usize) -> DensePolynomial<Fr> {
    let mut evals = vec![Fr::zero(); domain.size()];
    evals[n] = Fr::one();
    DensePolynomial::from_coefficients_vec(domain.ifft(&evals))
}

/// C_boundary = L_0(X) * (Z(X) - 1)：強制 Z(ω^0) = 1
fn compute_boundary_check_poly(
    domain: &Radix2EvaluationDomain<Fr>,
    z: &DensePolynomial<Fr>,
) -> DensePolynomial<Fr> {
    let domain_4n: Radix2EvaluationDomain<Fr> = EvaluationDomain::new(4 * domain.size()).unwrap();
    let l0 = lagrange_poly(domain, 0);

    let l0_e = domain_4n.fft(&l0.coeffs);
    let z_e = domain_4n.fft(&z.coeffs);

    let evals: Vec<Fr> = (0..domain_4n.size())
        .map(|i| l0_e[i] * (z_e[i] - Fr::one()))
        .collect();

    DensePolynomial::from_coefficients_vec(domain_4n.ifft(&evals))
}

/// C_transition = Z(X)*∏(col(X)+β·id_col(X)+γ) - Z(ωX)*∏(col(X)+β·σ_col(X)+γ)
///
/// 這是整個 permutation argument 的核心恆等式：如果對所有 domain 上的點都成立，
/// 代表 grand product 沿著 σ 走一圈之後真的繞回原點，也就是每個 copy constraint 群組
/// 內的 wire 值全部相等。
///
/// degree 估算：num(X)/den(X) 各是 4 個 degree ≤ N-1 多項式的乘積（4 欄），
/// 乘上 Z(X) 或 Z(ωX)（degree ≤ N-1）之後，總 degree 最高到 5*(N-1)。
/// N=64 時是 315，所以逐點相乘要在夠大的 domain 上做，這裡用 8N（512 點，
/// 而不是 gate identity 用的 4N）才不會 alias；Z(ωX) 對應 shift-by-8
/// （domain_8n 的步進 8 對應原 domain 的步進 1）。
fn compute_transition_check_poly(
    domain: &Radix2EvaluationDomain<Fr>,
    wires: &[&DensePolynomial<Fr>; NUM_COLUMNS],
    z: &DensePolynomial<Fr>,
    id: &[DensePolynomial<Fr>; NUM_COLUMNS],
    sigma: &[DensePolynomial<Fr>; NUM_COLUMNS],
    beta: Fr,
    gamma: Fr,
) -> DensePolynomial<Fr> {
    const SHIFT: usize = 8;
    let domain_big: Radix2EvaluationDomain<Fr> =
        EvaluationDomain::new(SHIFT * domain.size()).unwrap();
    let ev = |p: &DensePolynomial<Fr>| domain_big.fft(&p.coeffs);

    let wire_e: Vec<Vec<Fr>> = wires.iter().map(|p| ev(p)).collect();
    let id_e: Vec<Vec<Fr>> = id.iter().map(ev).collect();
    let sigma_e: Vec<Vec<Fr>> = sigma.iter().map(ev).collect();

    let mut z_e = ev(z);
    for i in 0..SHIFT {
        z_e.push(z_e[i]);
    }

    let evals: Vec<Fr> = (0..domain_big.size())
        .map(|i| {
            let mut num = Fr::one();
            let mut den = Fr::one();
            for col in 0..NUM_COLUMNS {
                num *= wire_e[col][i] + beta * id_e[col][i] + gamma;
                den *= wire_e[col][i] + beta * sigma_e[col][i] + gamma;
            }
            let z_next = z_e[i + SHIFT]; // Z(ωX)
            z_e[i] * num - z_next * den
        })
        .collect();

    DensePolynomial::from_coefficients_vec(domain_big.ifft(&evals))
}

/// 主函數：把 boundary check 跟 transition check 用挑戰 alpha 合成一個多項式，
/// 除以 vanishing polynomial，回傳 (quotient, remainder)。
/// 若所有 copy constraint 都成立，remainder 必須為零多項式。
pub fn compute(
    domain: &Radix2EvaluationDomain<Fr>,
    a_evals: &[Fr],
    b_evals: &[Fr],
    c_evals: &[Fr],
    d_evals: &[Fr],
    beta: Fr,
    gamma: Fr,
    alpha: Fr,
) -> (DensePolynomial<Fr>, DensePolynomial<Fr>) {
    let id_e = build_id_evaluations(domain);
    let sigma_e = build_sigma_evaluations(domain);

    let wires_evals: [&[Fr]; NUM_COLUMNS] = [a_evals, b_evals, c_evals, d_evals];
    let z_evals = compute_accumulator_values(&wires_evals, &id_e, &sigma_e, beta, gamma);

    let a = to_poly(domain, a_evals);
    let b = to_poly(domain, b_evals);
    let c = to_poly(domain, c_evals);
    let d = to_poly(domain, d_evals);
    let z = to_poly(domain, &z_evals);
    let id: [DensePolynomial<Fr>; NUM_COLUMNS] = std::array::from_fn(|i| to_poly(domain, &id_e[i]));
    let sigma: [DensePolynomial<Fr>; NUM_COLUMNS] =
        std::array::from_fn(|i| to_poly(domain, &sigma_e[i]));

    let wires: [&DensePolynomial<Fr>; NUM_COLUMNS] = [&a, &b, &c, &d];

    let boundary = compute_boundary_check_poly(domain, &z);
    let transition = compute_transition_check_poly(domain, &wires, &z, &id, &sigma, beta, gamma);

    // 用 alpha 把兩條恆等式合成一條：transition + alpha * boundary
    let boundary_scaled =
        DensePolynomial::from_coefficients_vec(boundary.coeffs.iter().map(|v| *v * alpha).collect());
    let combined = &transition + &boundary_scaled;
    combined.divide_by_vanishing_poly(*domain)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gate::layout::{build_gate_layout, build_z_copy_column, compute_wires};

    fn f(v: u64) -> Fr {
        Fr::from(v)
    }

    fn sample_domain() -> Radix2EvaluationDomain<Fr> {
        EvaluationDomain::new(NUM_GATES).unwrap()
    }

    fn sample_witness() -> (Vec<Fr>, Vec<Fr>, Vec<Fr>, Vec<Fr>) {
        let layout = build_gate_layout();
        let w = [
            [f(1), f(2), f(3)],
            [f(4), f(5), f(6)],
            [f(7), f(8), f(9)],
        ];
        let x = [f(1), f(1), f(1)];
        let b = [f(1), f(1), f(1)];
        let (a, b_wire, c) = compute_wires(&layout, &w, &x, &b);
        let d = build_z_copy_column(&c);
        (a, b_wire, c, d)
    }

    #[test]
    fn sigma_is_a_valid_permutation() {
        let sigma = build_sigma();
        let mut seen = vec![false; NUM_COLUMNS * NUM_GATES];
        for &target in sigma.iter() {
            assert!(!seen[target], "sigma 不是一一對應：{target} 被指到兩次");
            seen[target] = true;
        }
    }

    #[test]
    fn accumulator_closes_back_to_one_for_valid_witness() {
        let (a, b, c, d) = sample_witness();
        let domain = sample_domain();
        let id_e = build_id_evaluations(&domain);
        let sigma_e = build_sigma_evaluations(&domain);

        let beta = f(8);
        let gamma = f(13);
        let wires_evals: [&[Fr]; NUM_COLUMNS] = [&a, &b, &c, &d];

        let z = compute_accumulator_values(&wires_evals, &id_e, &sigma_e, beta, gamma);

        let n = a.len();
        let last = n - 1;
        let mut num = Fr::one();
        let mut den = Fr::one();
        for col in 0..NUM_COLUMNS {
            num *= wires_evals[col][last] + beta * id_e[col][last] + gamma;
            den *= wires_evals[col][last] + beta * sigma_e[col][last] + gamma;
        }
        let wrapped = *z.last().unwrap() * num / den;

        assert_eq!(wrapped, Fr::one());
    }

    #[test]
    fn valid_witness_gives_zero_remainder() {
        let (a, b, c, d) = sample_witness();
        let domain = sample_domain();

        let beta = f(8);
        let gamma = f(13);
        let alpha = f(21);

        let (_quotient, remainder) = compute(&domain, &a, &b, &c, &d, beta, gamma, alpha);
        assert!(remainder.is_zero());
    }

    #[test]
    fn wrong_wiring_is_rejected_by_permutation_check() {
        let mut layout = build_gate_layout();
        layout[14].a = crate::gate::layout::WireSource::B(2);

        let w = [
            [f(1), f(2), f(3)],
            [f(4), f(5), f(6)],
            [f(7), f(8), f(9)],
        ];
        let x = [f(1), f(1), f(1)];
        let b = [f(1), f(1), f(1)];
        let (a, b_wire, c) = compute_wires(&layout, &w, &x, &b);
        let d = build_z_copy_column(&c);

        let domain = sample_domain();
        let beta = f(8);
        let gamma = f(13);
        let alpha = f(21);

        let (_quotient, remainder) = compute(&domain, &a, &b_wire, &c, &d, beta, gamma, alpha);
        assert!(!remainder.is_zero());
    }

    #[test]
    fn swapped_x_usage_is_rejected() {
        let mut layout = build_gate_layout();
        layout[3].b = crate::gate::layout::WireSource::X(1);

        let w = [
            [f(1), f(2), f(3)],
            [f(4), f(5), f(6)],
            [f(7), f(8), f(9)],
        ];
        let x = [f(1), f(2), f(3)]; // x0 != x1，才能讓錯誤真的造成差異
        let b = [f(1), f(1), f(1)];
        let (a, b_wire, c) = compute_wires(&layout, &w, &x, &b);
        let d = build_z_copy_column(&c);

        let domain = sample_domain();
        let beta = f(8);
        let gamma = f(13);
        let alpha = f(21);

        let (_quotient, remainder) = compute(&domain, &a, &b_wire, &c, &d, beta, gamma, alpha);
        assert!(!remainder.is_zero());
    }

    #[test]
    fn tampered_z_copy_column_is_rejected() {
        // 直接竄改 D 欄位（不動 A/B/C），模擬「lookup 用的 z 副本」跟真正的 z 不一致。
        let (a, b, c, mut d) = sample_witness();
        d[1] += f(1);

        let domain = sample_domain();
        let beta = f(8);
        let gamma = f(13);
        let alpha = f(21);

        let (_quotient, remainder) = compute(&domain, &a, &b, &c, &d, beta, gamma, alpha);
        assert!(!remainder.is_zero());
    }
}
