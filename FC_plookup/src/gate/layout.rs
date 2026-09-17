use ark_bls12_381::Fr;
use ark_std::Zero;

/// 3x3 矩陣
pub const MATRIX_DIM: usize = 3;
/// gate 總數。跟 ReLU table 的大小（plookup::lookup::proof::TABLE_SIZE = 64）一致，
/// 這樣電路的 wire polynomial、z 的專用 copy 欄位（D column，見 gate::permutation）
/// 才能跟 plookup 的 witness polynomial 活在同一個 domain 上，
/// 不需要跨 domain 開值就能把 z、y 綁進 lookup 論證（見 FC_plookup 整體設計討論）。
pub const NUM_GATES: usize = 64;

/// z0, z1, z2 分別落在哪個 gate 的輸出（c 值）上，算術正確性由 M2 gate identity 保證，
/// wiring 正確性由 M3 permutation argument 保證。
pub const Z_GATE_INDICES: [usize; 3] = [11, 14, 17];

/// 一個 gate 是乘法還是加法，決定要驗證哪種恆等式
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GateType {
    Mul, // a * b = c
    Add, // a + b = c
}

/// 描述某個 wire (a 或 b) 的值應該從哪裡來
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WireSource {
    W(usize, usize),   // w_ij
    X(usize),          // x_j
    B(usize),          // b_i
    GateOutput(usize), // 引用某個更早的 gate 的 c 值
    Zero,              // padding 用的常數 0
}

/// 一個 gate 的完整描述：類型 + a、b 兩個輸入的來源
/// c 不需要來源描述，因為 c 是這個 gate 自己算出來的結果，
/// 之後其他 gate 要引用它，是透過 WireSource::GateOutput(這個 gate 的編號)
#[derive(Clone, Copy, Debug)]
pub struct GateLayout {
    pub gate_type: GateType,
    pub a: WireSource,
    pub b: WireSource,
}

fn mul(a: WireSource, b: WireSource) -> GateLayout {
    GateLayout { gate_type: GateType::Mul, a, b }
}

fn add(a: WireSource, b: WireSource) -> GateLayout {
    GateLayout { gate_type: GateType::Add, a, b }
}

/// 產生固定、公開的電路拓樸（跟 w, x, b 的具體數值完全無關）
///
/// 編號規則：
/// 0..=8   乘法 gate（p_ij = w_ij * x_j）
/// 9..=17  加法 gate（逐步把三個 partial product 跟 bias 加起來，算出 z_i）
/// 18..=63 padding gate（0 + 0 = 0，湊滿 64 個 gate 給後面的 FFT 用）
pub fn build_gate_layout() -> [GateLayout; NUM_GATES] {
    let padding = add(WireSource::Zero, WireSource::Zero);
    let mut layout = [padding; NUM_GATES];

    let real_gates: [GateLayout; 18] = [
        // 乘法 gate 0..=8：p_ij = w_ij * x_j
        mul(WireSource::W(0, 0), WireSource::X(0)),
        mul(WireSource::W(0, 1), WireSource::X(1)),
        mul(WireSource::W(0, 2), WireSource::X(2)),
        mul(WireSource::W(1, 0), WireSource::X(0)),
        mul(WireSource::W(1, 1), WireSource::X(1)),
        mul(WireSource::W(1, 2), WireSource::X(2)),
        mul(WireSource::W(2, 0), WireSource::X(0)),
        mul(WireSource::W(2, 1), WireSource::X(1)),
        mul(WireSource::W(2, 2), WireSource::X(2)),
        // 加法 gate 9..=17：逐步加總每列的 partial product 跟 bias
        add(WireSource::GateOutput(0), WireSource::GateOutput(1)),  // 9:  t01
        add(WireSource::GateOutput(9), WireSource::GateOutput(2)),  // 10: t02
        add(WireSource::GateOutput(10), WireSource::B(0)),          // 11: z0
        add(WireSource::GateOutput(3), WireSource::GateOutput(4)),  // 12: t11
        add(WireSource::GateOutput(12), WireSource::GateOutput(5)), // 13: t12
        add(WireSource::GateOutput(13), WireSource::B(1)),          // 14: z1
        add(WireSource::GateOutput(6), WireSource::GateOutput(7)),  // 15: t21
        add(WireSource::GateOutput(15), WireSource::GateOutput(8)), // 16: t22
        add(WireSource::GateOutput(16), WireSource::B(2)),          // 17: z2
    ];

    layout[..18].copy_from_slice(&real_gates);
    layout
}

/// 把某個 WireSource 解析成具體的 Fr 值。
/// `c_so_far` 是目前已經算出來的 c 值（依照 gate 編號順序遞增填入），
/// GateOutput(k) 只能引用 k < 目前正在算的 gate 編號，這個限制由呼叫端（compute_wires）保證。
fn resolve_source(
    source: WireSource,
    w: &[[Fr; MATRIX_DIM]; MATRIX_DIM],
    x: &[Fr; MATRIX_DIM],
    b: &[Fr; MATRIX_DIM],
    c_so_far: &[Fr],
) -> Fr {
    match source {
        WireSource::W(i, j) => w[i][j],
        WireSource::X(j) => x[j],
        WireSource::B(i) => b[i],
        WireSource::GateOutput(k) => c_so_far[k],
        WireSource::Zero => Fr::zero(),
    }
}

/// 給定具體的 w, x, b，依照 layout 算出三條長度 NUM_GATES 的 wire 向量 (a, b, c)
pub fn compute_wires(
    layout: &[GateLayout; NUM_GATES],
    w: &[[Fr; MATRIX_DIM]; MATRIX_DIM],
    x: &[Fr; MATRIX_DIM],
    b: &[Fr; MATRIX_DIM],
) -> (Vec<Fr>, Vec<Fr>, Vec<Fr>) {
    let mut a_vec = Vec::with_capacity(NUM_GATES);
    let mut b_vec = Vec::with_capacity(NUM_GATES);
    let mut c_vec = Vec::with_capacity(NUM_GATES);

    for gate in layout.iter() {
        let a_val = resolve_source(gate.a, w, x, b, &c_vec);
        let b_val = resolve_source(gate.b, w, x, b, &c_vec);
        let c_val = match gate.gate_type {
            GateType::Mul => a_val * b_val,
            GateType::Add => a_val + b_val,
        };

        a_vec.push(a_val);
        b_vec.push(b_val);
        c_vec.push(c_val);
    }

    (a_vec, b_vec, c_vec)
}

/// 從 c 向量取出 z0, z1, z2（固定位置由 Z_GATE_INDICES 決定）
pub fn extract_z_outputs(c: &[Fr]) -> [Fr; MATRIX_DIM] {
    [
        c[Z_GATE_INDICES[0]],
        c[Z_GATE_INDICES[1]],
        c[Z_GATE_INDICES[2]],
    ]
}

/// 建立 y 的私密欄位：長度 NUM_GATES，只有前三格是 y0,y1,y2，其餘補 0。
/// 跟 build_z_copy_column 是同樣的資料形狀，方便跟 D 欄位一起做同態組合（M4）。
pub fn build_y_column(y: &[Fr; MATRIX_DIM]) -> Vec<Fr> {
    let mut col = vec![Fr::zero(); NUM_GATES];
    col[0] = y[0];
    col[1] = y[1];
    col[2] = y[2];
    col
}

/// 建立餵給 lookup 論證用的「z 專用欄位」（D column，見 gate::permutation）：
/// 長度 NUM_GATES，只有前三格是真正的 z0,z1,z2，其餘補 0。
/// 這欄位透過 M3 的 permutation argument 被強制等於 Z_GATE_INDICES 算出來的值，
/// 是整個綁定機制（M4）的關鍵：z 全程不會以明文形式開值給 verifier。
pub fn build_z_copy_column(c: &[Fr]) -> Vec<Fr> {
    let mut d = vec![Fr::zero(); NUM_GATES];
    let z = extract_z_outputs(c);
    d[0] = z[0];
    d[1] = z[1];
    d[2] = z[2];
    d
}

/// 驗證每個 gate 自己的恆等式是否成立（a op b == c）
/// 注意：這只檢查「單一 gate 內部算對了」，不保證 wiring 接對（那是 M3 的責任）
pub fn verify_gate_identities(
    layout: &[GateLayout; NUM_GATES],
    a: &[Fr],
    b: &[Fr],
    c: &[Fr],
) -> bool {
    for k in 0..NUM_GATES {
        let expected = match layout[k].gate_type {
            GateType::Mul => a[k] * b[k],
            GateType::Add => a[k] + b[k],
        };
        if expected != c[k] {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn f(v: u64) -> Fr {
        Fr::from(v)
    }

    fn sample_wxb() -> (
        [[Fr; MATRIX_DIM]; MATRIX_DIM],
        [Fr; MATRIX_DIM],
        [Fr; MATRIX_DIM],
    ) {
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
    fn layout_has_expected_gate_counts() {
        let layout = build_gate_layout();
        assert_eq!(layout.len(), NUM_GATES);

        let mul_count = layout.iter().filter(|g| g.gate_type == GateType::Mul).count();
        let add_count = layout.iter().filter(|g| g.gate_type == GateType::Add).count();

        assert_eq!(mul_count, 9);
        // 9 個真的加法 gate + 46 個 padding
        assert_eq!(add_count, 55);
    }

    #[test]
    fn compute_wires_matches_manual_z() {
        let layout = build_gate_layout();
        let (w, x, b) = sample_wxb();

        // 手算：z_i = w_i0 + w_i1 + w_i2 + b_i（因為 x 全是 1）
        let expected_z = [f(1 + 2 + 3 + 1), f(4 + 5 + 6 + 1), f(7 + 8 + 9 + 1)];

        let (a_vec, b_vec, c_vec) = compute_wires(&layout, &w, &x, &b);

        assert_eq!(extract_z_outputs(&c_vec), expected_z);
        assert!(verify_gate_identities(&layout, &a_vec, &b_vec, &c_vec));
    }

    #[test]
    fn z_copy_column_matches_and_pads_with_zero() {
        let layout = build_gate_layout();
        let (w, x, b) = sample_wxb();
        let (_, _, c_vec) = compute_wires(&layout, &w, &x, &b);

        let d = build_z_copy_column(&c_vec);
        assert_eq!(d.len(), NUM_GATES);
        assert_eq!(d[0..3], extract_z_outputs(&c_vec));
        assert!(d[3..].iter().all(|v| v.is_zero()));
    }

    #[test]
    fn tampering_c_breaks_identity_check() {
        let layout = build_gate_layout();
        let (w, x, b) = sample_wxb();

        let (a_vec, b_vec, mut c_vec) = compute_wires(&layout, &w, &x, &b);
        c_vec[5] += f(1); // 故意改壞其中一個 gate 的輸出

        assert!(!verify_gate_identities(&layout, &a_vec, &b_vec, &c_vec));
    }

    #[test]
    fn wrong_wiring_can_still_pass_gate_identity_but_gives_wrong_z() {
        // gate identity 不保證 wiring 接對：只改接線描述，每個 gate 自己仍然成立，
        // 但 z1 會是錯的——這正是 M3 permutation argument 存在的理由。
        let mut layout = build_gate_layout();
        layout[14].a = WireSource::B(2); // 原本該接 GateOutput(13) + B(1)

        let (w, x, b) = sample_wxb();
        let (a_vec, b_vec, c_vec) = compute_wires(&layout, &w, &x, &b);
        let z = extract_z_outputs(&c_vec);

        assert!(verify_gate_identities(&layout, &a_vec, &b_vec, &c_vec));
        assert_ne!(z[1], f(4 + 5 + 6 + 1));
    }
}
