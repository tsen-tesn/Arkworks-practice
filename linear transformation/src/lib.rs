use ark_bls12_381::Fr;
use ark_relations::lc;
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystemRef, LinearCombination, SynthesisError, Variable,
};

pub type ConstraintF = Fr;

// 證明存在 w、x、b，使公開的 y 滿足線性轉換。
// 定義 circuit 資料
#[derive(Clone, Debug, Default)]
pub struct ScalarLinearCircuit {
    // Private witness: weight
    pub w: Option<ConstraintF>,

    // Private witness: input
    pub x: Option<ConstraintF>,

    // Private witness: bias
    pub b: Option<ConstraintF>,

    // Public input: output
    pub y: Option<ConstraintF>,
}

impl ConstraintSynthesizer<ConstraintF> for ScalarLinearCircuit {
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        // Private witness
        let w_var = cs.new_witness_variable(|| self.w.ok_or(SynthesisError::AssignmentMissing))?;

        // Private witness
        let x_var = cs.new_witness_variable(|| self.x.ok_or(SynthesisError::AssignmentMissing))?;

        // Private witness
        let b_var = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;

        // Public input
        let y_var = cs.new_input_variable(|| self.y.ok_or(SynthesisError::AssignmentMissing))?;

        // TODO: 由你建立核心 R1CS constraint。
        cs.enforce_constraint(lc!() + w_var, lc!() + x_var, lc!() + y_var - b_var)?;
        // 你可以使用的 circuit variables：
        // - w_var
        // - x_var
        // - b_var
        // - y_var
        //
        // 目標：限制它們符合 y = w * x + b。
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_relations::r1cs::ConstraintSystem;

    #[test]
    fn scalar_circuit_accepts_correct_output() {
        let circuit = ScalarLinearCircuit {
            w: Some(ConstraintF::from(3u64)),
            x: Some(ConstraintF::from(6u64)),
            b: Some(ConstraintF::from(3u64)),

            // TODO：請手動計算正確的 y，再填入這裡。
            y: Some(ConstraintF::from(21u64)),
        };

        // 建立一個測試用 constraint system
        let cs = ConstraintSystem::<ConstraintF>::new_ref();

        // 把 circuit 轉換成 constraints
        circuit
            .generate_constraints(cs.clone())
            .expect("建立 constraints 失敗");

        let is_satisfied = cs.is_satisfied().expect("檢查 constraint system 失敗");

        println!("is_satisfied = {is_satisfied}");
        assert!(is_satisfied);
    }

    #[test]
    fn scalar_circuit_rejects_wrong_output() {
        let circuit = ScalarLinearCircuit {
            w: Some(ConstraintF::from(3u64)),
            x: Some(ConstraintF::from(6u64)),
            b: Some(ConstraintF::from(3u64)),

            // TODO：故意填入一個錯誤的 y。
            y: Some(ConstraintF::from(10u64)),
        };

        let cs = ConstraintSystem::<ConstraintF>::new_ref();

        circuit
            .generate_constraints(cs.clone())
            .expect("建立 constraints 失敗");

        let is_satisfied = cs.is_satisfied().expect("檢查 constraint system 失敗");

        println!("is_satisfied = {is_satisfied}");
        assert!(!is_satisfied);
    }
}

// 定義 circuit 資料
#[derive(Clone, Debug, Default)]
pub struct MatrixLinearCircuit {
    // Private witness: 2 x 2 weight matrix
    pub w: [[Option<ConstraintF>; 2]; 2],

    // Private witness: input vector
    pub x: [Option<ConstraintF>; 2],

    // Private witness: bias vector
    pub b: [Option<ConstraintF>; 2],

    // Public input: output vector
    pub y: [Option<ConstraintF>; 2],
}

impl ConstraintSynthesizer<ConstraintF> for MatrixLinearCircuit {
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        // Private witness
        let mut w_vars = Vec::with_capacity(2);

        for i in 0..2 {
            let mut row_vars = Vec::with_capacity(2);

            for j in 0..2 {
                let value = self.w[i][j];

                let var =
                    cs.new_witness_variable(|| value.ok_or(SynthesisError::AssignmentMissing))?;

                row_vars.push(var);
            }

            w_vars.push(row_vars);
        }

        // Private witness
        let mut x_vars = Vec::with_capacity(2);

        for i in 0..2 {
            let value = self.x[i];

            let var = cs.new_witness_variable(|| value.ok_or(SynthesisError::AssignmentMissing))?;

            x_vars.push(var);
        }

        // Private witness
        let mut b_vars = Vec::with_capacity(2);
        for i in 0..2 {
            let value = self.b[i];

            let var = cs.new_witness_variable(|| value.ok_or(SynthesisError::AssignmentMissing))?;

            b_vars.push(var);
        }

        // Public input
        let mut y_vars = Vec::with_capacity(2);

        for i in 0..2 {
            let value = self.y[i];

            let var = cs.new_input_variable(|| value.ok_or(SynthesisError::AssignmentMissing))?;

            y_vars.push(var);
        }

        // 配置 intermediate witness
        let mut product_vars = Vec::with_capacity(2);
        for i in 0..2 {
            let mut row_vars = Vec::with_capacity(2);

            for j in 0..2 {
                let w_value = self.w[i][j];
                let x_value = self.x[j];

                let product_var = cs.new_witness_variable(|| {
                    let w = w_value.ok_or(SynthesisError::AssignmentMissing)?;
                    let x = x_value.ok_or(SynthesisError::AssignmentMissing)?;

                    Ok(w * x)
                })?;

                row_vars.push(product_var);
            }
            product_vars.push(row_vars);
        }

        // 乘法 constraints
        // W[i][j] x X[j] = product[i][j]
        for i in 0..2 {
            for j in 0..2 {
                let left_lc = LinearCombination::<ConstraintF>::from(w_vars[i][j]);
                let right_lc = LinearCombination::<ConstraintF>::from(x_vars[j]);
                let output_lc = LinearCombination::<ConstraintF>::from(product_vars[i][j]);

                cs.enforce_constraint(left_lc, right_lc, output_lc)?;
            }
        }

        // 加法 constraints
        // product[i][0] + product[i][1] + B[i] = Y[i]
        for i in 0..2 {
            let output_lc = LinearCombination::<ConstraintF>::from(product_vars[i][0])
                + product_vars[i][1]
                + b_vars[i];

            let one_lc = LinearCombination::<ConstraintF>::from(Variable::One);
            let y_lc = LinearCombination::<ConstraintF>::from(y_vars[i]);

            cs.enforce_constraint(output_lc, one_lc, y_lc)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod matrix_tests {
    use super::*;
    use ark_relations::r1cs::ConstraintSystem;

    fn f(value: u64) -> ConstraintF {
        ConstraintF::from(value)
    }

    #[test]
    fn matrix_circuit_accepts_correct_output() {
        let circuit = MatrixLinearCircuit {
            w: [[Some(f(2)), Some(f(1))], [Some(f(3)), Some(f(4))]],
            x: [Some(f(5)), Some(f(2))],
            b: [Some(f(3)), Some(f(1))],
            y: [Some(f(15)), Some(f(24))],
        };

        let cs = ConstraintSystem::<ConstraintF>::new_ref();

        circuit
            .generate_constraints(cs.clone())
            .expect("建立矩陣 constraints 失敗");

        let is_satisfied = cs.is_satisfied().expect("檢查 constraint system 失敗");

        assert!(is_satisfied);
    }

    #[test]
    fn matrix_circuit_rejects_wrong_output() {
        let circuit = MatrixLinearCircuit {
            w: [[Some(f(2)), Some(f(1))], [Some(f(3)), Some(f(4))]],
            x: [Some(f(5)), Some(f(2))],
            b: [Some(f(3)), Some(f(1))],

            // 第二個輸出故意由 24 改成 25
            y: [Some(f(15)), Some(f(25))],
        };

        let cs = ConstraintSystem::<ConstraintF>::new_ref();

        circuit
            .generate_constraints(cs.clone())
            .expect("建立矩陣 constraints 失敗");

        let is_satisfied = cs.is_satisfied().expect("檢查 constraint system 失敗");

        assert!(!is_satisfied);
    }
}
