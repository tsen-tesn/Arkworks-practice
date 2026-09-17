use ark_bls12_381::Fr;

pub struct ReLURow {
    pub input: i64,
    pub output: i64
}

/// 將 i64 轉成 Fr
pub fn encode_signed(value: i64) -> Fr {
    if value >= 0 {
        Fr::from(value as u64)
    } else {
        -Fr::from(value.unsigned_abs())
    }
}

pub struct ReLUTable {
    rows: Vec<ReLURow>
}

impl ReLUTable {
    /// 建立 relu table
    pub fn new() -> Self {
        let mut rows = Vec::new();
        for input in -32..=31 {
            let output;

            if input < 0 {
                output = 0;
            } else {
                output = input;
            }

            let row = ReLURow {
                input: input,
                output: output,
            };

            rows.push(row);
        }

        ReLUTable {
            rows: rows
        }
    }

    pub fn read(&self, input: i64) -> Option<i64> {
        for row in &self.rows {
            if row.input == input {
                return Some(row.output);
            }
        }

        None
    }

    /// 將表中原本的 i64 編碼的輸入輸出分別轉成 Fr 並存入 encode_input_column, encode_output_column
    pub fn encode_input_column(&self) -> Vec<Fr> {
        let mut column = Vec::new();
        for row in &self.rows {
            let encoded_input = encode_signed(row.input);
            column.push(encoded_input);
        }

        column
    }
    
    pub fn encode_output_column(&self) -> Vec<Fr> {
        let mut column = Vec::new();
        for row in &self.rows {
            let encoded_output = encode_signed(row.output);
            column.push(encoded_output);
        }

        column
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn negative_input() {
        let table = ReLUTable::new();

        assert_eq!(table.rows[8].input, -24);
        assert_eq!(table.rows[8].output, 0);
    }

    #[test]
    fn positive_input() {
        let table = ReLUTable::new();
        assert_eq!(table.rows[38].input, 6);
        assert_eq!(table.rows[38].output, 6);
    }

    #[test]
    fn the_number_of_row() {
        let table = ReLUTable::new();
        assert_eq!(table.rows.len(), 64);
    }

    #[test]
    fn encode_positive_integer() {
        let result = encode_signed(15);
        let expected = Fr::from(15u64);
        assert_eq!(result, expected);
    }

    #[test]
    fn encode_negative_integer() {
        let result = encode_signed(-25);
        let expected = -Fr::from(25u64);
        assert_eq!(result, expected);
    }

    #[test]
    fn encode_input_column() {
        let table = ReLUTable::new();
        let column = table.encode_input_column();

        assert_eq!(column.len(), 64);
        assert_eq!(column[0], encode_signed(-32));
        assert_eq!(column[8], encode_signed(-24));
        assert_eq!(column[32], encode_signed(0));
    }

    #[test]
    fn encode_output_column() {
        let table = ReLUTable::new();
        let column = table.encode_output_column();

        assert_eq!(column.len(), 64);
        assert_eq!(column[18], encode_signed(0));
        assert_eq!(column[32], encode_signed(0));
        assert_eq!(column[47], encode_signed(15));
    }

    #[test]
    fn read_valid_inputs() {
        let table = ReLUTable::new();

        assert_eq!(table.read(-24), Some(0));
        assert_eq!(table.read(24), Some(24));
    }
}