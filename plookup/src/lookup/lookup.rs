use crate::lookup::table::relu::ReLUTable;

pub struct LookUp {
    table: ReLUTable,
    input_wires: Vec<i64>,
    output_wires: Vec<i64>
}

impl LookUp {
    pub fn new(table: ReLUTable) -> Self {
        LookUp {
            table: table,
            input_wires: Vec::new(),
            output_wires: Vec::new()
        }
    }

    pub fn read(&mut self, input: i64) -> bool {
        let result = self.table.read(input);
        if result.is_none() {
            return false;
        }

        let output = result.unwrap();
        self.input_wires.push(input);
        self.output_wires.push(output);

        true
    }

    /// 測試用
    pub fn input_wires(&self) -> &[i64] {
        &self.input_wires
    }

    pub fn output_wires(&self) -> &[i64] {
        &self.output_wires
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_integer_relu_example() {
        let table = ReLUTable::new();
        let mut lookup = LookUp::new(table);

        let inputs = [2, 4, -1, 0, 15, 6, -24, -2, 15];
        let expected_outputs = [2, 4, 0, 0, 15, 6, 0, 0, 15];

        for input in inputs {
            let is_added = lookup.read(input);
            assert_eq!(is_added, true);
        }

        assert_eq!(lookup.input_wires(), inputs.as_slice());
        assert_eq!(lookup.output_wires(), expected_outputs.as_slice());
    }
}