use ark_bls12_381::Fr;
use crate::lookup::table::relu::encode_signed;

pub const TABLE_SIZE: usize = 64;
pub const WITNESS_SIZE: usize = TABLE_SIZE - 1;

pub fn encode_and_pad_witness(inputs: &[i64], outputs: &[i64]) -> (Vec<Fr>, Vec<Fr>) {
    let mut encoded_inputs = Vec::new();
    let mut encoded_outputs = Vec::new();

    for input in inputs {
        encoded_inputs.push(encode_signed(*input));
    }

    for output in outputs {
        encoded_outputs.push(encode_signed(*output));
    }

    while encoded_inputs.len() < WITNESS_SIZE {
        encoded_inputs.push(encode_signed(0));
        encoded_outputs.push(encode_signed(0));
    }

    (encoded_inputs, encoded_outputs)
}

pub fn compress_column(inputs: &[Fr], outputs: &[Fr], theta: Fr) -> Vec<Fr> {
    let mut compressed = Vec::new();

    for index in 0..inputs.len() {
        let value = inputs[index] + theta * outputs[index];
        compressed.push(value);
    }

    compressed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_and_pad_integer_witness() {
        let inputs = [2, 4, -1, 0, 15, 6, -24, -2, 15];
        let outputs = [2, 4, 0, 0, 15, 6, 0, 0, 15];

        let (encoded_inputs, encoded_outputs) = encode_and_pad_witness(&inputs, &outputs);

        assert_eq!(encoded_inputs[0], encode_signed(2));
        assert_eq!(encoded_outputs[0], encode_signed(2));

        assert_eq!(encoded_inputs[2], encode_signed(-1));
        assert_eq!(encoded_outputs[2], encode_signed(0));

        assert_eq!(encoded_inputs[9], encode_signed(0));
        assert_eq!(encoded_outputs[9], encode_signed(0));
    }

    #[test]
    fn compress_integer_witness() {
        let inputs = [2, -1];
        let outputs = [2, 0];

        let (encoded_inputs, encoded_outputs) = encode_and_pad_witness(&inputs, &outputs);

        let theta = Fr::from(7u64);
        let compressed = compress_column(&encoded_inputs, &encoded_outputs, theta);

        assert_eq!(compressed.len(), 63);

        let expected_first = encode_signed(2) + theta * encode_signed(2);
        let expected_second = encode_signed(-1) + theta * encode_signed(0);

        assert_eq!(compressed[0], expected_first);
        assert_eq!(compressed[1], expected_second);
    }

    #[test]
    fn compress_relu_table() {
        use crate::lookup::table::relu::ReLUTable;

        let table = ReLUTable::new();

        let inputs = table.encode_input_column();
        let outputs = table.encode_output_column();
        let theta = Fr::from(7u64);

        let compressed = compress_column(&inputs, &outputs, theta);

        let expected_first = encode_signed(-32) + theta * encode_signed(0);
        let expected_last = encode_signed(31) + theta * encode_signed(31);

        assert_eq!(compressed[0], expected_first);
        assert_eq!(compressed[63], expected_last);
    }
}