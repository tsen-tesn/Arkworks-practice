use ark_bls12_381::Fr;
use ark_poly::{polynomial::univariate::DensePolynomial as Polynomial, DenseUVPolynomial, EvaluationDomain};

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct MultiSet(pub Vec<Fr>);

impl MultiSet {
    // Creates an empty Multiset
    pub fn new() -> MultiSet {
        MultiSet(vec![])
    }

    pub fn from_slice(slice: &[Fr]) -> MultiSet {
        MultiSet(slice.to_vec())
    }

    // Pushes a value onto the end of the set
    pub fn push(&mut self, value: Fr) {
        self.0.push(value)
    }

    // Pushes 'n' elements into the multiset
    pub fn extend(&mut self, n: usize, value: Fr) {
        let elements = vec![value; n];
        self.0.extend(elements);
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn last(&self) -> Fr {
        *self.0.last().unwrap()
    }

    pub fn as_slice(&self) -> &[Fr] {
        &self.0
    }

    // Checks if an element is in the MultiSet
    pub fn contains(&self, element: &Fr) -> bool {
        self.0.contains(element)
    }

    // Checks whether self is a subset of other
    pub fn is_subset_of(&self, other: &MultiSet) -> bool {
        let mut is_subset = true;

        for x in self.0.iter() {
            is_subset = other.contains(x);
            if is_subset == false {
                break;
            }
        }
        is_subset
    }

    // Returns the position of the element in the Multiset
    // Returns None if the element is not in the Multiset
    pub fn position(&self, element: &Fr) -> usize {
        let index = self.0.iter().position(|&x| x == *element).unwrap();
        index
    }

    // s = sort(f || t)
    pub fn concatenate_and_sort(&self, t: &MultiSet) -> MultiSet {
        assert!(self.is_subset_of(t));
        let mut result = t.clone();

        for element in self.0.iter() {
            let index = result.position(element);
            result.0.insert(index, *element);
        }

        result
    }
    
    // Splits a multiset into halves as specified by the paper
    // If s = [1,2,3,4,5,6,7], we can deduce n using |s| = 2 * n + 1 = 7
    // n is therefore 3
    // We split s into two MultiSets of size n+1 each
    // s_0 = [1,2,3,4] ,|s_0| = n+1 = 4
    // s_1 = [4,5,6,7] , |s_1| = n+1 = 4
    // Notice that the last element of the first half equals the first element in the second half
    // This is specified in the paper
    pub fn halve(&self) -> (MultiSet, MultiSet) {
        let length = self.0.len();
        let first_half = MultiSet::from_slice(&self.0[0..=length / 2]);
        let second_half = MultiSet::from_slice(&self.0[length / 2..]);
        (first_half, second_half)
    }

    // Treats each element in the multiset as evaluation points
    // Computes IFFT of the set of evaluation points
    // and returns the coefficients as a Polynomial data structure
    pub fn to_polynomial<E: EvaluationDomain<Fr>>(&self, domain: &E) -> Polynomial<Fr> {
        Polynomial::from_coefficients_vec(domain.ifft(&self.0))
    }
}

#[cfg(test)]
mod test {
    use ark_poly::{Polynomial, Radix2EvaluationDomain};
    use crate::lookup::table::relu::ReLUTable;
    use crate::lookup::lookup::LookUp;
    use crate::lookup::proof::{compress_column, encode_and_pad_witness};

    use super::*;

    fn fr(values: &[i64]) -> Vec<Fr> {
        values
            .iter()
            .map(|v| {
                if *v >= 0 {
                    Fr::from(*v as u64)
                } else {
                    -Fr::from(v.unsigned_abs())
                }
            })
            .collect()
    }

    // ---- basic operations ----

    #[test]
    fn new_multiset_is_empty() {
        let ms = MultiSet::new();
        assert_eq!(ms.len(), 0);
        assert!(ms.is_empty());
    }

    #[test]
    fn push_appends_value() {
        let mut ms = MultiSet::new();
        ms.push(Fr::from(5u64));

        assert_eq!(ms.len(), 1);
        assert_eq!(ms.last(), Fr::from(5u64));
    }

    #[test]
    fn extend_appends_n_copies() {
        let mut ms = MultiSet::new();
        ms.extend(3, Fr::from(7u64));

        assert_eq!(ms.len(), 3);
        assert_eq!(ms.as_slice(), fr(&[7, 7, 7]).as_slice());
    }

    #[test]
    fn from_slice_matches_input() {
        let values = fr(&[1, 2, 3]);
        let ms = MultiSet::from_slice(&values);

        assert_eq!(ms.as_slice(), values.as_slice());
    }

    // ---- membership ----

    #[test]
    fn contains_detects_membership() {
        let ms = MultiSet::from_slice(&fr(&[1, 2, 3]));

        assert!(ms.contains(&Fr::from(2u64)));
        assert!(!ms.contains(&Fr::from(9u64)));
    }

    #[test]
    fn is_subset_of_true_when_every_element_exists() {
        let table = MultiSet::from_slice(&fr(&[1, 2, 3, 4, 5]));
        let witness = MultiSet::from_slice(&fr(&[2, 2, 5]));

        assert!(witness.is_subset_of(&table));
    }

    #[test]
    fn is_subset_of_false_when_an_element_is_missing() {
        let table = MultiSet::from_slice(&fr(&[1, 2, 3, 4, 5]));
        let witness = MultiSet::from_slice(&fr(&[2, 9]));

        assert!(!witness.is_subset_of(&table));
    }

    #[test]
    fn position_finds_index_of_first_match() {
        let ms = MultiSet::from_slice(&fr(&[1, 2, 3]));

        assert_eq!(ms.position(&Fr::from(2u64)), 1);
    }

    // ---- concatenate_and_sort ----

    #[test]
    fn concatenate_and_sort_preserves_lengths() {
        // table: n + 1 = 8 rows, witness: n = 7 rows
        let table = MultiSet::from_slice(&fr(&[1, 2, 3, 4, 5, 6, 7, 8]));
        let witness = MultiSet::from_slice(&fr(&[2, 2, 5, 8, 1, 1, 3]));

        let sorted = witness.concatenate_and_sort(&table);

        assert_eq!(table.len(), 8);
        assert_eq!(witness.len(), 7);
        assert_eq!(sorted.len(), 2 * 7 + 1);
    }

    #[test]
    fn concatenate_and_sort_preserves_multiplicity_and_order() {
        let table = MultiSet::from_slice(&fr(&[1, 2, 3, 4, 5, 6, 7, 8]));
        let witness = MultiSet::from_slice(&fr(&[2, 2, 5, 8, 1, 1, 3]));

        let sorted = witness.concatenate_and_sort(&table);

        // Hand-verified: each table row v that is queried k times in the
        // witness appears k + 1 times, grouped together, in table order.
        let expected = fr(&[1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 5, 6, 7, 8, 8]);
        assert_eq!(sorted.as_slice(), expected.as_slice());
    }

    // ---- halve ----

    #[test]
    fn halve_splits_with_overlap() {
        let table = MultiSet::from_slice(&fr(&[1, 2, 3, 4, 5, 6, 7, 8]));
        let witness = MultiSet::from_slice(&fr(&[2, 2, 5, 8, 1, 1, 3]));
        let sorted = witness.concatenate_and_sort(&table);

        let (h1, h2) = sorted.halve();

        assert_eq!(h1.len(), 8);
        assert_eq!(h2.len(), 8);
        // h1(omega^n) == h2(1)
        assert_eq!(h1.last(), h2.as_slice()[0]);
    }

    // ---- to_polynomial ----

    #[test]
    fn to_polynomial_round_trips_through_evaluation() {
        let values = fr(&[10, 20, 30, 40]);
        let ms = MultiSet::from_slice(&values);

        let domain = Radix2EvaluationDomain::<Fr>::new(values.len()).unwrap();
        let poly = ms.to_polynomial(&domain);

        for (i, expected) in values.iter().enumerate() {
            let point = domain.element(i);
            assert_eq!(poly.evaluate(&point), *expected);
        }
    }

    // ---- ReLU test ----
    #[test]
    fn relu_witness_concatenate_and_sort() {
        let table = ReLUTable::new();
        let mut lookup = LookUp::new(ReLUTable::new());

        let inputs = [2, 4, -1, 0, 15, 6, -24, -2, 15];
        for input in inputs {
            assert!(lookup.read(input));
        }

        let (encoded_inputs, encoded_outputs) = encode_and_pad_witness(lookup.input_wires(), lookup.output_wires());

        let theta = Fr::from(7u64);
        let compressed_witness = compress_column(&encoded_inputs, &encoded_outputs, theta);
        let compressed_table = compress_column(&table.encode_input_column(), &table.encode_output_column(), theta);
        
        let witness_ms = MultiSet::from_slice(&compressed_witness);
        let table_ms = MultiSet::from_slice(&compressed_table);

        let sorted = witness_ms.concatenate_and_sort(&table_ms);

        assert_eq!(witness_ms.len(), 63);
        assert_eq!(table_ms.len(), 64);
        assert_eq!(sorted.len(), 127);   // 2 * 63 + 1

        let (h1, h2) = sorted.halve();
        assert_eq!(h1.len(), 64);
        assert_eq!(h2.len(), 64);
        assert_eq!(h1.last(), h2.as_slice()[0]);
    }
}