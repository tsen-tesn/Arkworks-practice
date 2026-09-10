use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::PrimeField;
use ark_poly_commit::kzg10::Commitment;
use ark_serialize::CanonicalSerialize;
use merlin::Transcript;

pub trait TranscriptProtocol {
    /// Append a `commitment` with the given `label`.
    fn append_commitment(&mut self, label: &'static [u8], comm: &Commitment<Bls12_381>);

    /// Append a `Scalar` with the given `label`.
    fn append_scalar(&mut self, label: &'static [u8], s: &Fr);

    /// Compute a `label`ed challenge variable.
    fn challenge_scalar(&mut self, label: &'static [u8]) -> Fr;
}

impl TranscriptProtocol for Transcript {
    fn append_commitment(&mut self, label: &'static [u8], comm: &Commitment<Bls12_381>) {
        let mut bytes = Vec::new();
        comm.serialize_compressed(&mut bytes).unwrap();
        self.append_message(label, &bytes);
    }

    fn append_scalar(&mut self, label: &'static [u8], s: &Fr) {
        let mut bytes = Vec::new();
        s.serialize_compressed(&mut bytes).unwrap();
        self.append_message(label, &bytes);
    }

    fn challenge_scalar(&mut self, label: &'static [u8]) -> Fr {
        let mut buf = [0u8; 64];
        self.challenge_bytes(label, &mut buf);
        Fr::from_le_bytes_mod_order(&buf)
    }
}
