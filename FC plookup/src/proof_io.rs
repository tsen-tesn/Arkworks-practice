use ark_bls12_381::{Bls12_381, Fr};
use ark_poly_commit::kzg10::Commitment;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use plookup::multiset::proof::{Commitments, EqualityProof, Evaluations};
use std::io::{Read, Write};

/// `plookup::multiset::proof::EqualityProof` 沒有實作 CanonicalSerialize，
/// 這裡依欄位順序手動序列化/反序列化，讓 proof 可以存成檔案讓 verifier 獨立讀取。
pub fn write_equality_proof<W: Write>(
    proof: &EqualityProof,
    mut writer: W,
) -> Result<(), SerializationError> {
    proof.aggregate_witness_comm.serialize_compressed(&mut writer)?;
    proof
        .shifted_aggregate_witness_comm
        .serialize_compressed(&mut writer)?;

    let e = &proof.evaluations;
    e.f.serialize_compressed(&mut writer)?;
    e.t.serialize_compressed(&mut writer)?;
    e.t_omega.serialize_compressed(&mut writer)?;
    e.h_1.serialize_compressed(&mut writer)?;
    e.h_1_omega.serialize_compressed(&mut writer)?;
    e.h_2.serialize_compressed(&mut writer)?;
    e.h_2_omega.serialize_compressed(&mut writer)?;
    e.z.serialize_compressed(&mut writer)?;
    e.z_omega.serialize_compressed(&mut writer)?;

    let c = &proof.commitments;
    c.f.serialize_compressed(&mut writer)?;
    c.q.serialize_compressed(&mut writer)?;
    c.h_1.serialize_compressed(&mut writer)?;
    c.h_2.serialize_compressed(&mut writer)?;
    c.z.serialize_compressed(&mut writer)?;

    Ok(())
}

pub fn read_equality_proof<R: Read>(mut reader: R) -> Result<EqualityProof, SerializationError> {
    let aggregate_witness_comm = Commitment::<Bls12_381>::deserialize_compressed(&mut reader)?;
    let shifted_aggregate_witness_comm =
        Commitment::<Bls12_381>::deserialize_compressed(&mut reader)?;

    let f = Fr::deserialize_compressed(&mut reader)?;
    let t = Fr::deserialize_compressed(&mut reader)?;
    let t_omega = Fr::deserialize_compressed(&mut reader)?;
    let h_1 = Fr::deserialize_compressed(&mut reader)?;
    let h_1_omega = Fr::deserialize_compressed(&mut reader)?;
    let h_2 = Fr::deserialize_compressed(&mut reader)?;
    let h_2_omega = Fr::deserialize_compressed(&mut reader)?;
    let z = Fr::deserialize_compressed(&mut reader)?;
    let z_omega = Fr::deserialize_compressed(&mut reader)?;

    let c_f = Commitment::<Bls12_381>::deserialize_compressed(&mut reader)?;
    let c_q = Commitment::<Bls12_381>::deserialize_compressed(&mut reader)?;
    let c_h_1 = Commitment::<Bls12_381>::deserialize_compressed(&mut reader)?;
    let c_h_2 = Commitment::<Bls12_381>::deserialize_compressed(&mut reader)?;
    let c_z = Commitment::<Bls12_381>::deserialize_compressed(&mut reader)?;

    Ok(EqualityProof {
        aggregate_witness_comm,
        shifted_aggregate_witness_comm,
        evaluations: Evaluations {
            f,
            t,
            t_omega,
            h_1,
            h_1_omega,
            h_2,
            h_2_omega,
            z,
            z_omega,
        },
        commitments: Commitments {
            f: c_f,
            q: c_q,
            h_1: c_h_1,
            h_2: c_h_2,
            z: c_z,
        },
    })
}
