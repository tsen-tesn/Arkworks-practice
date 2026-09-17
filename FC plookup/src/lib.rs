pub mod arithmetic;
pub mod binding;
pub mod columns;
pub mod proof_io;

/// KZG trusted setup 支援的最大多項式次數，跟 plookup crate 既有測試的量級一致。
pub const MAX_DEGREE: usize = 256;
