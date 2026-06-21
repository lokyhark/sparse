//! Sparse library.
mod coo;
mod csc;
mod int;
mod scalar;

// Public API.
pub use coo::{CooMat, CooPat, CoordinateMatrix, CoordinatePattern};
pub use csc::{CompressedSparseColumnMatrix, CompressedSparseColumnPattern, CscMat, CscPat};
pub use scalar::Scalar;

// Re-export API.
#[doc(inline)]
pub use sparse_colamd as colamd;
