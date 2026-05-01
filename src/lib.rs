//! Sparse library.

mod array;
mod coo;
mod csc;
mod index;
mod scalar;

// Public API.
pub use coo::{CooMat, CooPat, CoordinateMatrix, CoordinatePattern};
pub use csc::{CompressedSparseColumnMatrix, CompressedSparseColumnPattern, CscMat, CscPat};
pub use index::Index;
pub use scalar::Scalar;
