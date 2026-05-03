//! Sparse library.

mod array;
mod coo;
mod csc;
mod dimension;
mod index;
mod int;
mod offset;
mod scalar;

// Public API.
pub use coo::{CooMat, CooPat, CoordinateMatrix, CoordinatePattern};
pub use csc::{CompressedSparseColumnMatrix, CompressedSparseColumnPattern, CscMat, CscPat};
pub use scalar::Scalar;
