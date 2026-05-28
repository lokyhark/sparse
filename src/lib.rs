//! Sparse library.

mod array;
mod colamd;
mod coo;
mod csc;
mod int;
mod scalar;

// Public API.
pub use colamd::{Colamd, colamd};
pub use coo::{CooMat, CooPat, CoordinateMatrix, CoordinatePattern};
pub use csc::{CompressedSparseColumnMatrix, CompressedSparseColumnPattern, CscMat, CscPat};
pub use scalar::Scalar;
