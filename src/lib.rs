//! Sparse library.

mod array;
mod coo;
mod index;
mod scalar;

// Public API.
pub use coo::{CooMat, CooPat, CoordinateMatrix, CoordinatePattern};
pub use index::Index;
pub use scalar::Scalar;
