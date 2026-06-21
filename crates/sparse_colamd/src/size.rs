use super::int::ColamdInt;

/// Validated dimensions of a sparse matrix for COLAMD algorithm.
///
/// Returned by [`Colamd::check`] after successfully validating the input arguments.
/// Guarantees that `nrows`, `ncols`, and `nnz` are each in the range `(0, isize::MAX]`.
#[derive(Debug)]
pub struct ColamdSize<I: ColamdInt> {
    /// The number of rows in the matrix.
    pub(super) nrows: I,
    /// The number of columns in the matrix.
    pub(super) ncols: I,
    /// The number of non-zero entries in the matrix.
    pub(super) nnz: I,
}

impl<I: ColamdInt> ColamdSize<I> {
    /// Creates a new `ColamdSize` with the given number of rows, columns and non-zero entries.
    pub fn new(nrows: I, ncols: I, nnz: I) -> Self {
        Self { nrows, ncols, nnz }
    }
}
