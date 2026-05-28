use super::int::ColamdInt;

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
