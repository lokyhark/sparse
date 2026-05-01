use crate::{Index, Scalar, array::Array};

/// Compressed Sparse Column (CSC) format pattern.
pub type CscPat<I> = CompressedSparseColumnPattern<I>;
/// Compressed Sparse Column (CSC) format matrix.
pub type CscMat<I, S> = CompressedSparseColumnMatrix<I, S>;

/// Compressed Sparse Column (CSC) format pattern.
#[derive(Clone, Debug)]
pub struct CompressedSparseColumnPattern<I: Index> {
    /// Number of rows.
    nrows: I,
    /// Number of columns.
    ncols: I,
    /// Column pointers.
    colptr: Array<I>,
    /// Row indices.
    rowind: Array<I>,
}

impl<I: Index> CompressedSparseColumnPattern<I> {
    /// Returns the number of rows.
    pub fn nrows(&self) -> I {
        self.nrows
    }

    /// Returns the number of columns.
    pub fn ncols(&self) -> I {
        self.ncols
    }

    pub fn colptr(&self) -> &[I] {
        &self.colptr
    }

    pub fn rowind(&self) -> &[I] {
        &self.rowind
    }
}

/// Compressed Sparse Column (CSC) format matrix.
#[derive(Clone, Debug)]
pub struct CompressedSparseColumnMatrix<I: Index, S: Scalar> {
    /// Number of rows.
    nrows: I,
    /// Number of columns.
    ncols: I,
    /// Column pointers.
    colptr: Array<I>,
    /// Row indices.
    rowind: Array<I>,
    /// Non-zero values.
    values: Array<S>,
}

impl<I: Index, S: Scalar> CompressedSparseColumnMatrix<I, S> {
    /// Returns the number of rows.
    pub fn nrows(&self) -> I {
        self.nrows
    }

    /// Returns the number of columns.
    pub fn ncols(&self) -> I {
        self.ncols
    }

    pub fn colptr(&self) -> &[I] {
        &self.colptr
    }

    pub fn rowind(&self) -> &[I] {
        &self.rowind
    }

    pub fn values(&self) -> &[S] {
        &self.values
    }
}
