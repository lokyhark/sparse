use crate::{Scalar, array::Array, dimension::Dimension, int::Int};

/// Compressed Sparse Column (CSC) format pattern.
pub type CscPat<I> = CompressedSparseColumnPattern<I>;
/// Compressed Sparse Column (CSC) format matrix.
pub type CscMat<I, S> = CompressedSparseColumnMatrix<I, S>;

/// Compressed Sparse Column (CSC) format pattern.
#[derive(Clone, Debug)]
pub struct CompressedSparseColumnPattern<I: Int> {
    /// Number of rows.
    nrows: Dimension<I>,
    /// Number of columns.
    ncols: Dimension<I>,
    /// Column pointers.
    colptr: Array<I>,
    /// Row indices.
    rowind: Array<I>,
}

impl<I: Int> CompressedSparseColumnPattern<I> {
    /// Returns the number of rows.
    pub fn nrows(&self) -> I {
        self.nrows.get()
    }

    /// Returns the number of columns.
    pub fn ncols(&self) -> I {
        self.ncols.get()
    }

    /// Returns the number of non-zero entries.
    pub fn nnz(&self) -> I {
        self.colptr[self.ncols.index()]
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
pub struct CompressedSparseColumnMatrix<I: Int, S: Scalar> {
    /// Number of rows.
    nrows: Dimension<I>,
    /// Number of columns.
    ncols: Dimension<I>,
    /// Column pointers.
    colptr: Array<I>,
    /// Row indices.
    rowind: Array<I>,
    /// Non-zero values.
    values: Array<S>,
}

impl<I: Int, S: Scalar> CompressedSparseColumnMatrix<I, S> {
    /// Returns the number of rows.
    pub fn nrows(&self) -> I {
        self.nrows.get()
    }

    /// Returns the number of columns.
    pub fn ncols(&self) -> I {
        self.ncols.get()
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
