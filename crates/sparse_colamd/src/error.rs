use std::{error::Error, fmt::Display};

use sparse_array::ArrayError;

use super::int::ColamdInt;

/// COLAMD ordering error.
///
/// # See Also
///
/// [`colamd`](crate::colamd)
#[derive(Debug)]
pub struct ColamdError<I: ColamdInt> {
    /// Error kind.
    kind: ColamdErrorKind<I>,
}

/// COLAMD ordering error kind.
#[derive(Debug)]
pub(super) enum ColamdErrorKind<I: ColamdInt> {
    /// Number of rows must be strictly positive.
    InvalidNumberOfRows { nrows: I },
    /// Number of columns must be strictly positive.
    InvalidNumberOfColumns { ncols: I },
    /// Empty matrix.
    EmptyMatrix,
    /// Elbow room too short.
    InvalidElbowRoom { size: usize, ncols: usize },
    /// Column pointers slice must be non empty.
    EmptyColptrSlice,
    /// Column pointers slice must have a length of `ncols + 1`.
    InvalidColptrSliceLength { actual: usize, expected: usize },
    /// First element of column pointer must be `0`.
    InvalidFirstColPtr { ptr: I },
    /// Column pointers must be in `[0;isize::MAX]`.
    InvalidColPtr { index: usize, ptr: I },
    /// Each column must have a length `>= 0`.
    DecreasingColptrSlice { col: usize, start: I, stop: I },
    /// Row indices slice must have a length of `colptr[ncols]`.
    InvalidRowindSliceLength { actual: usize, expected: usize },
    /// Row indices must be in `0..nrows`.
    OutOfBoundRowIndex { index: usize, row: I, nrows: I },
    /// Row indices must be ordered by column.
    UnorderedRowIndSlice { col: I, ptr: usize },
    /// Array error.
    ArrayError { inner: ArrayError },
}

impl<I: ColamdInt> Display for ColamdError<I> {
    fn fmt(&self, fmt: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.kind {
            ColamdErrorKind::InvalidNumberOfRows { nrows } => write!(fmt, "invalid number of rows; must be > 0 but got {nrows}"),
            ColamdErrorKind::InvalidNumberOfColumns { ncols } => write!(fmt, "invalid number of columns; must be > 0 but got {ncols}"),
            ColamdErrorKind::EmptyMatrix => write!(fmt, "matrix is empty"),
            ColamdErrorKind::InvalidElbowRoom { size, ncols } => write!(fmt, "elbow room too short; must be > {ncols} but got {size}"),
            ColamdErrorKind::EmptyColptrSlice => write!(fmt, "column pointers slice is empty"),
            ColamdErrorKind::InvalidColptrSliceLength { actual, expected } => write!(fmt, "invalid column pointers slice length: expected `{expected}` (`= ncols + 1`) and got `{actual}`"),
            ColamdErrorKind::InvalidFirstColPtr { ptr } => write!(fmt, "invalid column pointers slice first entry: expected `0` and got `{ptr}`"),
            ColamdErrorKind::InvalidColPtr { index, ptr } => write!(fmt, "invalid column pointer at `{index}`; must be >= {} and < {} but got `{ptr}`", super::MIN, super::MAX),
            ColamdErrorKind::DecreasingColptrSlice { col, start, stop } => write!(fmt, "decreasing column pointers at column `{col}` (start = {start}, stop = {stop}"),
            ColamdErrorKind::InvalidRowindSliceLength { actual, expected } => write!(fmt, "invalid row indices slice length: expected `{expected}` (= colptr[ncols]) and got `{actual}`"),
            ColamdErrorKind::OutOfBoundRowIndex { index, row, nrows } => write!(fmt, "out of bounds row index `{row}` at rowind[`{index}`]; must be in `0..{nrows}(=nrows)`"),
            ColamdErrorKind::UnorderedRowIndSlice { col, ptr } => write!(fmt, "row indices in column `{col}` are unordered at index `{ptr}`"),
            ColamdErrorKind::ArrayError { inner } => write!(fmt, "{inner}"),
        }
    }
}

impl<I: ColamdInt> Error for ColamdError<I> {}

impl<I: ColamdInt> From<ColamdErrorKind<I>> for ColamdError<I> {
    fn from(kind: ColamdErrorKind<I>) -> Self {
        Self { kind }
    }
}

impl<I: ColamdInt> From<ArrayError> for ColamdError<I> {
    fn from(value: ArrayError) -> Self {
        Self {
            kind: ColamdErrorKind::ArrayError { inner: value },
        }
    }
}
