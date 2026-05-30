use crate::{array::Array, colamd::int::ColamdInt};

/// COLAMD ordering result.
///
/// # See Also
///
/// [`colamd`]
#[derive(Clone, Debug)]
pub struct ColamdResult<I: ColamdInt> {
    /// Number of columns alive.
    pub(super) cols: I,
    /// Number of rows alive.
    pub(super) rows: I,
    /// Maximum row degree.
    pub(super) max_degree: I,
    /// Minimum column score.
    pub(super) min_score: I,
    /// Column ordering.
    pub(super) order: Array<I>,
}

impl<I: ColamdInt> ColamdResult<I> {
    pub fn cols(&self) -> I {
        self.cols
    }

    pub fn rows(&self) -> I {
        self.rows
    }

    pub fn max_degree(&self) -> I {
        self.max_degree
    }

    pub fn min_score(&self) -> I {
        self.min_score
    }

    pub fn order(&self) -> &[I] {
        &self.order
    }
}
