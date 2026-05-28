use super::int::ColamdInt;

/// Configuration parameters for the COLAMD algorithm.
///
/// # Parameters
/// - `dense_row_control`: Controls the threshold for considering a row as dense.
///   A row is considered dense if it has more than `dense_row_control * sqrt(ncols)` non-zero entries
/// - `dense_col_control`: Controls the threshold for considering a column as dense.
///   A column is considered dense if it has more than `dense_col_control * sqrt(min(nrows, ncols))` non-zero entries
///
/// # See Also
/// [`colamd`]
#[derive(Clone, Debug)]
pub struct ColamdConfig<I: ColamdInt> {
    /// The dense *row* control parameter.
    pub(super) dense_row_control: I,
    /// The dense *column* control parameter.
    pub(super) dense_col_control: I,
}

impl<I: ColamdInt> ColamdConfig<I> {
    /// Create a new `ColamdConfig` with default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the dense row control parameter.
    ///
    /// # See Also
    /// - [`ColamdConfig`]
    /// - [`with_dense_col_control`]
    pub fn with_dense_row_control(mut self, value: I) -> Self {
        match value.try_into() {
            Result::<isize, _>::Err(_) => panic!("invalid dense row control value: must bein (0;{}]", isize::MAX),
            Ok(value) if value <= 0 => panic!("dense row control must be positive"),
            _ => (),
        }
        self.dense_row_control = value;
        self
    }

    /// Set the dense column control parameter.
    ///
    /// # See Also
    /// - [`ColamdConfig`]
    /// - [`with_dense_row_control`]
    pub fn with_dense_col_control(mut self, value: I) -> Self {
        match value.try_into() {
            Result::<isize, _>::Err(_) => panic!("invalid dense col control value: must bein (0;{}]", isize::MAX),
            Ok(value) if value <= 0 => panic!("dense col control must be positive"),
            _ => (),
        }
        self.dense_col_control = value;
        self
    }
}

impl<I: ColamdInt> Default for ColamdConfig<I> {
    fn default() -> Self {
        Self {
            dense_row_control: I::TEN,
            dense_col_control: I::TEN,
        }
    }
}
