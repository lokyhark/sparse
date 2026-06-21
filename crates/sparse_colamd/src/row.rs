use super::int::ColamdInt;

/// Represents a row in the COLAMD algorithm.
#[derive(Clone, Debug)]
pub struct ColamdRow<I: ColamdInt> {
    /// The starting index of the row in the compressed sparse row format.
    pub(super) start: I,
    /// The number of non-zero entries in the row.
    pub(super) length: I,
    /// The degree of the row (i.e the number of columns in which the row is present).
    pub(super) degree: I,
    /// The mark used to compute set differences.
    pub(super) mark: I,
}

impl<I: ColamdInt> ColamdRow<I> {
    pub fn alive(&self) -> bool {
        self.mark >= I::ZERO
    }

    pub fn dead(&self) -> bool {
        self.mark < I::ZERO
    }

    pub fn kill(&mut self) {
        self.mark = -I::ONE;
    }
}
