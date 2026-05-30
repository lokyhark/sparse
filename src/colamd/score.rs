use super::int::ColamdInt;

#[derive(Debug)]
pub struct ColamdScore<I: ColamdInt> {
    /// Number of columns alive.
    pub(super) cols: I,
    /// Number of rows alive.
    pub(super) rows: I,
    /// Maximum row degree.
    pub(super) max_degree: I,
    /// Minimum column score.
    pub(super) min_score: I,
}
