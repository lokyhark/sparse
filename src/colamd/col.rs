use super::int::ColamdInt;
/// Represents a column in the COLAMD algorithm.
#[derive(Clone, Debug)]
pub struct ColamdCol<I: ColamdInt> {
    /// The starting index of the column in the compressed sparse column format.
    pub(super) start: I,
    /// The number of non-zero entries in the column.
    pub(super) length: I,
    /// The weight of the column (i.e its thickness).
    pub(super) weight: I,
    /// The rank of the column which is either the score of the column or the actual order of the column in the final permutation.
    pub(super) rank: I,
    /// The doubly linked list prev column pointer.
    pub(super) prev: I,
    /// The doubly linked list next column pointer.
    pub(super) next: I,
}

impl<I: ColamdInt> ColamdCol<I> {
    pub fn alive(&self) -> bool {
        self.start >= I::ZERO
    }

    pub fn dead(&self) -> bool {
        self.start < I::ZERO
    }

    pub fn kill(&mut self) {
        self.start = -I::ONE;
    }
}
