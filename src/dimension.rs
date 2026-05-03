use std::usize;

use crate::int::Int;

#[derive(Clone, Debug)]
#[repr(transparent)]
pub struct Dimension<I: Int>(I);

/// A dimension.
///
/// It is a non-zero positive integer that is less than or equal to `isize::MAX`.
impl<I: Int> Dimension<I> {
    pub const MAX: usize = isize::MAX as usize;
    pub fn checked(dimension: I) -> Option<Self> {
        let value = dimension.try_into().unwrap_or(usize::MAX);
        if value > 0 && value <= Self::MAX { Some(Self(dimension)) } else { None }
    }

    pub unsafe fn unchecked(&self, dimension: I) -> Self {
        Self(dimension)
    }

    pub fn get(&self) -> I {
        self.0
    }

    pub fn index(&self) -> usize {
        // SAFETY: The caller must ensure that the dimension is valid and it leads to a valid index.
        unsafe { self.0.index() }
    }
}
