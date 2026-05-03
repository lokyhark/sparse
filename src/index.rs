use crate::int::Int;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct Index<I: Int>(I);

impl<I: Int> Index<I> {
    pub fn checked(index: I) -> Option<Self> {
        let value = index.try_into().unwrap_or(usize::MAX);
        if value <= I::IDXMAX { Some(Self(index)) } else { None }
    }

    pub unsafe fn unchecked(&self, index: I) -> Self {
        Self(index)
    }

    pub fn get(&self) -> I {
        self.0
    }

    /// Return an index.
    pub fn index(&self) -> usize {
        // SAFETY: The inner value is guaranteed to be a valid index.
        unsafe { self.0.index() }
    }

    /// Return an offset.
    pub fn offset(&self) -> isize {
        // SAFETY: The inner value is guaranteed to be a valid offset.
        unsafe { self.0.offset() }
    }
}
