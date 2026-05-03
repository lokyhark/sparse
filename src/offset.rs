use crate::int::Int;

#[derive(Copy, Clone, Debug)]
pub struct Offset<I: Int>(I);

impl<I: Int> Offset<I> {
    pub fn checked(offset: I) -> Option<Self> {
        let value = offset.try_into().unwrap_or(isize::MIN);
        if value > isize::MIN && value <= I::IDXMAX as isize { Some(Self(offset)) } else { None }
    }

    pub unsafe fn unchecked(&self, offset: I) -> Self {
        Self(offset)
    }

    pub fn get(&self) -> I {
        self.0
    }

    pub unsafe fn index(&self) -> usize {
        // SAFETY: The caller must ensure that the offset is valid.
        unsafe { self.0.index() }
    }
}
