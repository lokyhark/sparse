use std::fmt::Display;

/// Seal trait for index types.
pub trait Seal {}
impl Seal for u8 {}
impl Seal for u16 {}
impl Seal for u32 {}
impl Seal for u64 {}
impl Seal for u128 {}
impl Seal for usize {}
impl Seal for i8 {}
impl Seal for i16 {}
impl Seal for i32 {}
impl Seal for i64 {}
impl Seal for i128 {}
impl Seal for isize {}

/// A trait for index types that can be used in sparse matrices.
pub trait Int: Sized + Copy + Display + Eq + Ord + Send + Sync + 'static + Seal + TryInto<isize> + TryInto<usize> {
    const IDXMAX: usize;
    fn zero() -> Self;
    unsafe fn index(self) -> usize;
    unsafe fn offset(self) -> isize;
}

impl Int for u8 {
    const IDXMAX: usize = i8::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for u16 {
    const IDXMAX: usize = i16::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for u32 {
    #[cfg(target_pointer_width = "16")]
    const IDXMAX: usize = i16::MAX as usize;
    #[cfg(target_pointer_width = "32")]
    const IDXMAX: usize = i32::MAX as usize;
    #[cfg(target_pointer_width = "64")]
    const IDXMAX: usize = i32::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for u64 {
    #[cfg(target_pointer_width = "16")]
    const IDXMAX: usize = i16::MAX as usize;
    #[cfg(target_pointer_width = "32")]
    const IDXMAX: usize = i32::MAX as usize;
    #[cfg(target_pointer_width = "64")]
    const IDXMAX: usize = i64::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for u128 {
    #[cfg(target_pointer_width = "16")]
    const IDXMAX: usize = i16::MAX as usize;
    #[cfg(target_pointer_width = "32")]
    const IDXMAX: usize = i32::MAX as usize;
    #[cfg(target_pointer_width = "64")]
    const IDXMAX: usize = i64::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for usize {
    const IDXMAX: usize = isize::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for i8 {
    const IDXMAX: usize = i8::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for i16 {
    const IDXMAX: usize = i16::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for i32 {
    #[cfg(target_pointer_width = "16")]
    const IDXMAX: usize = i16::MAX as usize;
    #[cfg(target_pointer_width = "32")]
    const IDXMAX: usize = i32::MAX as usize;
    #[cfg(target_pointer_width = "64")]
    const IDXMAX: usize = i32::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for i64 {
    #[cfg(target_pointer_width = "16")]
    const IDXMAX: usize = i16::MAX as usize;
    #[cfg(target_pointer_width = "32")]
    const IDXMAX: usize = i32::MAX as usize;
    #[cfg(target_pointer_width = "64")]
    const IDXMAX: usize = i64::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for i128 {
    #[cfg(target_pointer_width = "16")]
    const IDXMAX: usize = i16::MAX as usize;
    #[cfg(target_pointer_width = "32")]
    const IDXMAX: usize = i32::MAX as usize;
    #[cfg(target_pointer_width = "64")]
    const IDXMAX: usize = i64::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}

impl Int for isize {
    const IDXMAX: usize = isize::MAX as usize;

    fn zero() -> Self {
        0
    }

    unsafe fn index(self) -> usize {
        self as usize
    }

    unsafe fn offset(self) -> isize {
        self as isize
    }
}
