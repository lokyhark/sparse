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
pub trait Index: Sized + Copy + Display + Eq + Ord + Send + Sync + 'static + Seal {
    fn zero() -> Self;
}
impl Index for u8 {
    fn zero() -> Self {
        0
    }
}
impl Index for u16 {
    fn zero() -> Self {
        0
    }
}
impl Index for u32 {
    fn zero() -> Self {
        0
    }
}
impl Index for u64 {
    fn zero() -> Self {
        0
    }
}
impl Index for u128 {
    fn zero() -> Self {
        0
    }
}
impl Index for usize {
    fn zero() -> Self {
        0
    }
}
impl Index for i8 {
    fn zero() -> Self {
        0
    }
}
impl Index for i16 {
    fn zero() -> Self {
        0
    }
}
impl Index for i32 {
    fn zero() -> Self {
        0
    }
}
impl Index for i64 {
    fn zero() -> Self {
        0
    }
}
impl Index for i128 {
    fn zero() -> Self {
        0
    }
}
impl Index for isize {
    fn zero() -> Self {
        0
    }
}
