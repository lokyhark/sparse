use std::{
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Rem, RemAssign, Sub, SubAssign},
};

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
pub trait Int:
    Sized
    + Copy
    + Debug
    + Display
    + Eq
    + Ord
    + Send
    + Sync
    + 'static
    + Seal
    + TryFrom<usize>
    + TryInto<usize>
    + TryInto<isize>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Rem<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + RemAssign
{
    /// The zero value.
    const ZERO: Self;
    /// The one value.
    const ONE: Self;
    /// The ten value.
    const TEN: Self;

    /// Convert to `usize`.
    fn as_usize(self) -> usize;
    /// Return an iterator over the range `[start, stop)`.
    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self>;
    /// Return the integer square root of `self`.
    fn isqrt(self) -> Self;
}

impl Int for u8 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for u16 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for u32 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for u64 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    #[allow(clippy::cast_possible_truncation)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for u128 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for usize {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    fn as_usize(self) -> usize {
        self
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for i8 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    #[allow(clippy::cast_sign_loss)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for i16 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    #[allow(clippy::cast_sign_loss)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for i32 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    #[allow(clippy::cast_sign_loss)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for i64 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    #[allow(clippy::cast_sign_loss)]
    #[allow(clippy::cast_possible_truncation)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for i128 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    #[allow(clippy::cast_sign_loss)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl Int for isize {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TEN: Self = 10;

    #[allow(clippy::cast_sign_loss)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}
