use std::{
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign},
};

mod private {
    pub trait Seal {}

    impl Seal for i16 {}
    impl Seal for i32 {}
    impl Seal for i64 {}
    impl Seal for isize {}
}

/// COLAMD algorithm supported integer types.
///
/// This trait is sealed to prevent downstream implementations.
pub trait ColamdInt:
    Sized
    + Copy
    + Debug
    + Display
    + Eq
    + Ord
    + Send
    + Sync
    + 'static
    + TryFrom<usize, Error: Debug>
    + TryInto<isize>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Rem<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + RemAssign
    + private::Seal
{
    /// The zero value.
    #[doc(hidden)]
    const ZERO: Self;
    /// The one value.
    #[doc(hidden)]
    const ONE: Self;
    /// The two value.
    #[doc(hidden)]
    const TWO: Self;
    /// The minus one value.
    #[doc(hidden)]
    const NEG: Self;
    /// The ten value.
    #[doc(hidden)]
    const TEN: Self;
    /// The maximum value.
    #[doc(hidden)]
    const MAX: Self;
    /// Convert to `usize`.
    #[doc(hidden)]
    fn as_usize(self) -> usize;
    /// Convert to `f64`.
    #[doc(hidden)]
    fn as_f64(self) -> f64;
    /// Convert from `f64`.
    #[doc(hidden)]
    fn from_f64(value: f64) -> Self;
    /// Return an iterator over the range `[start, stop)`.
    #[doc(hidden)]
    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self>;
    /// Return the integer square root of `self`.
    #[doc(hidden)]
    #[must_use]
    fn isqrt(self) -> Self;
}

impl ColamdInt for i16 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TWO: Self = 2;
    const NEG: Self = -1;
    const TEN: Self = 10;
    const MAX: Self = Self::MAX;

    #[cfg(debug_assertions)]
    fn as_usize(self) -> usize {
        self.try_into().unwrap_or_else(|_| panic!("failed to convert to usize: {}", self))
    }

    #[cfg(not(debug_assertions))]
    #[allow(clippy::cast_sign_loss)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn as_f64(self) -> f64 {
        self as f64
    }

    fn from_f64(value: f64) -> Self {
        value as Self
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl ColamdInt for i32 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TWO: Self = 2;
    const NEG: Self = -1;
    const TEN: Self = 10;
    const MAX: Self = Self::MAX;

    #[cfg(debug_assertions)]
    fn as_usize(self) -> usize {
        self.try_into().unwrap_or_else(|_| panic!("failed to convert to usize: {}", self))
    }

    #[cfg(not(debug_assertions))]
    #[allow(clippy::cast_sign_loss)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn as_f64(self) -> f64 {
        self as f64
    }

    fn from_f64(value: f64) -> Self {
        value as Self
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl ColamdInt for i64 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TWO: Self = 2;
    const NEG: Self = -1;
    const TEN: Self = 10;
    const MAX: Self = Self::MAX;

    #[cfg(debug_assertions)]
    fn as_usize(self) -> usize {
        self.try_into().unwrap_or_else(|_| panic!("failed to convert to usize: {}", self))
    }

    #[cfg(not(debug_assertions))]
    #[allow(clippy::cast_sign_loss)]
    #[allow(clippy::cast_possible_truncation)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn as_f64(self) -> f64 {
        self as f64
    }

    fn from_f64(value: f64) -> Self {
        value as Self
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}

impl ColamdInt for isize {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const TWO: Self = 2;
    const NEG: Self = -1;
    const TEN: Self = 10;
    const MAX: Self = Self::MAX;

    #[cfg(debug_assertions)]
    fn as_usize(self) -> usize {
        self.try_into().unwrap_or_else(|_| panic!("failed to convert to usize: {}", self))
    }

    #[cfg(not(debug_assertions))]
    #[allow(clippy::cast_sign_loss)]
    fn as_usize(self) -> usize {
        self as usize
    }

    fn as_f64(self) -> f64 {
        self as f64
    }

    fn from_f64(value: f64) -> Self {
        value as Self
    }

    fn range(start: Self, stop: Self) -> impl DoubleEndedIterator<Item = Self> {
        start..stop
    }

    fn isqrt(self) -> Self {
        self.isqrt()
    }
}
