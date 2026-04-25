//! Scalar types for sparse matrices.

// Following types are supported :
// - bool
// - signed integers: i8, i16, i32, i64, i128, isize
// - unsigned integers: u8, u16, u32, u64, u128, usize
// - floating point numbers: f32, f64
// Following types should be added when stabilized.
// - floating point numbers: f16, f128
// - complex numbers: Complex<f16>, Complex<f32>, Complex<f64>, Complex<f128>

/// Seal trait for scalar types.
pub trait Seal {}
impl Seal for bool {}
impl Seal for i8 {}
impl Seal for i16 {}
impl Seal for i32 {}
impl Seal for i64 {}
impl Seal for i128 {}
impl Seal for isize {}
impl Seal for u8 {}
impl Seal for u16 {}
impl Seal for u32 {}
impl Seal for u64 {}
impl Seal for u128 {}
impl Seal for usize {}
impl Seal for f32 {}
impl Seal for f64 {}

/// A trait for scalar types that can be used in sparse matrices.
pub trait Scalar: Sized + Copy + Send + Sync + 'static + Seal {}
impl Scalar for bool {}
impl Scalar for i8 {}
impl Scalar for i16 {}
impl Scalar for i32 {}
impl Scalar for i64 {}
impl Scalar for i128 {}
impl Scalar for isize {}
impl Scalar for u8 {}
impl Scalar for u16 {}
impl Scalar for u32 {}
impl Scalar for u64 {}
impl Scalar for u128 {}
impl Scalar for usize {}
impl Scalar for f32 {}
impl Scalar for f64 {}
