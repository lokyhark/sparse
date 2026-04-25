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
pub trait Index: Sized + Copy + Send + Sync + 'static + Seal {}
impl Index for u8 {}
impl Index for u16 {}
impl Index for u32 {}
impl Index for u64 {}
impl Index for u128 {}
impl Index for usize {}
impl Index for i8 {}
impl Index for i16 {}
impl Index for i32 {}
impl Index for i64 {}
impl Index for i128 {}
impl Index for isize {}
