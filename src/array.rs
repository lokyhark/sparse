use std::{
    alloc::{Layout, alloc, dealloc},
    error::Error,
    fmt::Display,
    ops::{Deref, DerefMut},
};

/// A contiguous, homogeneous, invariant and heap allocated collections of values.
#[derive(Debug)]
pub struct Array<T> {
    /// Pointer to the heap allocated array.
    ptr: *mut T,
    /// Length of the array.
    len: usize,
    /// Capacity of the array.
    cap: usize,
}

/// Errors that can occur when creating or modifying an array.
#[derive(Debug)]
pub struct ArrayError {
    kind: ArrayErrorKind,
}

/// Kinds of errors that can occur when creating or modifying an array.
#[derive(Debug)]
enum ArrayErrorKind {
    /// Zero capacity requested.
    ZeroCapacity,
    /// Zero sized type requested.
    ZeroSizedType,
    /// Type needs drop.
    TypeNeedsDrop,
    /// Layout error.
    LayoutError,
    /// Alloc error.
    AllocError,
    /// Capacity overflow.
    CapacityOverflow,
}

impl Display for ArrayError {
    fn fmt(&self, fmt: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.kind {
            ArrayErrorKind::ZeroCapacity => {
                write!(fmt, "Empty arrays are not supported.")
            }
            ArrayErrorKind::ZeroSizedType => {
                write!(fmt, "Zero sized types are not supported.")
            }
            ArrayErrorKind::TypeNeedsDrop => {
                write!(fmt, "Types that need drop are not supported.")
            }
            ArrayErrorKind::LayoutError => {
                write!(fmt, "Invalid memory layout.")
            }
            ArrayErrorKind::AllocError => {
                write!(fmt, "Memory allocation failed.")
            }
            ArrayErrorKind::CapacityOverflow => {
                write!(fmt, "Array capacity overflow.")
            }
        }
    }
}

impl Error for ArrayError {}

impl From<ArrayErrorKind> for ArrayError {
    fn from(kind: ArrayErrorKind) -> Self {
        Self { kind }
    }
}

impl<T: Clone> Array<T> {
    /// Creates a new array with the given capacity.
    pub fn new(capacity: usize) -> Result<Self, ArrayError> {
        // Empty arrays not supported.
        if capacity == 0 {
            return Err(ArrayErrorKind::ZeroCapacity.into());
        }
        // Zero sized types not supported.
        if std::mem::size_of::<T>() == 0 {
            return Err(ArrayErrorKind::ZeroSizedType.into());
        }
        // Needs drop type not supported.
        if std::mem::needs_drop::<T>() {
            return Err(ArrayErrorKind::TypeNeedsDrop.into());
        }
        // Determine allocation layout.
        let Ok(layout) = Layout::array::<T>(capacity) else {
            return Err(ArrayErrorKind::LayoutError.into());
        };
        // Allocate memory.
        let ptr = unsafe { alloc(layout) };
        // Check ptr for allocation error.
        if ptr.is_null() {
            return Err(ArrayErrorKind::AllocError.into());
        }

        let array = Self { ptr: ptr.cast(), len: 0, cap: capacity };
        Ok(array)
    }

    /// Returns array length.
    pub fn length(&self) -> usize {
        self.len
    }

    /// Returns array capacity.
    pub fn capacity(&self) -> usize {
        self.cap
    }

    /// Clears the array.
    pub fn clear(&mut self) {
        // SAFETY: Entries do not need drop.
        self.len = 0;
    }

    /// Pushes a value to the array.
    ///
    /// # Errors
    ///
    /// Returns an error if the array is full.
    pub fn push(&mut self, element: T) -> Result<(), ArrayError> {
        if self.len == self.cap {
            return Err(ArrayErrorKind::CapacityOverflow.into());
        }
        // SAFETY: We have exclusive access to the array and we are within bounds.
        unsafe {
            self.ptr.add(self.len).write(element);
        }
        self.len += 1;
        Ok(())
    }

    pub fn resize(&mut self, length: usize, value: T) -> Result<(), ArrayError> {
        if length > self.cap {
            return Err(ArrayErrorKind::CapacityOverflow.into());
        }
        self.len = length;
        self.fill(value);
        Ok(())
    }
}

impl<T: Clone> Clone for Array<T> {
    fn clone(&self) -> Self {
        let mut clone = Self::new(self.cap).expect("failed to clone array");
        // SAFETY: We have exclusive access to the clone and capacity >= length.
        unsafe {
            std::ptr::copy_nonoverlapping(self.ptr, clone.ptr, self.len);
        }
        clone.len = self.len;
        clone
    }
}

impl<T> Drop for Array<T> {
    fn drop(&mut self) {
        // SAFETY: ptr is non-null and valid for len elements within a single allocation.
        unsafe {
            let layout = Layout::array::<T>(self.cap).unwrap();
            dealloc(self.ptr.cast(), layout);
        }
    }
}

impl<T> Deref for Array<T> {
    type Target = [T];

    fn deref(&self) -> &Self::Target {
        // SAFETY: ptr is non-null and valid for len elements within a single allocation.
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }
}

impl<T> DerefMut for Array<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        // SAFETY: ptr is non-null and valid for len elements within a single allocation.
        unsafe { std::slice::from_raw_parts_mut(self.ptr, self.len) }
    }
}
