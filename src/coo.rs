use crate::{Index, Scalar, array::Array};

/// Coordinate format (COO) sparse pattern.
pub type CooPat<I> = CoordinatePattern<I>;
/// Coordinate format (COO) sparse matrix.
pub type CooMat<I, S> = CoordinateMatrix<I, S>;

/// Coordinate format (COO) sparse pattern.
///
/// # Format
///
/// The coordinate format (COO) stores tuples `(row, col)` for each non zero entry in the pattern.
///
/// # Properties
///
/// - The coordinate format is intended for incremental pattern construction.
/// - The coordinate format *allows* duplicate entries.
///
/// # Storage
///
/// Let `A` a 3-by-4 pattern with 6 non-zero entries:
///
/// ```text
///     | 0 X 0 X |
/// A = | 0 0 X 0 |
///     | X X X 0 |
/// ```
/// In arbitrary order, 'A' may be stored as follows:
///
/// ```text
/// (row, col) idx
/// (  0,   1) [0]
/// (  1,   2) [1]
/// (  2,   1) [2]
/// (  0,   3) [3]
/// (  2,   0) [4]
/// (  2,   2) [5]
/// ```
///
/// In row major order, `A` is stored as follows:
///
/// ```text
/// (row, col) idx
/// (  0,  1)  [0]
/// (  0,  3)  [1]
/// (  1,  2)  [2]
/// (  2,  0)  [3]
/// (  2,  1)  [4]
/// (  2,  2)  [5]
///    ^
///    | major
/// ```
///
/// In column major order, `A` is stored as follows:
///
/// ```text
/// (row, col) idx
/// (  2,   0) [0]
/// (  0,   1) [1]
/// (  2,   1) [2]
/// (  1,   2) [3]
/// (  2,   2) [4]
/// (  0,   3) [5]
///    ^
///    | major
/// ```
///
/// # Methods overview
///
/// | Method                                                  | Description                                                             |
/// |---------------------------------------------------------|-------------------------------------------------------------------------|
/// | [`new(nrows, ncols, capacity)`](CoordinatePattern::new) | Creates a new coordinate pattern with the given dimensions and capacity.|
/// | [`push(row, col)`](CoordinatePattern::push)             | Push a new entry to the pattern.                                        |
/// | [`iter()`](CoordinatePattern::iter)                     | Returns an iterator over all entries in the pattern.                    |
#[derive(Clone, Debug)]
pub struct CoordinatePattern<I: Index> {
    /// Number of rows.
    nrows: I,
    /// Number of columns.
    ncols: I,
    /// Nonzero entries.
    entries: Array<(I, I)>,
}

impl<I: Index> CoordinatePattern<I> {
    /// Creates a new empty coordinate pattern with the given dimensions and capacity.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinatePattern;
    ///
    /// // Create a new coordinate pattern with 3 rows, 4 columns and capacity for 10 entries.
    /// let mut pat = CoordinatePattern::new(3, 4, 10);
    /// ```
    ///
    /// # Panics
    /// Panics if:
    /// - `nrows` is not positive.
    /// - `ncols` is not positive.
    /// - `capacity` is not positive.
    /// - memory allocation fails.
    pub fn new(nrows: I, ncols: I, capacity: usize) -> Self {
        if nrows <= I::zero() {
            panic!("number of rows must be positive");
        }
        if ncols <= I::zero() {
            panic!("number of columns must be positive");
        }
        if capacity == 0 {
            panic!("capacity must be positive");
        }
        match Array::new(capacity) {
            Ok(entries) => Self { nrows, ncols, entries },
            Err(err) => panic!("failed to create coordinate pattern: {err}"),
        }
    }

    /// Returns the number of rows.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinatePattern;
    ///
    /// let pat = CoordinatePattern::new(3, 4, 10);
    /// assert_eq!(pat.nrows(), 3);
    /// ```
    pub fn nrows(&self) -> I {
        self.nrows
    }

    /// Returns the number of columns.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinatePattern;
    ///
    /// let pat = CoordinatePattern::new(3, 4, 10);
    /// assert_eq!(pat.ncols(), 4);
    /// ```
    pub fn ncols(&self) -> I {
        self.ncols
    }

    /// Returns the number of entries in the pattern.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinatePattern;
    ///
    /// let mut pat = CoordinatePattern::new(3, 4, 10);
    /// assert_eq!(pat.length(), 0);
    /// pat.push(0, 1);
    /// pat.push(1, 2);
    /// assert_eq!(pat.length(), 2);
    /// ```
    pub fn length(&self) -> usize {
        self.entries.length()
    }

    /// Returns the capacity of the pattern.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinatePattern;
    ///
    /// let pat = CoordinatePattern::new(3, 4, 10);
    /// assert_eq!(pat.capacity(), 10);
    /// ```
    pub fn capacity(&self) -> usize {
        self.entries.capacity()
    }

    /// Pushes a new entry to the pattern.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinatePattern;
    ///
    /// let mut pat = CoordinatePattern::new(3, 4, 10);
    /// pat.push(0, 1);
    /// pat.push(1, 2);
    /// assert_eq!(pat.length(), 2);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - `row` is out of bounds.
    /// - `col` is out of bounds.
    /// - capacity overflow.
    pub fn push(&mut self, row: I, col: I) {
        if row >= self.nrows {
            panic!("row index `{}` out of bounds; must be in [{};{})", row, 0, self.nrows());
        }
        if col >= self.ncols {
            panic!("column index `{}` out of bounds; must be in [{};{})", col, 0, self.ncols());
        }
        match self.entries.push((row, col)) {
            Ok(()) => (),
            Err(_) => panic!("failed to push entry to pattern"),
        }
    }

    /// Clears the pattern.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinatePattern;
    ///
    /// let mut pat = CoordinatePattern::new(3, 4, 10);
    /// pat.push(0, 1);
    /// pat.push(1, 2);
    /// assert_eq!(pat.length(), 2);
    /// pat.clear();
    /// assert_eq!(pat.length(), 0);
    /// ```
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// Returns an iterator over all entries in the pattern.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinatePattern;
    ///
    /// let mut pat = CoordinatePattern::new(3, 4, 10);
    /// pat.push(0, 1);
    /// pat.push(1, 2);
    ///
    /// let mut entries = pat.iter();
    /// assert_eq!(entries.next(), Some((0, 1)));
    /// assert_eq!(entries.next(), Some((1, 2)));
    /// assert_eq!(entries.next(), None);
    /// ```
    pub fn iter(&self) -> impl Iterator<Item = (I, I)> {
        self.entries.iter().map(|&(row, col)| (row, col))
    }
}

/// Coordinate format (COO) sparse matrix.
///
/// # Format
///
/// The coordinate format (COO) stores tuples `(row, col, value)` for each non zero entry in the matrix.
///
/// # Properties
///
/// - The coordinate format is intended for incremental matrix construction.
/// - The coordinate format *allows* duplicate entries.
/// - When converted to other compressed formats, duplicate entries are summed together.
///
/// # Storage
///
/// Let `A` a 3-by-4 matrix with 6 non-zero entries:
///
/// ```text
///     | 0 1 0 2 |
/// A = | 0 0 3 0 |
///     | 4 5 6 0 |
/// ```
///
/// In arbitrary order, 'A' may be stored as follows:
///
/// ```text
/// (row, col, val) idx
/// (  0,   1,  1.) [0]
/// (  1,   2,  3.) [1]
/// (  2,   1,  5.) [2]
/// (  0,   3,  2.) [3]
/// (  2,   0,  4.) [4]
/// (  2,   2,  6.) [5]
/// ```
///
/// In row major order, `A` is stored as follows:
///
/// ```text
/// (row, col, val) idx
/// (  0,   1,  1.) [0]
/// (  0,   3,  2.) [1]
/// (  1,   2,  3.) [2]
/// (  2,   0,  4.) [3]
/// (  2,   1,  5.) [4]
/// (  2,   2,  6.) [5]
///    ^
///    | major
/// ```
///
/// In column major order, `A` is stored as follows:
///
/// ```text
/// (row, col, val) idx
/// (  2,   0, 4.0) [0]
/// (  0,   1, 1.0) [1]
/// (  2,   1, 5.0) [2]
/// (  1,   2, 3.0) [3]
/// (  2,   2, 6.0) [4]
/// (  0,   3, 2.0) [5]
///         ^
///         | major
/// ```
///
/// # Methods overview
///
/// | Method                                                  | Description                                                                             |
/// |---------------------------------------------------------|-----------------------------------------------------------------------------------------|
/// | [`new(nrows, ncols, capacity)`](CoordinateMatrix::new)  | Creates a new coordinate matrix with the given dimensions and capacity.                 |
/// | [`push(row, col, value)`](CoordinateMatrix::push)       | Push a new entry to the matrix.                                                         |
/// | [`iter()`](CoordinateMatrix::iter)                      | Returns an iterator over all entries in the matrix with shared reference to values.     |
/// | [`iter_mut()`](CoordinateMatrix::iter_mut)              | Returns an iterator over all entries in the matrix with exclusive reference to values.  |
#[derive(Clone, Debug)]
pub struct CoordinateMatrix<I: Index, S: Scalar> {
    /// Number of rows.
    nrows: I,
    /// Number of columns.
    ncols: I,
    /// Nonzero entries.
    entries: Array<(I, I, S)>,
}

impl<I: Index, S: Scalar> CoordinateMatrix<I, S> {
    /// Creates a new empty coordinate matrix with the given dimensions and capacity.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinateMatrix;
    ///
    /// // Create a new coordinate matrix with 3 rows, 4 columns and capacity for 10 entries.
    /// let mut mat: CoordinateMatrix<_, f64> = CoordinateMatrix::new(3, 4, 10);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - `nrows` is not positive.
    /// - `ncols` is not positive.
    /// - `capacity` is not positive.
    /// - Memory allocation fails.
    pub fn new(nrows: I, ncols: I, capacity: usize) -> Self {
        if nrows <= I::zero() {
            panic!("number of rows must be positive");
        }
        if ncols <= I::zero() {
            panic!("number of columns must be positive");
        }
        if capacity == 0 {
            panic!("capacity must be positive");
        }
        match Array::new(capacity) {
            Ok(entries) => Self { nrows, ncols, entries },
            Err(err) => panic!("failed to create coordinate pattern: {err}"),
        }
    }

    /// Returns the number of rows.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinateMatrix;
    ///
    /// let mat: CoordinateMatrix<_, f64> = CoordinateMatrix::new(3, 4, 10);
    /// assert_eq!(mat.nrows(), 3);
    /// ```
    pub fn nrows(&self) -> I {
        self.nrows
    }

    /// Returns the number of columns.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinateMatrix;
    ///
    /// let mat: CoordinateMatrix<_, f64> = CoordinateMatrix::new(3, 4, 10);
    /// assert_eq!(mat.ncols(), 4);
    /// ```
    pub fn ncols(&self) -> I {
        self.ncols
    }

    /// Returns the number of entries in the matrix.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinateMatrix;
    ///
    /// let mut mat = CoordinateMatrix::new(3, 4, 10);
    /// assert_eq!(mat.length(), 0);
    /// mat.push(0, 1, 1.0);
    /// mat.push(1, 2, 3.0);
    /// assert_eq!(mat.length(), 2);
    /// ```
    pub fn length(&self) -> usize {
        self.entries.len()
    }

    /// Returns the capacity of the matrix.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinateMatrix;
    ///
    /// let mat: CoordinateMatrix<_, f64> = CoordinateMatrix::new(3, 4, 10);
    /// assert_eq!(mat.capacity(), 10);
    /// ```
    pub fn capacity(&self) -> usize {
        self.entries.capacity()
    }

    /// Clears the matrix.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinateMatrix;
    ///
    /// let mut mat = CoordinateMatrix::new(3, 4, 10);
    /// mat.push(0, 1, 1.0);
    /// mat.push(1, 2, 3.0);
    /// assert_eq!(mat.length(), 2);
    /// mat.clear();
    /// assert_eq!(mat.length(), 0);
    /// ```
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// Push new entry to the matrix.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinateMatrix;
    ///
    /// let mut mat = CoordinateMatrix::new(3, 4, 10);
    /// mat.push(0, 1, 1.0);
    /// mat.push(1, 2, 3.0);
    /// assert_eq!(mat.length(), 2);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - `row` is out of bounds.
    /// - `col` is out of bounds.
    /// - capacity overflow.
    pub fn push(&mut self, row: I, col: I, value: S) {
        if row >= self.nrows {
            panic!("row index `{}` out of bounds; must be in [{};{})", row, 0, self.nrows());
        }
        if col >= self.ncols {
            panic!("column index `{}` out of bounds; must be in [{};{})", col, 0, self.ncols());
        }
        match self.entries.push((row, col, value)) {
            Ok(()) => (),
            Err(_) => panic!("failed to push entry to matrix"),
        }
    }

    /// Returns an iterator over all entries in the matrix with shared reference to values.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinateMatrix;
    ///
    /// let mut mat = CoordinateMatrix::new(3, 4, 10);
    /// mat.push(0, 1, 1.0);
    /// mat.push(1, 2, 3.0);
    ///
    /// let mut entries = mat.iter();
    /// assert_eq!(entries.next(), Some((0, 1, &1.0)));
    /// assert_eq!(entries.next(), Some((1, 2, &3.0)));
    /// assert_eq!(entries.next(), None);
    /// ```
    pub fn iter(&self) -> impl Iterator<Item = (I, I, &'_ S)> {
        self.entries.iter().map(|(row, col, val)| (*row, *col, val))
    }

    /// Returns an iterator over all entries in the matrix with exclusive reference to values.
    ///
    /// # Examples
    ///
    /// ```
    /// use sparse::CoordinateMatrix;
    ///
    /// let mut mat = CoordinateMatrix::new(3, 4, 10);
    /// mat.push(0, 1, 1.0);
    /// mat.push(1, 2, 3.0);
    ///
    /// let mut entries = mat.iter_mut();
    /// assert_eq!(entries.next(), Some((0, 1, &mut 1.0)));
    /// assert_eq!(entries.next(), Some((1, 2, &mut 3.0)));
    /// assert_eq!(entries.next(), None);
    /// ```
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (I, I, &'_ mut S)> {
        self.entries.iter_mut().map(|(row, col, val)| (*row, *col, val))
    }
}
