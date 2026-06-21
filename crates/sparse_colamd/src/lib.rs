//! COLAMD Algorithm module.
//!
//! The main entry point of the COLAMD algorithm is the [`colamd`] function.
//!
//! ```
//! # use sparse_colamd::colamd;
//!
//! // Matrix  5-by-4
//! // x 0 x 0
//! // x 0 x x
//! // 0 x x 0
//! // 0 0 x x
//! // x x 0 0
//! let rowind = [0, 1, 4, 2, 4, 0, 1, 2, 3, 1, 3];
//! let colptr = [0, 3, 5, 9, 11];
//! // Execute COLAMD algorithm.
//! let result = colamd(5, 4, &colptr, &rowind).unwrap();
//! // Retrieve ordering.
//! let order = result.order();
//! assert_eq!(order[0], 1);
//! assert_eq!(order[1], 0);
//! assert_eq!(order[2], 2);
//! assert_eq!(order[3], 3);
//! ```
//!
//! If full control over allocation and potential workspace reuse is possible,
//! maximum performances can be reached with [`Colamd`] workspace:
//!
//! ```
//! # use sparse_colamd::Colamd;
//! # use sparse_colamd::ColamdError;
//! # fn main() -> Result<(), ColamdError<i32>> {
//! // Allocate workspace for 5 rows, 4 columns, 11 non zero entries and 10 elbow room capacity.
//! let mut colamd = Colamd::alloc(5, 4, 11, 10)?;
//! // Determine ordering.
//! let order = colamd.run(5, 4, &[0, 3, 5, 9, 11], &[0, 1, 4, 2, 4, 0, 1, 2, 3, 1, 3])?;
//! // The workspace can be reused without any reallocation.
//! let order = colamd.run(5, 4, &[0, 3, 5, 9, 11], &[0, 1, 4, 2, 4, 0, 1, 2, 3, 1, 3])?;
//! # Ok(())
//! # }
//! ```
//!
//! # References
//!
//! - T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng,<br />
//!   *An approximate column minimum degree ordering algorithm*,<br />
//!   ACM Transactions on Mathematical Software, vol. 30, no. 3., pp. 353-376, 2004.<br />
//!   <https://doi.org/10.1145/1024074.1024079>
//! - T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng,<br />
//!   *Algorithm 836: COLAMD, an approximate column minimum degree ordering algorithm*,<br />
//!   ACM Transactions on Mathematical Software, vol. 30, no. 3., pp. 377-380, 2004.<br />
//!   <https://doi.org/10.1145/1024074.1024080>
use sparse_array::Array;

mod col;
mod config;
mod error;
mod int;
mod result;
mod row;
mod score;
mod size;

use col::ColamdCol;
use config::ColamdConfig;
use error::ColamdErrorKind;
use row::ColamdRow;
use score::ColamdScore;
use size::ColamdSize;

pub use error::ColamdError;
pub use int::ColamdInt;
pub use result::ColamdResult;

/// COLAMD algorithm.
///
/// # Errors
///
/// Return a [`ColamdError`] if one of the following conditions hold:
/// - Column pointers slice is empty (`colptr.is_empty() == true`),
/// - Column pointers length is invalid (`colptr.len() != ncols + 1`),
/// - Column pointers slice first element is not zero (`colptr[0] != 0`),
/// - Column pointers is not in increasing order,
///
/// # References
///
/// - T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng,<br />
///   *An approximate column minimum degree ordering algorithm*,<br />
///   ACM Transactions on Mathematical Software, vol. 30, no. 3., pp. 353-376, 2004.<br />
///   <https://doi.org/10.1145/1024074.1024079>
/// - T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng,<br />
///   *Algorithm 836: COLAMD, an approximate column minimum degree ordering algorithm*,<br />
///   ACM Transactions on Mathematical Software, vol. 30, no. 3., pp. 377-380, 2004.<br />
///   <https://doi.org/10.1145/1024074.1024080>
pub fn colamd<I: ColamdInt>(nrows: I, ncols: I, colptr: &[I], rowind: &[I]) -> Result<ColamdResult<I>, ColamdError<I>> {
    // Check inputs.
    let size = Colamd::check(nrows, ncols, colptr, rowind)?;
    // Allocate algorithm.
    let elbow = size.ncols.as_usize().saturating_add(size.nnz.as_usize() / 5);
    let mut colamd: Colamd<I> = Colamd::alloc(size.nrows.as_usize(), size.ncols.as_usize(), size.nnz.as_usize(), elbow)?;
    // Execute algorithm.
    colamd.exec(size.nrows, size.ncols, size.nnz, colptr, rowind)
}

/// COLAMD algorithm workspace.
///
/// # Examples
///
/// ```
/// # use sparse_colamd::Colamd;
/// # use sparse_colamd::ColamdError;
/// # fn main() -> Result<(), ColamdError<i32>> {
/// let mut colamd = Colamd::alloc(5, 4, 11, 10)?;
/// let order = colamd.run(5, 4, &[0, 3, 5, 9, 11], &[0, 1, 4, 2, 4, 0, 1, 2, 3, 1, 3])?;
/// # Ok(())
/// # }
/// ```
///
/// # See Also
///
/// [`colamd`]
#[derive(Debug)]
pub struct Colamd<I: ColamdInt> {
    /// COLAMD algorithm configuration.
    config: ColamdConfig<I>,
    /// Row form matrix representation.
    rows: Array<ColamdRow<I>>,
    /// Column form matrix representation.
    cols: Array<ColamdCol<I>>,
    /// Row/Column/Pivot indices.
    inds: Array<I>,
    /// Degree list.
    degree: Array<I>,
}

/// Minimum index value.
const MIN: usize = 0;
/// Maximum index value.
const MAX: usize = isize::MAX as usize;

impl<I: ColamdInt> Colamd<I> {
    const EMPTY: I = I::NEG;

    /// Allocate COLAMD workspace for a matrix of size `nrows` by `ncols` with `nnz` non zeros.
    ///
    /// # Errors
    ///
    /// Return a [`ColamdError`] if one of the following conditions hold:
    /// - Number of rows is zero (`nrows == 0`),
    /// - Number of columns is zero (`ncols == 0`).
    /// - Elbow room is too short (`elbow < ncols`).
    /// - Workspace memory allocation fails.
    pub fn alloc(nrows: usize, ncols: usize, nnz: usize, elbow: usize) -> Result<Self, ColamdError<I>> {
        // Check number of rows.
        if nrows == 0 {
            return Err(ColamdErrorKind::InvalidNumberOfRows { nrows: I::ZERO }.into());
        }
        // Check number of columns.
        if ncols == 0 {
            return Err(ColamdErrorKind::InvalidNumberOfColumns { ncols: I::ZERO }.into());
        }
        // Check number of non-zeros.
        if nnz == 0 {
            return Err(ColamdErrorKind::EmptyMatrix.into());
        }
        // Check elbow room.
        if elbow < ncols {
            return Err(ColamdErrorKind::InvalidElbowRoom { size: elbow, ncols }.into());
        }

        // Compute indices capacity.
        let cap = nnz.saturating_mul(2).saturating_add(elbow);

        // Allocate workspace.
        let rows = Array::new(nrows)?;
        let cols = Array::new(ncols)?;
        let inds = Array::new(cap)?;
        let degree = Array::new(ncols + 1)?;

        Ok(Self {
            config: ColamdConfig::default(),
            rows,
            cols,
            inds,
            degree,
        })
    }

    /// Set the COLAMD configuration.
    #[must_use]
    pub fn with_config(mut self, config: ColamdConfig<I>) -> Self {
        self.config = config;
        self
    }

    /// Return a shared reference to the column at `index`.
    #[inline]
    fn col(&self, index: I) -> &ColamdCol<I> {
        &self.cols[index.as_usize()]
    }

    /// Return an exclusive reference to the column at `index`.
    #[inline]
    fn col_mut(&mut self, index: I) -> &mut ColamdCol<I> {
        &mut self.cols[index.as_usize()]
    }

    /// Return a shared reference to the row at `index`.
    #[inline]
    fn row(&self, index: I) -> &ColamdRow<I> {
        &self.rows[index.as_usize()]
    }

    /// Return an exclusive reference to the row at `index`.
    #[inline]
    fn row_mut(&mut self, index: I) -> &mut ColamdRow<I> {
        &mut self.rows[index.as_usize()]
    }

    /// Return the row/column index stored at `index` in the indices workspace.
    #[inline]
    fn ind(&self, index: I) -> &I {
        &self.inds[index.as_usize()]
    }

    /// Return an exclusive reference to the index at `index` in the indices workspace.
    #[inline]
    fn ind_mut(&mut self, index: I) -> &mut I {
        &mut self.inds[index.as_usize()]
    }

    /// Return a shared reference to the degree list head at `score`.
    #[inline]
    fn deg(&self, score: I) -> &I {
        &self.degree[score.as_usize()]
    }

    /// Return an exclusive reference to the degree list head at `score`.
    #[inline]
    fn deg_mut(&mut self, score: I) -> &mut I {
        &mut self.degree[score.as_usize()]
    }

    /// Check COLAMD algorithm input arguments and return the required workspace size.
    ///
    /// Following proporties are checked:
    /// - `0 <= nrows < isize::MAX + 1`
    /// - `0 <= ncols < isize::MAX + 1`
    /// - `colptr.len() = ncols + 1`
    /// - `colptr[0] = 0`
    /// - `colptr[ncols + 1] = nnz`
    /// - `0 <= nnz < isize::MAX + 1`
    /// - `rowind.len() = nnz`
    ///
    /// Following properties still need to be checked:
    /// - colptr is a non-decreasing sequence
    /// - colptr pointers are valid indices
    /// - rowind indices are valid indices
    /// - rowind indices are sorted by column
    fn check(nrows: I, ncols: I, colptr: &[I], rowind: &[I]) -> Result<ColamdSize<I>, ColamdError<I>> {
        // Check matrix dimensions:
        // If dim <= 0, matrix is mathematically invalid.
        // If dim > isize::MAX, algorithm will exhaust memory: arrays ensure they never allocate more than isize::MAX bytes
        // for non zero sized types thus return instant error if future rows/cols array would exceed isize::MAX bytes.
        // Check number of rows.
        match nrows.try_into() {
            Result::<isize, _>::Err(_) => {
                return Err(ColamdErrorKind::InvalidNumberOfRows { nrows }.into());
            }
            Ok(value) if value <= 0 => {
                return Err(ColamdErrorKind::InvalidNumberOfRows { nrows }.into());
            }
            Ok(_) => (),
        }
        // Check number of columns.
        match ncols.try_into() {
            Result::<isize, _>::Err(_) => {
                return Err(ColamdErrorKind::InvalidNumberOfColumns { ncols }.into());
            }
            Ok(value) if value <= 0 => {
                return Err(ColamdErrorKind::InvalidNumberOfColumns { ncols }.into());
            }
            Ok(_) => (),
        }

        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        debug_assert!(nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}");
        debug_assert!(ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");

        // Check colptr.
        // colptr must hold multiple invariants:
        // 1. colptr length is ncols + 1 (i.e. colptr is not empty because ncols >= MIN = 0)
        // 2. colptr[0] = 0
        // 3. colptr is a non-decreasing sequence
        // 4. colptr is a memory contiguous allocation
        // 5. colptr pointers are valid indices
        // 6. colptr[ncols] = nnz (number of non zero in matrix)
        // Following code check (1.) and (2.), assumes (4.), relies on (6.), and (3.) and (5.) will be checked later on.
        if colptr.is_empty() {
            return Err(ColamdErrorKind::EmptyColptrSlice.into());
        }
        // SAFETY: ncols <= MAX = isize::MAX (see previous checkpoint) so ncols + 1 cannot overflow usize.
        if colptr.len() != ncols.as_usize() + 1 {
            return Err(ColamdErrorKind::InvalidColptrSliceLength {
                actual: colptr.len(),
                expected: ncols.as_usize() + 1,
            }
            .into());
        }

        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // colptr.len() = ncols + 1
        debug_assert!(nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}");
        debug_assert!(ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert_eq!(colptr.len(), ncols.as_usize() + 1, "invalid col pointers slice length: actual={}, expect={}", colptr.len(), ncols.as_usize() + 1);

        // Check first column pointer.
        // SAFETY: colptr is not empty (see previous checkpoint).
        if colptr[0] != I::ZERO {
            return Err(ColamdErrorKind::InvalidFirstColPtr { ptr: colptr[0] }.into());
        }

        // Check last column pointer which is the number on non zeros.
        // SAFETY: colptr has correct length (see previous checkpoint).
        let nnz = colptr[ncols.as_usize()];
        match nnz.try_into() {
            Result::<isize, _>::Err(_) => {
                return Err(ColamdErrorKind::InvalidColPtr {
                    index: ncols.as_usize(),
                    ptr: colptr[ncols.as_usize()],
                }
                .into());
            }
            Ok(0) => return Err(ColamdErrorKind::EmptyMatrix.into()),
            Ok(value) if value < 0 => {
                return Err(ColamdErrorKind::InvalidColPtr {
                    index: ncols.as_usize(),
                    ptr: colptr[ncols.as_usize()],
                }
                .into());
            }
            Ok(_) => (),
        }

        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // colptr.len() = ncols + 1
        // colptr[0] = 0
        // colptr[ncols + 1] = nnz
        // 0 <= nnz < isize::MAX + 1
        debug_assert!(nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}");
        debug_assert!(ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert_eq!(colptr.len(), ncols.as_usize() + 1, "invalid col pointers slice length: actual={}, expect={}", colptr.len(), ncols.as_usize() + 1);
        debug_assert_eq!(colptr[0], I::ZERO, "invalid first column pointer: actual={}, expect=0", colptr[0]);
        debug_assert_eq!(colptr[ncols.as_usize()], nnz, "invalid number of nnz: actual={}, expect={nnz}", colptr[ncols.as_usize()]);
        debug_assert!(MIN < nnz.as_usize() && nnz.as_usize() <= MAX, "invalid number of non zeros: nnz={nnz}, min={MIN}, max={MAX}");

        // Check rowind.
        // rowind must hold multiple invariants:
        // 1. rowind length is nnz(=colptr[ncols])
        // 2. rowind indices are valid indices
        // 3. rowind is a memory continuous allocation
        // Following code checks (1.), assumes (3.) and (2.) will be checked later on.
        if rowind.len() != nnz.as_usize() {
            return Err(ColamdErrorKind::InvalidRowindSliceLength {
                actual: rowind.len(),
                expected: nnz.as_usize(),
            }
            .into());
        }

        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // colptr.len() = ncols + 1
        // colptr[0] = 0
        // colptr[ncols + 1] = nnz
        // 0 <= nnz < isize::MAX + 1
        // rowind.len() = nnz
        debug_assert!(nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}");
        debug_assert!(ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert_eq!(colptr.len(), ncols.as_usize() + 1, "invalid col pointers slice length: actual={}, expect={}", colptr.len(), ncols.as_usize() + 1);
        debug_assert_eq!(colptr[0], I::ZERO, "invalid first column pointer: actual={}, expect=0", colptr[0]);
        debug_assert_eq!(colptr[ncols.as_usize()], nnz, "invalid number of nnz: actual={}, expect={nnz}", colptr[ncols.as_usize()]);
        debug_assert!(MIN < nnz.as_usize() && nnz.as_usize() <= MAX, "invalid number of non zeros: nnz={nnz}, min={MIN}, max={MAX}");
        debug_assert_eq!(rowind.len(), nnz.as_usize(), "invalid row indices slice length: actual={}, expect={nnz}", rowind.len());

        // At this point, following still need to be checked:
        // - colptr is a non-decreasing sequence
        // - colptr pointers are valid indices
        // - rowind indices are valid indices
        // - rowind indices are sorted by column

        Ok(ColamdSize::new(nrows, ncols, nnz))
    }

    /// Run the COLAMD algorithm on a matrix of size `nrows` by `ncols` with `nnz` non zeros represented
    /// in column compressed form by `colptr` and `rowind`.
    ///
    /// # Errors
    ///
    /// On input following invariants must hold or an error is returned:
    /// - `0 <= nrows < isize::MAX + 1`
    /// - `0 <= ncols < isize::MAX + 1`
    /// - `colptr.len() = ncols + 1`
    /// - `colptr[0] = 0`
    /// - `colptr[ncols + 1] = nnz`
    /// - `0 <= nnz < isize::MAX + 1`
    /// - `rowind.len() = nnz`
    pub fn run(&mut self, nrows: I, ncols: I, colptr: &[I], rowind: &[I]) -> Result<ColamdResult<I>, ColamdError<I>> {
        // Check input.
        let size = Colamd::check(nrows, ncols, colptr, rowind)?;
        // Execute algorithm.
        self.exec(size.nrows, size.ncols, size.nnz, colptr, rowind)
    }

    /// Execute the COLAMD algorithm.
    fn exec(&mut self, nrows: I, ncols: I, nnz: I, colptr: &[I], rowind: &[I]) -> Result<ColamdResult<I>, ColamdError<I>> {
        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // colptr.len() = ncols + 1
        // colptr[0] = 0
        // colptr[ncols + 1] = nnz
        // 0 <= nnz < isize::MAX + 1
        // rowind.len() = nnz
        debug_assert!(nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}");
        debug_assert!(ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert_eq!(colptr.len(), ncols.as_usize() + 1, "invalid col pointers slice length: actual={}, expect={}", colptr.len(), ncols.as_usize() + 1);
        debug_assert_eq!(colptr[0], I::ZERO, "invalid first column pointer: actual={}, expect=0", colptr[0]);
        debug_assert_eq!(colptr[ncols.as_usize()], nnz, "invalid number of nnz: actual={}, expect={nnz}", colptr[ncols.as_usize()]);
        debug_assert!(MIN < nnz.as_usize() && nnz.as_usize() <= MAX, "invalid number of non zeros: nnz={nnz}, min={MIN}, max={MAX}");
        debug_assert_eq!(rowind.len(), nnz.as_usize(), "invalid row indices slice length: actual={}, expect={nnz}", rowind.len());

        // Initialize rows and cols.
        self.init(nrows, ncols, nnz, colptr, rowind)?;

        // Score rows and cols.
        let score = self.score(nrows, ncols, nnz)?;

        // Find order.
        self.find(nrows, ncols, nnz, score.cols, score.max_degree)?;

        // Order.
        let order = self.order(ncols)?;

        Ok(ColamdResult {
            cols: score.cols,
            rows: score.rows,
            max_degree: score.max_degree,
            min_score: score.min_score,
            order,
        })
    }

    /// Initialize the COLAMD algorithm.
    ///
    /// On input following invariants must hold:
    /// - `0 <= nrows < isize::MAX + 1`
    /// - `0 <= ncols < isize::MAX + 1`
    /// - `colptr.len() = ncols + 1`
    /// - `colptr[0] = 0`
    /// - `colptr[ncols + 1] = nnz`
    /// - `0 <= nnz < isize::MAX + 1`
    /// - `rowind.len() = nnz`
    ///
    /// On output following invariants hold:
    /// - `0 <= nrows < isize::MAX + 1`
    /// - `0 <= ncols < isize::MAX + 1`
    /// - `colptr.len() = ncols + 1`
    /// - `colptr[0] = 0`
    /// - `colptr[ncols + 1] = nnz`
    /// - `0 <= nnz < isize::MAX + 1`
    /// - `rowind.len() = nnz`
    /// - colptr pointers are valid indices
    /// - colptr is ordered
    /// - rowind indices are valid indices
    /// - rowind is ordered by column
    /// - self.rows is allocated and initialized and contains row compressed form
    /// - self.cols is allocated, initialized and contains column compressed form
    /// - self.rowind is allocated and initialized with column compressed indices
    /// - self.colind is allocated and initialized with row compressed indices
    fn init(&mut self, nrows: I, ncols: I, nnz: I, colptr: &[I], rowind: &[I]) -> Result<(), ColamdError<I>> {
        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // colptr.len() = ncols + 1
        // colptr[0] = 0
        // colptr[ncols + 1] = nnz
        // 0 <= nnz < isize::MAX + 1
        // rowind.len() = nnz
        debug_assert!(nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}");
        debug_assert!(ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert_eq!(colptr.len(), ncols.as_usize() + 1, "invalid col pointers slice length: actual={}, expect={}", colptr.len(), ncols.as_usize() + 1);
        debug_assert_eq!(colptr[0], I::ZERO, "invalid first column pointer: actual={}, expect=0", colptr[0]);
        debug_assert_eq!(colptr[ncols.as_usize()], nnz, "invalid number of nnz: actual={}, expect={nnz}", colptr[ncols.as_usize()]);
        debug_assert!(MIN < nnz.as_usize() && nnz.as_usize() <= MAX, "invalid number of non zeros: nnz={nnz}, min={MIN}, max={MAX}");
        debug_assert_eq!(rowind.len(), nnz.as_usize(), "invalid row indices slice length: actual={}, expect={nnz}", rowind.len());

        // Columns do not need to be initialized because they will be accessed in increasing order below
        // and pushed one by one.
        // The user may reuse a previously allocated workspace so clear the columns array.
        self.cols.clear();

        // Initialize columns.
        // SAFETY: colptr[0] = 0 previously checked
        let mut start = I::ZERO;
        for col in 0..ncols.as_usize() {
            // SAFETY: colptr length was previously checked.
            let stop = colptr[col + 1];
            // Check if the column pointer is valid. It must have following properties:
            // 1. <= nnz
            // 2. >= 0
            // 3. non decreasing
            // (2.) is implied by (3.) because colptr[0] = 0 === start
            // Check maximum.
            if stop > nnz {
                return Err(ColamdErrorKind::InvalidColPtr { index: col + 1, ptr: stop }.into());
            }
            // Check increasing order.
            if stop < start {
                return Err(ColamdErrorKind::DecreasingColptrSlice { col, start, stop }.into());
            }
            // Column length.
            // SAFETY: cannot overflow because stop >= start and both are valid index (previously checked).
            let length = stop - start;
            self.cols.push(ColamdCol {
                start,
                length,
                weight: I::ONE,
                rank: I::ZERO,
                prev: Self::EMPTY,
                next: Self::EMPTY,
            })?;
            // Set start for next iteration.
            start = stop;
        }

        // Resize and initialize rows.
        // Rows must be initialized because they will be accessed in random order below.
        let default = ColamdRow {
            start,
            length: I::ZERO,
            degree: Self::EMPTY,
            mark: I::ZERO,
        };
        self.rows.resize(nrows.as_usize(), default)?;
        // Use column form.
        for j in I::range(I::ZERO, ncols) {
            let col = self.col(j);
            let start = col.start;
            // SAFETY: column length was previously checked.
            let stop = start + col.length;
            // Cache previous row to detect unsorted columns.
            let mut prev = I::NEG;
            // At this point, row indices may be invalid (out of bounds), out of order or present multiple times.
            for ptr in I::range(start, stop) {
                // SAFETY: rowind and colptr lengths were previously checked.
                let ptr = ptr.as_usize();
                let row = rowind[ptr];
                // Check bounds.
                if row >= nrows || row < I::ZERO {
                    return Err(ColamdErrorKind::OutOfBoundRowIndex { index: ptr, row, nrows }.into());
                }
                // Check if row indices in column is ordered and unique.
                if row <= prev {
                    return Err(ColamdErrorKind::UnorderedRowIndSlice { col: j, ptr }.into());
                }
                prev = row;

                // SAFETY: row is valid index because >= 0 and < nrows which is a valid index.
                let row = self.row_mut(row);
                // SAFETY: cannot overflow nrows which is a valid index.
                row.length += I::ONE;
            }
        }

        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // colptr.len() = ncols + 1
        // colptr[0] = 0
        // colptr[ncols + 1] = nnz
        // 0 <= nnz < isize::MAX + 1
        // rowind.len() = nnz
        // colptr pointers are valid indices
        // colptr is ordered
        // rowind indices are valid indices
        // rowind is ordered by column
        debug_assert!(nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}");
        debug_assert!(ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert_eq!(colptr.len(), ncols.as_usize() + 1, "invalid col pointers slice length: actual={}, expect={}", colptr.len(), ncols.as_usize() + 1);
        debug_assert_eq!(colptr[0], I::ZERO, "invalid first column pointer: actual={}, expect=0", colptr[0]);
        debug_assert_eq!(colptr[ncols.as_usize()], nnz, "invalid number of nnz: actual={}, expect={nnz}", colptr[ncols.as_usize()]);
        debug_assert!(MIN < nnz.as_usize() && nnz.as_usize() <= MAX, "invalid number of non zeros: nnz={nnz}, min={MIN}, max={MAX}");
        debug_assert_eq!(rowind.len(), nnz.as_usize(), "invalid row indices slice length: actual={}, expect={nnz}", rowind.len());
        debug_assert!(colptr.iter().all(|&i| I::ZERO <= i && i <= nnz), "invalid column pointers");
        debug_assert!(colptr.is_sorted(), "unsorted column pointers");
        debug_assert!(rowind.iter().all(|&i| I::ZERO <= i && i < nrows));
        debug_assert!(colptr.array_windows().all(|[a, b]| rowind[a.as_usize()..b.as_usize()].is_sorted_by(|a, b| a < b)), "unsorted rows");

        // From this point all the invariants on input are satisfied and we can safely access cols, rows, colind and rowind.

        // Create row pointers.
        let mut ptr = nnz;
        for row in self.rows.iter_mut() {
            // Set row pointer.
            row.start = ptr;
            // Degree is used to store the last index of row in colind.
            row.degree = ptr;
            // Add row length to pointer.
            // SAFETY: ptr cannot overflow nnz and nnz is a valid index (previously checked).
            ptr += row.length;
            // Clear mark.
            row.mark = Self::EMPTY;
        }

        // Initialize column/row indices.
        self.inds.resize(2 * nnz.as_usize(), Self::EMPTY)?;
        // Initialize column indices.
        for j in I::range(I::ZERO, ncols) {
            let col = self.col(j);
            let start = col.start;
            let stop = col.start + col.length;
            for ptr in I::range(start, stop) {
                // SAFETY: row indices are valid (checked previously).
                let row = rowind[ptr.as_usize()];
                // Retrieve last index of row in colind.
                let idx = &mut self.rows[row.as_usize()].degree;
                // Store column.
                self.inds[idx.as_usize()] = j;
                // Increment last index of row.
                *idx += I::ONE;
            }
        }

        // At this point, self.inds is allocated and first nnz entries are initialized and contains row form indices.
        debug_assert_eq!(self.inds.length(), 2 * nnz.as_usize());
        debug_assert!(self.inds[nnz.as_usize()..2 * nnz.as_usize()].iter().all(|&j| j >= I::ZERO && j < ncols), "invalid column indices");

        // Copy row indices.
        // SAFETY: rowind length was previously checked.
        self.inds[0..nnz.as_usize()].copy_from_slice(rowind);

        // At this point, self.rowind is allocated, initialized and contains row form indices.
        debug_assert_eq!(self.inds.length(), 2 * nnz.as_usize());
        debug_assert!(self.inds[nnz.as_usize()..2 * nnz.as_usize()].iter().all(|&i| i >= I::ZERO && i < nrows), "invalid row indices");

        // Set row degree.
        for row in self.rows.iter_mut() {
            row.degree = row.length;
            row.mark = I::ZERO;
        }

        // At this point both row and column form is available :
        // Row form:
        // - row pointers in `self.rows`
        // - row indices in `self.rowind`
        // Column form:
        // - column pointers in `self.cols`
        // - column indices in `self.colind`

        Ok(())
    }

    /// Initialize scores of rows and columns.
    fn score(&mut self, nrows: I, ncols: I, nnz: I) -> Result<ColamdScore<I>, ColamdError<I>> {
        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // 0 <= nnz < isize::MAX + 1
        debug_assert!(nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}");
        debug_assert!(ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert!(MIN < nnz.as_usize() && nnz.as_usize() <= MAX, "invalid number of non zeros: nnz={nnz}, min={MIN}, max={MAX}");

        // Compute dense row and column thresholds.
        // SAFETY:
        // - 16 is representable by all unsigned integer types (ColamdInt which is sealed).
        // - ncols and nrows are representable by I (checked previously).
        // - dense_row_control and dense_col_control are non negative and less than nrows and ncols which are representable by I.
        let dense_row_count = Self::dense(self.config.dense_row_control, ncols);
        let dense_col_count = Self::dense(self.config.dense_col_control, nrows.min(ncols));
        // Number of columns alive.
        let mut cols = ncols;

        // Kill empty columns.
        for col in self.cols.iter_mut().rev() {
            if col.length == I::ZERO {
                cols -= I::ONE;
                col.kill();
                col.rank = cols;
            }
        }
        // Kill dense columns.
        for j in I::range(I::ZERO, ncols).rev() {
            let col = self.col(j);
            // Ignore dead column.
            if col.dead() {
                continue;
            }
            if col.length > dense_col_count {
                // Column is dense so update row degrees accordingly and kill it.
                let start = col.start;
                let stop = start + col.length;
                // Update row degrees.
                for ptr in I::range(start, stop) {
                    // SAFETY: column pointers are valid indices (checked previously).
                    let i = *self.ind(ptr);
                    // SAFETY: row indices are valid (checked previously).
                    let row = self.row_mut(i);
                    // Decrement row degree.
                    row.degree -= I::ONE;
                }
                // Kill column and place it last.
                cols -= I::ONE;
                let col = self.col_mut(j);
                col.kill();
                col.rank = cols;
            }
        }

        // Number of rows alive.
        let mut rows = nrows;
        // Maximum degree of alive rows.
        let mut max_degree = I::ZERO;

        // Kill dense and empty rows.
        for row in self.rows.iter_mut() {
            let degree = row.degree;
            debug_assert!(degree >= I::ZERO && degree <= ncols, "out of bounds degree: degree={}, min={}, max={}", degree, I::ZERO, ncols);
            if degree == I::ZERO || degree > dense_row_count {
                // Row is empty or dense so kill it.
                rows -= I::ONE;
                row.kill();
            } else {
                // Keep track of maximum row degree.
                max_degree = max_degree.max(degree);
            }
        }

        // Compute initial column scores.
        for j in I::range(I::ZERO, ncols).rev() {
            let col = self.col(j);
            // Ignore dead column.
            if col.dead() {
                continue;
            }
            let mut score = I::ZERO;
            let start = col.start;
            let stop = col.start + col.length;
            // Pointer to compact the column.
            let mut pos = col.start;
            for ptr in I::range(start, stop) {
                // SAFETY: column pointers are valid indices (checked previously).
                let i = *self.ind(ptr);
                // Ignore dead rows.
                if self.row(i).dead() {
                    continue;
                }
                // Compact the column.
                *self.ind_mut(pos) = i;
                // Update pointer.
                pos += I::ONE;
                // Prevent overflow by saturating at ncols which is valid.
                let degree = self.row(i).degree;
                if ncols - (degree - I::ONE) < score {
                    score = ncols;
                } else {
                    score += degree - I::ONE;
                }
            }
            // Compute pruned column length.
            let length = pos - start;
            // Kill column and place it last if it is now empty.
            let col = self.col_mut(j);
            if length == I::ZERO {
                cols -= I::ONE;
                col.kill();
                col.rank = cols;
            } else {
                debug_assert!(score >= I::ZERO && score <= ncols, "invalid column score: score={}, min={}, max={}", score, I::ZERO, ncols);
                col.length = length;
                col.rank = score;
            }
        }

        // Initialize degree lists.
        let mut min_score = ncols;
        // SAFETY: ncols + 1 cannot overflow usize because ncols <= MAX = isize::MAX.
        self.degree.resize(ncols.as_usize() + 1, Self::EMPTY)?;
        for j in I::range(I::ZERO, ncols).rev() {
            let col = self.col(j);
            // Ignore dead column.
            if col.dead() {
                continue;
            }
            let score = col.rank;
            debug_assert!(min_score >= I::ZERO && min_score <= ncols, "invalid minimum column score: score={}, min={}, max={}", score, I::ZERO, ncols);
            debug_assert!(score >= I::ZERO && score <= ncols, "invalid column score: score={}, min={}, max={}", score, I::ZERO, ncols);
            let head = *self.deg(score);
            let col = self.col_mut(j);
            col.prev = Self::EMPTY;
            col.next = head;
            if head != Self::EMPTY {
                self.col_mut(head).prev = j;
            }
            *self.deg_mut(score) = j;
            // Keep track of minimum score.
            min_score = min_score.min(score);
        }

        Ok(ColamdScore { cols, rows, max_degree, min_score })
    }

    fn find(&mut self, nrows: I, ncols: I, nnz: I, cols: I, mut max_degree: I) -> Result<(), ColamdError<I>> {
        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // 0 <= nnz < isize::MAX + 1
        // cols <= ncols
        // 0 < min_score <= ncols
        debug_assert!(nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}");
        debug_assert!(ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert!(MIN < nnz.as_usize() && nnz.as_usize() <= MAX, "invalid number of non zeros: nnz={nnz}, min={MIN}, max={MAX}");
        debug_assert!(I::ZERO < cols && cols <= ncols, "invalid number of alive columns: cols={}, min={}, max={}", cols, I::ZERO, ncols);

        let mut tag = I::ONE;
        let mut min_score = I::ZERO;
        let mut k = I::ZERO;

        while k < cols {
            // Choose pivot column with minimum score.
            let pivot_col_j = loop {
                let head = *self.deg(min_score);
                if head != Self::EMPTY {
                    break head;
                }
                min_score += I::ONE;
            };
            debug_assert!(pivot_col_j >= I::ZERO && pivot_col_j <= ncols);
            // Remove pivot column from degree list by placing next column as head.
            let next = self.col(pivot_col_j).next;
            *self.deg_mut(min_score) = next;
            if next != Self::EMPTY {
                self.col_mut(next).prev = Self::EMPTY;
            }
            // Retrieve pivot column.
            let pivot_col = self.col_mut(pivot_col_j);
            // Retrieve score and thickness of pivot column.
            let pivot_col_thickness = pivot_col.weight;

            // Retrieve pivot col score.
            let pivot_col_score = pivot_col.rank;
            // Place column.
            pivot_col.rank = k;
            // Increment counter by column thickness.
            k += pivot_col_thickness;

            // Compact if need room for pivot row.
            let memory = pivot_col_score.min(ncols - k);
            if memory.as_usize() > (self.inds.capacity() - self.inds.length()) {
                self.compact(nrows, ncols);
            }

            // Pivot row start point.
            let pivot_row_start: I = self.inds.length().try_into().unwrap();

            // Initialize pivot row degree.
            let mut pivot_row_degree = I::ZERO;

            // Mark pivot column.
            self.col_mut(pivot_col_j).weight = -pivot_col_thickness;

            // Create the pivot row which is the union of all rows in the pivot column.
            let start = self.col(pivot_col_j).start;
            let stop = start + self.col(pivot_col_j).length;
            for ptr in I::range(start, stop) {
                let i = *self.ind(ptr);
                let row = self.row_mut(i);
                // Ignore dead rows.
                if row.alive() {
                    let start = row.start;
                    let stop = start + row.length;
                    for ptr in I::range(start, stop) {
                        let j = *self.ind(ptr);
                        let col = self.col_mut(j);
                        let col_thickness = col.weight;
                        // Ignore dead columns and non already visited columns (thickness > 0).
                        if col.alive() && col_thickness > I::ZERO {
                            // Mark column as visited.
                            col.weight = -col_thickness;
                            // Place column in pivot row.
                            self.inds.push(j)?;
                            // Update pivot row degree.
                            pivot_row_degree += col_thickness;
                        }
                    }
                }
            }
            // Clear pivot column.
            self.col_mut(pivot_col_j).weight = pivot_col_thickness;
            max_degree = max_degree.max(pivot_row_degree);

            // Kill rows used to construct pivot row.
            let start = self.col(pivot_col_j).start;
            let stop = start + self.col(pivot_col_j).length;
            for ptr in I::range(start, stop) {
                let i = *self.ind(ptr);
                let row = self.row_mut(i);
                row.kill();
            }

            // Retrieve pivot row length.
            let pivot_row_length: I = <usize as TryInto<I>>::try_into(self.inds.length()).unwrap() - pivot_row_start;
            // Select arbitrarily pivot row.
            let pivot_row_i = if pivot_row_length > I::ZERO {
                let pos = self.col(pivot_col_j).start;
                *self.ind(pos)
            } else {
                Self::EMPTY
            };

            // Approximate degree computation.
            let start = pivot_row_start;
            let stop = pivot_row_start + pivot_row_length;
            for ptr in I::range(start, stop) {
                // Retrieve column index.
                let j = *self.ind(ptr);
                // Retrieve thickness, score and pointers.
                let weight = -self.col(j).weight;
                let score = self.col(j).rank;
                let prev = self.col(j).prev;
                let next = self.col(j).next;
                debug_assert!(weight > I::ZERO, "weight={weight}");
                debug_assert!(score >= I::ZERO);
                debug_assert!(score <= ncols);
                // Clear column.
                self.col_mut(j).weight = weight;
                // Remove column from degree list.
                if prev == Self::EMPTY {
                    *self.deg_mut(score) = next;
                } else {
                    self.col_mut(prev).next = next;
                }
                if next != Self::EMPTY {
                    self.col_mut(next).prev = prev;
                }

                // Scan the column.
                let col = self.col(j);
                let start = col.start;
                let stop = start + col.length;
                for ptr in I::range(start, stop) {
                    let i = *self.ind(ptr);
                    let row = self.row(i);
                    // Ignore dead rows.
                    if row.dead() {
                        continue;
                    }
                    debug_assert!(i != pivot_row_i);
                    // Compute set differences.
                    let mut diff = row.mark - tag;
                    // Check if row already seen.
                    if diff < I::ZERO {
                        debug_assert!(row.degree <= max_degree);
                        diff = row.degree;
                    }
                    // Substract column thickness.
                    diff -= self.col(j).weight;
                    debug_assert!(diff >= I::ZERO);
                    // Kill row if difference is zero.
                    let row = self.row_mut(i);
                    if diff == I::ZERO {
                        row.kill();
                    } else {
                        row.mark = diff + tag;
                    }
                }
            }

            // Compute score for each column.
            let start = pivot_row_start;
            let stop = pivot_row_start + pivot_row_length;
            for ptr in I::range(start, stop) {
                // Retrieve column index.
                let j = *self.ind(ptr);
                // Retrieve column.
                let col = self.col(j);
                debug_assert_ne!(j, pivot_col_j);
                debug_assert!(col.alive());
                let start = col.start;
                let stop = start + col.length;
                let mut pos = start;
                let mut hash = 0usize;
                let mut score = I::ZERO;
                for ptr in I::range(start, stop) {
                    let i = *self.ind(ptr);
                    // Ignore dead rows.
                    if self.row(i).dead() {
                        continue;
                    }
                    // Compact the column.
                    *self.ind_mut(pos) = i;
                    pos += I::ONE;
                    hash = hash.wrapping_add(i.as_usize());
                    // Increment score and prevent overflow.
                    let row = self.row(i);
                    if ncols - (row.mark - tag) < score {
                        score = ncols;
                    } else {
                        score += row.mark - tag;
                    }
                }
                // Update column length.
                let col = self.col_mut(j);
                col.length = pos - start;
                if col.length == I::ZERO {
                    col.kill();
                    pivot_row_degree -= col.weight;
                    col.rank = k;
                    k += col.weight;
                } else {
                    col.rank = score;
                    hash %= ncols.as_usize() + 1;
                    let head = self.degree[hash];
                    let first = if head > Self::EMPTY {
                        let first = self.col(head).prev;
                        self.col_mut(head).prev = j;
                        first
                    } else {
                        self.degree[hash] = -(j + I::TWO);
                        -(head + I::TWO)
                    };
                    self.col_mut(j).next = first;
                    self.col_mut(j).prev = I::try_from(hash).unwrap();
                }
            }
            // Detect super columns.
            self.detect(pivot_row_start, pivot_row_length);

            // Kill pivot column.
            self.col_mut(pivot_col_j).kill();

            // Clear marks.
            tag = self.clear(tag + max_degree + I::ONE, ncols);

            // Finalize scores.
            let start = pivot_row_start;
            let stop = pivot_row_start + pivot_row_length;
            let mut pos = start;
            for ptr in I::range(start, stop) {
                // Retrieve column index.
                let j = *self.ind(ptr);
                // Ignore dead columns.
                if self.col(j).dead() {
                    continue;
                }
                *self.ind_mut(pos) = j;
                pos += I::ONE;
                // Add pivot row to column.
                let col = self.col(j);
                let idx = col.start + col.length;
                *self.ind_mut(idx) = pivot_row_i;
                let col = self.col_mut(j);
                col.length += I::ONE;
                // Retrieve score.
                let mut score = col.rank + pivot_row_degree;
                let max = ncols - k - col.weight;
                score -= col.weight;
                score = score.min(max);
                debug_assert!(score >= I::ZERO);
                // Update score.
                col.rank = score;
                // Place column in degree list.
                let next = *self.deg(score);
                let col = self.col_mut(j);
                col.next = next;
                col.prev = Self::EMPTY;
                if next != Self::EMPTY {
                    self.col_mut(next).prev = j;
                }
                *self.deg_mut(score) = j;
                min_score = min_score.min(score);
            }

            if pivot_row_degree > I::ZERO {
                let row = self.row_mut(pivot_row_i);
                row.start = pivot_row_start;
                row.length = pos - pivot_row_start;
                row.degree = pivot_row_degree;
                row.mark = I::ZERO;
            }
        }

        Ok(())
    }

    fn order(&mut self, ncols: I) -> Result<Array<I>, ColamdError<I>> {
        for k in I::range(I::ZERO, ncols) {
            let col = self.col(k);
            if !col.dead_principal() && col.rank == Self::EMPTY {
                let mut parent = k;
                loop {
                    parent = self.col(parent).weight;
                    if self.col(parent).dead_principal() {
                        break;
                    }
                }
                let mut j = k;
                let mut order = self.col(parent).rank;
                loop {
                    self.col_mut(j).rank = order;
                    order += I::ONE;
                    self.col_mut(j).weight = parent;
                    j = self.col(j).weight;
                    if self.col(j).rank != Self::EMPTY {
                        break;
                    }
                }
                self.col_mut(parent).rank = order;
            }
        }

        let mut order = Array::new(ncols.as_usize())?;
        order.resize(ncols.as_usize(), Self::EMPTY)?;
        for j in I::range(I::ZERO, ncols) {
            let col = self.col(j);
            order[col.rank.as_usize()] = j;
        }
        Ok(order)
    }

    fn detect(&mut self, start: I, length: I) {
        for ptr in I::range(start, start + length) {
            // Retrieve column index.
            let col = *self.ind(ptr);
            if self.col(col).dead() {
                continue;
            }
            let hash = self.col(col).prev;
            let head = *self.deg(hash);
            let mut ptr = if head > Self::EMPTY { self.col(head).prev } else { -(head + I::TWO) };
            let mut prev;
            while ptr != Self::EMPTY {
                let length = self.col(ptr).length;
                let score = self.col(ptr).rank;
                prev = ptr;
                let mut next = self.col(ptr).next;
                while next != Self::EMPTY {
                    if self.col(next).length != length || self.col(next).rank != score {
                        prev = next;
                        next = self.col(next).next;
                        continue;
                    }
                    let mut p1 = self.col(ptr).start;
                    let mut p2 = self.col(next).start;
                    let mut identical = true;
                    for _ in I::range(I::ZERO, length) {
                        if self.ind(p1) != self.ind(p2) {
                            identical = false;
                            break;
                        }
                        p1 += I::ONE;
                        p2 += I::ONE;
                    }
                    if !identical {
                        prev = next;
                        next = self.col(next).next;
                        continue;
                    }
                    let weigth = self.col(next).weight;
                    self.col_mut(ptr).weight += weigth;
                    self.col_mut(next).weight = ptr;
                    self.col_mut(next).hide();
                    self.col_mut(next).rank = Self::EMPTY;
                    self.col_mut(prev).next = self.col(next).next;
                    next = self.col(next).next;
                }
                ptr = self.col(ptr).next;
            }
            if head > Self::EMPTY {
                self.col_mut(head).prev = Self::EMPTY;
            } else {
                *self.deg_mut(hash) = Self::EMPTY;
            }
        }
    }

    fn compact(&mut self, nrows: I, ncols: I) {
        // Defragment the columns.
        let mut dst = I::ZERO;
        for j in I::range(I::ZERO, ncols) {
            if self.col(j).alive() {
                let mut src = self.col(j).start;
                debug_assert!(dst <= src);
                self.col_mut(j).start = dst;
                for _ in I::range(I::ZERO, self.col(j).length) {
                    let i = *self.ind(src);
                    src += I::ONE;
                    if self.row(i).alive() {
                        self.inds[dst.as_usize()] = i;
                        dst += I::ONE;
                    }
                }
                self.col_mut(j).length = dst - self.col(j).start;
            }
        }

        // Defragment the rows.
        for i in I::range(I::ZERO, nrows) {
            let row = &mut self.rows[i.as_usize()];
            if row.dead() || row.length == I::ZERO {
                row.kill();
            } else {
                let src = row.start;
                row.mark = self.inds[src.as_usize()];
                self.inds[src.as_usize()] = -i - I::ONE;
            }
        }

        let mut src = dst;
        while src.as_usize() < self.inds.length() {
            if *self.ind(src) < I::ZERO {
                let i = -*self.ind(src) - I::ONE;
                debug_assert!(i >= I::ZERO && i <= nrows, "invalid row index `{i}`");
                *self.ind_mut(src) = self.row(i).mark;
                self.row_mut(i).start = dst;
                let length = self.row(i).length;
                for _ in I::range(I::ZERO, length) {
                    let j = *self.ind(src);
                    src += I::ONE;
                    let col = self.col(j);
                    if col.alive() {
                        *self.ind_mut(dst) = j;
                        dst += I::ONE;
                    }
                }
                self.row_mut(i).length = dst - self.row(i).start;
            } else {
                src += I::ONE;
            }
        }

        self.inds.truncate(dst.as_usize());
    }

    fn dense(control: I, size: I) -> I {
        I::from_f64((control.as_f64() * ((size.as_f64()).sqrt())).max(16.0))
    }

    fn clear(&mut self, mut tag: I, ncols: I) -> I {
        if tag <= I::ZERO || tag >= (I::MAX - ncols) {
            for row in self.rows.iter_mut() {
                if row.alive() {
                    row.mark = I::ZERO;
                }
            }
            tag = I::ONE;
        }
        tag
    }
}
