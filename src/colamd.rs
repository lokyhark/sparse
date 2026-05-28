use crate::array::Array;

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
use error::{ColamdError, ColamdErrorKind};
use int::ColamdInt;
use result::ColamdResult;
use row::ColamdRow;
use score::ColamdScore;
use size::ColamdSize;

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
/// - T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng,
///   *An approximate column minimum degree ordering algorithm*,
///   ACM Transactions on Mathematical Software, vol. 30, no. 3., pp. 353-376, 2004.
///   <https://doi.org/10.1145/1024074.1024079>
/// - T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng,
///   *Algorithm 836: COLAMD, an approximate column minimum degree ordering algorithm*,
///   ACM Transactions on Mathematical Software, vol. 30, no. 3., pp. 377-380, 2004.
///   <https://doi.org/10.1145/1024074.1024080>
pub fn colamd<I: ColamdInt>(nrows: I, ncols: I, colptr: &[I], rowind: &[I]) -> Result<ColamdResult, ColamdError<I>> {
    // Check inputs.
    let size = Colamd::check(nrows, ncols, colptr, rowind)?;
    // Allocate algorithm.
    let elbow = size.ncols.as_usize().saturating_add(size.nnz.as_usize() / 5);
    let mut colamd: Colamd<I> = Colamd::alloc(size.nrows.as_usize(), size.ncols.as_usize(), size.nnz.as_usize(), elbow)?;
    // Execute algorithm.
    colamd.exec(size.nrows, size.ncols, size.nnz, colptr, rowind)
}

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
    pub fn with_config(mut self, config: ColamdConfig<I>) -> Self {
        self.config = config;
        self
    }

    /// Check COLAMD algorithm input arguments and return the required workspace size.
    ///
    /// Following proporties are checked:
    /// - 0 <= nrows < isize::MAX + 1
    /// - 0 <= ncols < isize::MAX + 1
    /// - colptr.len() = ncols + 1
    /// - colptr[0] = 0
    /// - colptr[ncols + 1] = nnz
    /// - 0 <= nnz < isize::MAX + 1
    /// - rowind.len() = nnz
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
        debug_assert!(MIN <= nrows.as_usize() && nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}",);
        debug_assert!(MIN <= ncols.as_usize() && ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");

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
        debug_assert!(MIN <= nrows.as_usize() && nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}",);
        debug_assert!(MIN <= ncols.as_usize() && ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert_eq!(colptr.len(), ncols.as_usize() + 1, "invalid col pointers slice length: actual={}, expect={}", colptr.len(), ncols.as_usize() + 1);

        // SAFETY: colptr is not empty (see previous checkpoint).
        if colptr[0] != I::ZERO {
            return Err(ColamdErrorKind::InvalidFirstColPtr { ptr: colptr[0] }.into());
        }

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
        debug_assert!(MIN <= nrows.as_usize() && nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}",);
        debug_assert!(MIN <= ncols.as_usize() && ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
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
        debug_assert!(MIN <= nrows.as_usize() && nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}",);
        debug_assert!(MIN <= ncols.as_usize() && ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
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
    pub fn run(&mut self, nrows: I, ncols: I, colptr: &[I], rowind: &[I]) -> Result<ColamdResult, ColamdError<I>> {
        // Check input.
        let size = Colamd::check(nrows, ncols, colptr, rowind)?;
        // Execute algorithm.
        self.exec(size.nrows, size.ncols, size.nnz, colptr, rowind)
    }

    /// Execute the COLAMD algorithm.
    ///
    /// On input following invariants must hold:
    /// - 0 <= nrows < isize::MAX + 1
    /// - 0 <= ncols < isize::MAX + 1
    /// - colptr.len() = ncols + 1
    /// - colptr[0] = 0
    /// - colptr[ncols + 1] = nnz
    /// - 0 <= nnz < isize::MAX + 1
    /// - rowind.len() = nnz
    fn exec(&mut self, nrows: I, ncols: I, nnz: I, colptr: &[I], rowind: &[I]) -> Result<ColamdResult, ColamdError<I>> {
        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // colptr.len() = ncols + 1
        // colptr[0] = 0
        // colptr[ncols + 1] = nnz
        // 0 <= nnz < isize::MAX + 1
        // rowind.len() = nnz
        debug_assert!(MIN <= nrows.as_usize() && nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}",);
        debug_assert!(MIN <= ncols.as_usize() && ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
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
        self.find(nrows, ncols, nnz, score.cols, score.min_score, score.max_degree)?;

        // Order.
        self.order(ncols);

        Ok(ColamdResult {})
    }

    /// Initialize the COLAMD algorithm.
    ///
    /// On input following invariants must hold:
    /// - 0 <= nrows < isize::MAX + 1
    /// - 0 <= ncols < isize::MAX + 1
    /// - colptr.len() = ncols + 1
    /// - colptr[0] = 0
    /// - colptr[ncols + 1] = nnz
    /// - 0 <= nnz < isize::MAX + 1
    /// - rowind.len() = nnz
    ///
    /// On output following invariants hold:
    /// - 0 <= nrows < isize::MAX + 1
    /// - 0 <= ncols < isize::MAX + 1
    /// - colptr.len() = ncols + 1
    /// - colptr[0] = 0
    /// - colptr[ncols + 1] = nnz
    /// - 0 <= nnz < isize::MAX + 1
    /// - rowind.len() = nnz
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
        debug_assert!(MIN <= nrows.as_usize() && nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}",);
        debug_assert!(MIN <= ncols.as_usize() && ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
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

        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // colptr.len() = ncols + 1
        // colptr[0] = 0
        // colptr[ncols + 1] = nnz
        // 0 <= nnz < isize::MAX + 1
        // rowind.len() = nnz
        // colptr is ordered
        // colptr pointers are valid indices
        debug_assert!(MIN <= nrows.as_usize() && nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}",);
        debug_assert!(MIN <= ncols.as_usize() && ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert_eq!(colptr.len(), ncols.as_usize() + 1, "invalid col pointers slice length: actual={}, expect={}", colptr.len(), ncols.as_usize() + 1);
        debug_assert_eq!(colptr[0], I::ZERO, "invalid first column pointer: actual={}, expect=0", colptr[0]);
        debug_assert_eq!(colptr[ncols.as_usize()], nnz, "invalid number of nnz: actual={}, expect={nnz}", colptr[ncols.as_usize()]);
        debug_assert!(MIN < nnz.as_usize() && nnz.as_usize() <= MAX, "invalid number of non zeros: nnz={nnz}, min={MIN}, max={MAX}");
        debug_assert_eq!(rowind.len(), nnz.as_usize(), "invalid row indices slice length: actual={}, expect={nnz}", rowind.len());
        debug_assert!(colptr.is_sorted(), "unsorted column pointers");
        debug_assert!(colptr.iter().all(|&i| I::ZERO <= i && i <= nnz), "invalid column pointers");

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
            let col = &self.cols[j.as_usize()];
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
                } else {
                    prev = row;
                }
                // SAFETY: row is valid index because >= 0 and < nrows which is a valid index.
                let idx = row.as_usize();
                let row = &mut self.rows[idx];
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
        debug_assert!(MIN <= nrows.as_usize() && nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}",);
        debug_assert!(MIN <= ncols.as_usize() && ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
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
            let col = &self.cols[j.as_usize()];
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
        debug_assert!(MIN <= nrows.as_usize() && nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}",);
        debug_assert!(MIN <= ncols.as_usize() && ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert!(MIN < nnz.as_usize() && nnz.as_usize() <= MAX, "invalid number of non zeros: nnz={nnz}, min={MIN}, max={MAX}");

        // Compute dense row and column thresholds.
        // SAFETY:
        // - 16 is representable by all unsigned integer types (ColamdInt which is sealed).
        // - ncols and nrows are representable by I (checked previously).
        // - dense_row_control and dense_col_control are non negative and less than nrows and ncols which are representable by I.
        let dense_row_count = self.config.dense_row_control * ncols.isqrt().max(unsafe { 16.try_into().unwrap_unchecked() });
        let dense_col_count = self.config.dense_col_control * nrows.min(ncols).isqrt().max(unsafe { 16.try_into().unwrap_unchecked() });

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
        for col in self.cols.iter_mut().rev() {
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
                    let i = self.inds[ptr.as_usize()];
                    // SAFETY: row indices are valid (checked previously).
                    let row = &mut self.rows[i.as_usize()];
                    // Decrement row degree.
                    row.degree -= I::ONE;
                }
                // Kill column and place it last.
                cols -= I::ONE;
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
        for col in self.cols.iter_mut().rev() {
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
                let i = self.inds[ptr.as_usize()];
                // SAFETY: row indices are valid (checked previously).
                let row = &self.rows[i.as_usize()];
                // Ignore dead rows.
                if row.dead() {
                    continue;
                }
                // Compact the column.
                self.inds[pos.as_usize()] = i;
                // Update pointer.
                pos += I::ONE;
                // Prevent overflow by saturating at ncols which is valid.
                debug_assert!(row.degree > I::ZERO, "invalid alive row degree: degree={}", row.degree);
                if ncols - (row.degree - I::ONE) < score {
                    score = ncols;
                } else {
                    score += row.degree - I::ONE;
                }
            }
            // Compute pruned column length.
            let length = pos - start;
            // Kill column and place it last if it is now empty.
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
            let col = &mut self.cols[j.as_usize()];
            // Ignore dead column.
            if col.dead() {
                continue;
            }
            let score = col.rank;
            debug_assert!(min_score >= I::ZERO && min_score <= ncols, "invalid minimum column score: score={}, min={}, max={}", score, I::ZERO, ncols);
            debug_assert!(score >= I::ZERO && score <= ncols, "invalid column score: score={}, min={}, max={}", score, I::ZERO, ncols);
            let head = &mut self.degree[score.as_usize()];
            col.prev = Self::EMPTY;
            col.next = *head;
            if *head != Self::EMPTY {
                self.cols[head.as_usize()].prev = j;
            }
            *head = j;
            // Keep track of minimum score.
            min_score = min_score.min(score);
        }

        Ok(ColamdScore { cols, rows, max_degree, min_score })
    }

    pub fn find(&mut self, nrows: I, ncols: I, nnz: I, cols: I, mut min_score: I, mut max_degree: I) -> Result<(), ColamdError<I>> {
        // CHECKPOINT:
        // 0 <= nrows < isize::MAX + 1
        // 0 <= ncols < isize::MAX + 1
        // 0 <= nnz < isize::MAX + 1
        // cols <= ncols
        // 0 < min_score <= ncols
        debug_assert!(MIN <= nrows.as_usize() && nrows.as_usize() <= MAX, "invalid number of rows: nrows={nrows}, min={MIN}, max={MAX}",);
        debug_assert!(MIN <= ncols.as_usize() && ncols.as_usize() <= MAX, "invalid number of columns: ncols={ncols}, min={MIN}, max={MAX}");
        debug_assert!(MIN < nnz.as_usize() && nnz.as_usize() <= MAX, "invalid number of non zeros: nnz={nnz}, min={MIN}, max={MAX}");
        debug_assert!(I::ZERO < cols && cols <= ncols, "invalid number of alive columns: cols={}, min={}, max={}", cols, I::ZERO, ncols);
        debug_assert!(min_score > I::ZERO && min_score <= ncols, "invalid minimum score");

        let mut tag = I::ONE;
        let mut score = min_score;
        let mut k = I::ZERO;
        while k < ncols {
            // Choose pivot column with minimum score.
            let pivot_col_j = loop {
                let head = self.degree[score.as_usize()];
                if head != Self::EMPTY {
                    break head;
                } else {
                    score += I::ONE;
                }
            };
            debug_assert!(pivot_col_j >= I::ZERO && pivot_col_j <= ncols);
            // Remove pivot column from degree list by placing next column as head.
            let next = self.cols[pivot_col_j.as_usize()].next;
            self.degree[score.as_usize()] = next;
            if next != Self::EMPTY {
                self.cols[next.as_usize()].prev = Self::EMPTY;
            }
            // Retrieve pivot column.
            let pivot_col = &mut self.cols[pivot_col_j.as_usize()];
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
                self.compact(nrows, ncols)
            }

            // Pivot row start point.
            let pivot_row_start: I = self.inds.length().try_into().unwrap();

            // Initialize pivot row degree.
            let mut pivot_row_degree = I::ZERO;

            // Mark pivot column.
            self.cols[pivot_col_j.as_usize()].weight = -pivot_col_thickness;

            // Create the pivot row which is the union of all rows in the pivot column.
            let start = self.cols[pivot_col_j.as_usize()].start;
            let stop = start + self.cols[pivot_col_j.as_usize()].length;
            for ptr in I::range(start, stop) {
                let i = self.inds[ptr.as_usize()];
                let row = &mut self.rows[i.as_usize()];
                // Ignore dead rows.
                if row.alive() {
                    let start = row.start;
                    let stop = start + row.length;
                    for ptr in I::range(start, stop) {
                        let j = self.inds[ptr.as_usize()];
                        let col = &mut self.cols[j.as_usize()];
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
            self.cols[pivot_col_j.as_usize()].weight = pivot_col_thickness;
            max_degree = max_degree.max(pivot_row_degree);

            // Kill rows used to construct pivot row.
            let start = self.cols[pivot_col_j.as_usize()].start;
            let stop = start + self.cols[pivot_col_j.as_usize()].length;
            for ptr in I::range(start, stop) {
                let i = self.inds[ptr.as_usize()];
                let row = &mut self.rows[i.as_usize()];
                row.kill();
            }

            // Retrieve pivot row length.
            let pivot_row_length: I = <usize as TryInto<I>>::try_into(self.inds.length()).unwrap() - pivot_row_start;
            // Select arbitrarily pivot row.
            let pivot_row_i = if pivot_row_length > I::ZERO {
                let pos = self.cols[pivot_col_j.as_usize()].start;
                self.inds[pos.as_usize()]
            } else {
                Self::EMPTY
            };

            // Approximate degree computation.
            let start = pivot_row_start;
            let stop = pivot_row_start + pivot_row_length;
            for ptr in I::range(start, stop) {
                // Retrieve column index.
                let j = self.inds[ptr.as_usize()];
                // Retrieve thickness, score and pointers.
                let weight = -self.cols[j.as_usize()].weight;
                let score = self.cols[j.as_usize()].rank;
                let prev = self.cols[j.as_usize()].prev;
                let next = self.cols[j.as_usize()].next;
                debug_assert!(weight > I::ZERO, "weight={}", weight);
                debug_assert!(score >= I::ZERO);
                debug_assert!(score <= ncols);
                // Clear column.
                self.cols[j.as_usize()].weight = weight;
                // Remove column from degree list.
                if prev == Self::EMPTY {
                    self.degree[score.as_usize()] = next;
                } else {
                    self.cols[prev.as_usize()].next = next;
                }
                if next != Self::EMPTY {
                    self.cols[next.as_usize()].prev = prev;
                }

                // Scan the column.
                let col = &self.cols[j.as_usize()];
                let start = col.start;
                let stop = start + col.length;
                for ptr in I::range(start, stop) {
                    let i = self.inds[ptr.as_usize()];
                    let row = &mut self.rows[i.as_usize()];
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
                    diff -= col.weight;
                    debug_assert!(diff >= I::ZERO);
                    // Kill row if difference is zero.
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
                let j = self.inds[ptr.as_usize()];
                // Retrieve column.
                let col = &mut self.cols[j.as_usize()];
                debug_assert_ne!(j, pivot_col_j);
                debug_assert!(col.alive());
                let start = col.start;
                let stop = start + col.length;
                let mut pos = start;
                let mut hash = I::ZERO;
                let mut score = I::ZERO;
                for ptr in I::range(start, stop) {
                    let i = self.inds[ptr.as_usize()];
                    let row = &self.rows[i.as_usize()];
                    // Ignore dead rows.
                    if row.dead() {
                        continue;
                    }
                    // Compact the column.
                    self.inds[pos.as_usize()] = i;
                    pos += I::ONE;
                    hash += i;
                    // Increment score and prevent overflow.
                    if ncols - (row.mark - tag) < score {
                        score = ncols;
                    } else {
                        score += row.mark - tag;
                    }
                }
                // Update column length.
                col.length = pos - start;
                if col.length == I::ZERO {
                    col.kill();
                    pivot_row_degree -= col.weight;
                    col.rank = k;
                    k += col.weight;
                } else {
                    col.rank = score;
                    hash %= ncols + I::ONE;
                    let head = self.degree[hash.as_usize()];
                    let first = if head != Self::EMPTY {
                        std::mem::replace(&mut col.prev, j)
                    } else {
                        self.degree[hash.as_usize()] = -(j + I::TWO);
                        -(head + I::TWO)
                    };
                    col.next = first;
                    col.prev = hash;
                }
            }

            // Detect super columns.
            self.detect(ncols, nnz, pivot_row_start, pivot_row_length);

            // Kill pivot column.
            self.cols[pivot_col_j.as_usize()].kill();

            // Clear marks.
            tag = self.clear(tag + max_degree + I::ONE, ncols);

            // Finalize scores.
            let start = pivot_row_start;
            let stop = pivot_row_start + pivot_row_length;
            let mut pos = start;
            for ptr in I::range(start, stop) {
                // Retrieve column index.
                let j = self.inds[ptr.as_usize()];
                // Retrieve column.
                let col = &mut self.cols[j.as_usize()];
                // Ignore dead columns.
                if col.dead() {
                    continue;
                }
                self.inds[pos.as_usize()] = j;
                // Add pivot row to column.
                let idx = col.start + col.length;
                self.inds[idx.as_usize()] = pivot_row_i;
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
                let next = self.degree[score.as_usize()];
                col.next = next;
                col.prev = Self::EMPTY;
                if next != Self::EMPTY {
                    self.cols[next.as_usize()].prev = j;
                }
                self.degree[score.as_usize()] = j;
                min_score = min_score.min(score);
            }

            if pivot_row_degree > I::ZERO {
                let row = &mut self.rows[pivot_row_i.as_usize()];
                row.start = pivot_row_start;
                row.length = pos - pivot_row_start;
                row.degree = pivot_row_degree;
                row.mark = I::ZERO;
            }
            for (r, row) in self.rows.iter().enumerate() {
                println!("row = {}, start = {}, length = {}, degree = {}, mark = {}", r, row.start, row.length, row.degree, row.mark);
            }
            for (c, col) in self.cols.iter().enumerate() {
                println!(
                    "col = {}, start = {}, length = {}, weight = {}, rank = {}, prev = {}, next = {}",
                    c, col.start, col.length, col.weight, col.rank, col.prev, col.next
                );
            }
        }

        Ok(())
    }

    fn order(&mut self, ncols: I) {
        for i in I::range(I::ZERO, ncols) {
            let mut parent = i;
            loop {
                parent = self.cols[parent.as_usize()].weight;
                break;
            }
            let mut order = self.cols[parent.as_usize()].rank;
            let mut col = i;
            loop {
                self.cols[col.as_usize()].rank = order;
                order += I::ONE;
                self.cols[col.as_usize()].weight = parent;
                col = self.cols[col.as_usize()].weight;
                if col == parent {
                    break;
                }
            }
            self.cols[parent.as_usize()].rank = order;
        }
    }

    fn detect(&mut self, ncols: I, nnz: I, start: I, length: I) {
        for ptr in I::range(start, start + length) {
            // Retrieve column index.
            let col = self.inds[ptr.as_usize()];
            if self.cols[col.as_usize()].dead() {
                continue;
            }
            let hash = self.cols[col.as_usize()].prev;
            let head = self.degree[hash.as_usize()];
            let mut ptr = if head >= I::ZERO { self.cols[head.as_usize()].prev } else { -(head + I::TWO) };
            let mut prev;
            while ptr != Self::EMPTY {
                let length = self.cols[ptr.as_usize()].length;
                let score = self.cols[ptr.as_usize()].rank;
                prev = ptr;
                let mut next = self.cols[col.as_usize()].next;
                while next != Self::EMPTY {
                    if self.cols[next.as_usize()].length != length || self.cols[next.as_usize()].rank != score {
                        prev = next;
                        continue;
                    }
                    let p1 = self.cols[ptr.as_usize()].start;
                    let p2 = self.cols[next.as_usize()].start;
                    let mut identical = true;
                    for _ in I::range(I::ZERO, length) {
                        if self.inds[p1.as_usize()] != self.inds[p2.as_usize()] {
                            identical = false;
                            break;
                        }
                    }
                    if !identical {
                        prev = next;
                        continue;
                    }
                    let weigth = self.cols[next.as_usize()].weight;
                    self.cols[ptr.as_usize()].weight += weigth;
                    self.cols[next.as_usize()].weight = ptr;
                    self.cols[next.as_usize()].kill();
                    self.cols[next.as_usize()].rank = Self::EMPTY;
                    self.cols[prev.as_usize()].next = self.cols[next.as_usize()].next;
                    next = self.cols[next.as_usize()].next;
                }
                ptr = self.cols[ptr.as_usize()].next;
            }
            if head > Self::EMPTY {
                self.cols[head.as_usize()].prev = Self::EMPTY;
            } else {
                self.degree[hash.as_usize()] = Self::EMPTY;
            }
        }
    }

    fn compact(&mut self, nrows: I, ncols: I) {
        // Defragment the columns.
        let mut pos = I::ZERO;
        for j in I::range(I::ZERO, ncols) {
            let col = &mut self.cols[j.as_usize()];
            if col.alive() {
                let mut ptr = col.start;
                col.start = pos;
                for _ in I::range(I::ZERO, col.length) {
                    let i = self.inds[ptr.as_usize()];
                    ptr += I::ONE;
                    let row = &mut self.rows[i.as_usize()];
                    if row.alive() {
                        self.inds[pos.as_usize()] = i;
                        pos += I::ONE;
                    }
                }
                col.length = pos - col.start;
            }
        }

        // Defragment the rows.
        for i in I::range(I::ZERO, nrows) {
            let row = &mut self.rows[i.as_usize()];
            if row.dead() || row.length == I::ZERO {
                row.kill();
            } else {
                let ptr = row.start;
                row.mark = ptr;
                self.inds[ptr.as_usize()] = -ptr - I::ONE;
            }
        }
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
