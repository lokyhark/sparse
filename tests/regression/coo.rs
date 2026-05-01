use sparse::CoordinateMatrix;

/// Test for detecting double free and aliasing issues when cloning pattern.
#[test]
fn clone_pattern() {
    let pat = sparse::CooPat::new(2, 3, 10);
    let clone = pat.clone();
    assert_eq!(pat.nrows(), clone.nrows());
    assert_eq!(pat.ncols(), clone.ncols());
    assert_eq!(pat.length(), clone.length());
    assert_eq!(pat.capacity(), clone.capacity());
}

/// Test for detecting double free and aliasing issues when cloning matrix.
#[test]
fn clone_matrix() {
    let mat: CoordinateMatrix<i32, f64> = sparse::CooMat::new(2, 3, 10);
    let clone = mat.clone();
    assert_eq!(mat.nrows(), clone.nrows());
    assert_eq!(mat.ncols(), clone.ncols());
    assert_eq!(mat.length(), clone.length());
    assert_eq!(mat.capacity(), clone.capacity());
}
