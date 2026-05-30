use sparse;

#[test]
fn colamd() {
    // Matrix  5-by-4
    // x 0 x 0
    // x 0 x x
    // 0 x x 0
    // 0 0 x x
    // x x 0 0
    let rowind = [0, 1, 4, 2, 4, 0, 1, 2, 3, 1, 3];
    let colptr = [0, 3, 5, 9, 11];
    let result = sparse::colamd(5, 4, &colptr, &rowind).unwrap();
    let order = result.order();
    assert_eq!(order[0], 1);
    assert_eq!(order[1], 0);
    assert_eq!(order[2], 2);
    assert_eq!(order[3], 3);
}
