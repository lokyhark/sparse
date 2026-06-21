# COLAMD algorithm Rust implementation

This crate contains a Rust implementation of the COLAMD algorithm.

## Folder architecture

`mem`: Memory harness
`ref`: Reference C implementation
`src`: Rust implementation
`valid`: Correctness harness

## References
- T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng,<br />
  *An approximate column minimum degree ordering algorithm*,<br />
  ACM Transactions on Mathematical Software, vol. 30, no. 3., pp. 353-376, 2004.<br />
  <https://doi.org/10.1145/1024074.1024079>
- T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng,<br />
  *Algorithm 836: COLAMD, an approximate column minimum degree ordering algorithm*,<br />
  ACM Transactions on Mathematical Software, vol. 30, no. 3., pp. 377-380, 2004.<br />
  <https://doi.org/10.1145/1024074.1024080>
