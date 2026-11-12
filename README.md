# Numerical methods in Pascal

Two interactive Pascal programs demonstrate elementary operations on 3 × 3 matrices:

- `gauss.pas` solves a linear system using Gaussian elimination with partial pivoting.
- `inversemat.pas` computes a matrix inverse using the determinant and adjugate.

## Build and run

The source is written in Pascal. With a compatible Free Pascal installation, try compiling each program separately:

```sh
fpc gauss.pas
fpc inversemat.pas
./gauss
./inversemat
```

Enter matrix elements when prompted; the linear-system program also asks for the right-hand-side vector. These are legacy educational examples with fixed-size arrays, not general-purpose numerical libraries. Compiler compatibility and behavior for singular or poorly conditioned matrices should be checked before use.
