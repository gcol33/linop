# Kronecker product of two operators

\`A (x) B\`, in the ordering of base R's \[kronecker()\]: the block in
row block \`i\` and column block \`j\` is \`A\[i, j\] \* B\`. It is
applied without forming either factor's entries, as \`vec(B X A^T)\`, at
one apply of each factor per block of right-hand sides. Hermitian,
symmetric, definite, unitary, diagonal, triangular and real carry over
when both factors have them.

## Usage

``` r
linop_kron(A, B)
```

## Arguments

- A, B:

  \`linop\`s or matrices.

## Value

A \`linop\` of dimension \`nrow(A) nrow(B)\` by \`ncol(A) ncol(B)\`.

## Examples

``` r
A <- matrix(c(2, 1, 0, 3), 2)
B <- matrix(1:6, 3)
K <- linop_kron(A, B)
max(abs(as.matrix(K) - kronecker(A, B)))
#> [1] 0
## the 2-D Laplacian on a 4 x 4 grid, from the 1-D one
L1 <- diag(2, 4)
L1[cbind(1:3, 2:4)] <- -1
L1[cbind(2:4, 1:3)] <- -1
L2 <- linop_kron(L1, linop_eye(4)) + linop_kron(linop_eye(4), linop(L1))
dim(L2)
#> [1] 16 16
```
