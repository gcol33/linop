# Low-rank operator

\`U diag(d) V^H\`, applied as three thin products, \`(m + n) r\` work
per column for rank \`r\`. With \`V\` left at \`U\` and \`d\` real, it
is hermitian by construction, and positive semidefinite when \`d\` is
also non-negative. A rank of zero is the zero operator.

## Usage

``` r
linop_lowrank(U, V = U, d = rep(1, ncol(U)))
```

## Arguments

- U:

  An \`m x r\` matrix.

- V:

  An \`n x r\` matrix, by default \`U\`.

- d:

  The \`r\` weights, by default all one.

## Value

A \`linop\`.

## Examples

``` r
U <- matrix(rnorm(20), 10, 2)
P <- linop_lowrank(U, d = c(3, 1))
max(abs(as.matrix(P) - U %*% diag(c(3, 1)) %*% t(U)))
#> [1] 3.552714e-15
```
