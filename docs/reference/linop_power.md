# Integer power of an operator

`A^k`, applied as `k` applies of `A`. `k = 0` is the identity; a
negative `k` is a power of the exact inverse, so it needs an operator
`solve(A)` can invert. A hermitian operator to an even power is positive
semidefinite, and positive definite when it is also established of full
rank.

## Usage

``` r
linop_power(A, k)
```

## Arguments

- A:

  A square `linop` or matrix.

- k:

  A whole number.

## Value

A `linop`.

## Examples

``` r
A <- linop(matrix(c(2, 1, 1, 3), 2))
as.matrix(linop_power(A, 3))
#>      [,1] [,2]
#> [1,]   15   20
#> [2,]   20   35
```
