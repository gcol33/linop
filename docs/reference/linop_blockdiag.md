# Block diagonal operator

The operator acting on each block of its input with its own operand,
blocks in the order given. The blocks may be rectangular. Every
capability the blocks share carries over, and full rank does when the
blocks are all tall or all wide.

## Usage

``` r
linop_blockdiag(...)
```

## Arguments

- ...:

  \`linop\`s or matrices, at least one.

## Value

A \`linop\`.

## Examples

``` r
D <- linop_blockdiag(matrix(c(2, 1, 1, 2), 2), linop_scaling(c(3, 4, 5)))
as.matrix(D)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    2    1    0    0    0
#> [2,]    1    2    0    0    0
#> [3,]    0    0    3    0    0
#> [4,]    0    0    0    4    0
#> [5,]    0    0    0    0    5
```
