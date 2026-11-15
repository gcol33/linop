# Permutation operator

The operator taking \`x\` to \`x\[p\]\`. Its transpose is its inverse,
so it is unitary by construction, and it is symmetric exactly when \`p\`
is an involution.

## Usage

``` r
linop_perm(p)
```

## Arguments

- p:

  A permutation of \`1:n\`.

## Value

A \`linop\`.

## Examples

``` r
P <- linop_perm(c(3, 1, 2))
P %*% c(10, 20, 30)
#> [1] 30 10 20
```
