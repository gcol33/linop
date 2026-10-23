# Materialise an operator as a dense matrix

An explicit request. Guarded by size, because the whole point of the
package is that this normally never happens.

## Usage

``` r
# S3 method for class 'linop'
as.matrix(x, max_entries = 1e+07, ...)
```

## Arguments

- x:

  A \`linop\`.

- max_entries:

  Guard on \`nrow \* ncol\`. Raise it deliberately.

- ...:

  Unused.

## Value

A base matrix.
