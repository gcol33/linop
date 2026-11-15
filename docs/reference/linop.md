# Build a linear operator

One generic, dispatching on its first argument. A matrix or \`Matrix\`
becomes an operator wrapping that storage; a function becomes a callback
operator. Third-party packages implement \`linop.myclass()\` and nothing
else.

## Usage

``` r
# S3 method for class 'DelayedMatrix'
linop(x, properties = NULL, ...)

# S3 method for class 'torch_tensor'
linop(x, properties = NULL, ...)

# S3 method for class 'big.matrix'
linop(x, properties = NULL, ...)

linop(x, ...)

# S3 method for class 'linop'
linop(x, ...)

# S3 method for class 'matrix'
linop(x, properties = NULL, ...)

# S3 method for class 'numeric'
linop(x, properties = NULL, ...)

# S3 method for class '`function`'
linop(
  x,
  adjoint = NULL,
  dim = NULL,
  dtype = "double",
  properties = NULL,
  cost = NULL,
  ...
)

# Default S3 method
linop(x, properties = NULL, ...)
```

## Arguments

- x:

  A matrix, a \`Matrix\`, a function implementing the forward action, or
  any object with a \`linop()\` method.

- properties:

  A named logical vector of declared capabilities, for example
  \`c(hermitian = TRUE, positive_definite = TRUE)\`. A bare character
  vector is accepted as shorthand for all-\`TRUE\`. Declared properties
  are recorded with \`source = "user_declaration"\`.

- ...:

  Passed to methods.

- adjoint:

  For the callback form, a function implementing the action of the
  conjugate transpose.

- dim:

  For the callback form, \`c(nrow, ncol)\`.

- dtype:

  \`"double"\` or \`"complex"\`.

- cost:

  Estimated flops per column.

## Value

A \`linop\`.

## Details

A \`DelayedMatrix\` (from the DelayedArray package, which also covers
HDF5-backed matrices from HDF5Array) becomes a callback leaf applied by
DelayedArray's block-processed product, so the matrix is never read
whole.

A two-dimensional \`torch_tensor\` becomes a callback leaf applied by
torch's matrix product on the tensor's device, with each block copied to
the device and back. Only \`float64\` and \`complex128\` tensors are
taken: a single-precision apply carries an error near 1e-7 that the
certificates' arithmetic floor, written for double precision, would not
see.

A \`big.matrix\` (bigmemory, including a file-backed, memory-mapped one)
becomes a callback leaf that reads it a block of rows at a time, so at
most one block is in memory. bigmemory stores real numbers only.

A \`linop\` behaves like a matrix: use ordinary matrix operations.

## Examples

``` r
A <- linop(diag(c(3, 2, 1)))
dim(A)
#> [1] 3 3
A %*% c(1, 1, 1)
#> [1] 3 2 1

n <- 5
S <- linop(function(X) rbind(X[-1, , drop = FALSE], 0), dim = c(n, n),
           adjoint = function(X) rbind(0, X[-n, , drop = FALSE]))
verify(S)$overall
#> [1] "qualified"
```
