# Build a linear operator

One generic, dispatching on its first argument. A matrix or \`Matrix\`
becomes an operator wrapping that storage; a function becomes a callback
operator. Third-party packages implement \`linop.myclass()\` and nothing
else.

## Usage

``` r
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

- ...:

  Passed to methods.

- properties:

  A named logical vector of declared capabilities, for example
  \`c(hermitian = TRUE, positive_definite = TRUE)\`. A bare character
  vector is accepted as shorthand for all-\`TRUE\`. Declared properties
  are recorded with \`source = "user_declaration"\`.

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

A \`linop\` behaves like a matrix: use ordinary matrix operations.

## Examples

``` r
A <- linop(diag(c(3, 2, 1)))
dim(A)
A %*% c(1, 1, 1)

n <- 5
S <- linop(function(X) rbind(X[-1, , drop = FALSE], 0), dim = c(n, n),
           adjoint = function(X) rbind(0, X[-n, , drop = FALSE]))
verify(S)$overall
```
