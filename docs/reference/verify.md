# Check that a claim is valid

Called on an operator, asks whether it satisfies the contract. Called on
a result together with the operator it came from, asks whether the
result is valid. Internally these run different checks; publicly they
answer one question.

## Usage

``` r
# S3 method for class 'linsolve'
verify(x, tol = 1e-08, n_probe = 5L, seed = 1L, ...)

verify(x, ...)

# S3 method for class 'linop'
verify(x, tol = 1e-08, n_probe = 10L, seed = 1L, block = c(1L, 3L), ...)
```

## Arguments

- x:

  A \`linop\`, or a solver result.

- tol:

  Relative tolerance.

- n_probe:

  Number of random probe vectors.

- seed:

  Random seed, so a report is reproducible.

- ...:

  Passed to methods.

- block:

  Block widths to test.

## Value

An object of class \`linop_certificate\`.
