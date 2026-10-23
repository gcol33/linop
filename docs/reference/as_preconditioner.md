# Treat an approximate inverse as a preconditioner

A \`linop\` \`M ~ A\` is applied by solving \`M z = r\`, so it needs a
solver; the inverse is never formed here. A \`linsolve\` already applies
\`A^-1\`, so it needs nothing further, and its contract (section 4.2)
determines whether the resulting preconditioner is fixed.

## Usage

``` r
as_preconditioner(x, ...)

# S3 method for class 'linop'
as_preconditioner(x, solver, ...)

# S3 method for class 'linsolve'
as_preconditioner(x, side = NULL, ...)

# Default S3 method
as_preconditioner(x, ...)
```

## Arguments

- x:

  A \`linop\` approximating \`A\`, or a \`linsolve\`.

- ...:

  Passed to \[preconditioner()\].

- solver:

  \`function(M, R)\` performing the solve. \`linop\` method only.

- side:

  \`"left"\`, \`"right"\` or \`"split"\`. Defaults by contract; see
  details.

## Value

A \`preconditioner\`.
