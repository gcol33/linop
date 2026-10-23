# Treat an operator P ~ A^-1 as a preconditioner

Applied by computing \`z = P r\`.

## Usage

``` r
as_preconditioner_inverse(P, ...)
```

## Arguments

- P:

  A \`linop\` approximating \`A^-1\`.

- ...:

  Passed to \[preconditioner()\].

## Value

A \`preconditioner\`.
