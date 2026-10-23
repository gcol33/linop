# Build a preconditioner

Build a preconditioner

## Usage

``` r
preconditioner(
  apply_inverse,
  fixed = TRUE,
  hermitian = NA,
  positive_definite = NA,
  side = "left",
  dim = NULL,
  apply_inverse_adjoint = NULL
)
```

## Arguments

- apply_inverse:

  \`function(R)\` returning \`z = M^-1 r\`, always this direction.

- fixed:

  Is the map fixed across calls? A flexible preconditioner is rejected
  by solvers that require a fixed one.

- hermitian, positive_definite:

  Declared properties, three-valued.

- side:

  \`"left"\`, \`"right"\` or \`"split"\`.

- dim:

  \`c(n, n)\`.

- apply_inverse_adjoint:

  \`function(R)\` returning \`M^-H r\`. Needed only by a method that
  runs on \`A M^-1\` and has to apply its adjoint, which in v0.1 is LSQR
  and LSMR. A hermitian \`M\` supplies it for nothing, since \`M^-H =
  M^-1\` there, so this is for the case where the two differ.

## Value

An object of class \`preconditioner\`.
