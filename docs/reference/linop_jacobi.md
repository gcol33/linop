# A Jacobi operator on the sequence space

A self-adjoint operator on `l^2(Z)`, tridiagonal in the standard basis,
given by its diagonal and off-diagonal sequences. Both are real and
eventually constant: they take the values you supply on a finite window
and their limiting values everywhere else.

## Usage

``` r
linop_jacobi(
  diagonal,
  offdiagonal = 1,
  from = NULL,
  diagonal_limit = 0,
  offdiagonal_limit = NULL
)
```

## Arguments

- diagonal:

  The diagonal entries on the window, as a real vector. A discrete
  Schrodinger operator is this argument alone: it is the potential.

- offdiagonal:

  The couplings, as a real vector of length 1 for a constant
  off-diagonal, or `length(diagonal) + 1` for a perturbed one. Zero is
  not allowed as the limiting value.

- from:

  The index of `diagonal[1]`. The default centres the window on the
  origin.

- diagonal_limit:

  The value the diagonal takes outside the window, or one period of
  values: outside the window `a_j = diagonal_limit[(j mod P) + 1]`. With
  periodic limits the essential spectrum is a union of bands with gaps
  between them, and the certificate measures isolation against the
  nearest band.

- offdiagonal_limit:

  The value the off-diagonal takes outside the window, or one period of
  values, anchored the same way. Defaults to `offdiagonal` when that is
  a single number.

## Value

A `linop` of dimension `Inf x Inf`.

## Details

The operator has infinite dimensions, so it has no matrix and nothing
numeric runs on it directly.
[`finite_section()`](https://gillescolling.com/linop/reference/finite_section.md)
truncates it to an operator that does, and the certificate that comes
back describes the operator here rather than the truncation.

The window is `from` to `from + length(diagonal) - 1`. `offdiagonal[i]`
couples index `from + i - 2` to `from + i - 1`, so the couplings
incident to the diagonal window run from the one on its left edge to the
one on its right, and there are one more of them than there are diagonal
entries.

A constructor knows the operator is self-adjoint the way
[`linop_eye()`](https://gillescolling.com/linop/reference/linop_eye.md)
knows it about the identity, so `hermitian` and `symmetric` are recorded
with `source = "construction"` rather than as declarations.

## Examples

``` r
## the potential v at the origin, whose top eigenvalue is sign(v) sqrt(v^2 + 4)
H <- linop_jacobi(diagonal = 1)
dim(H)
#> [1] Inf Inf
F <- finite_section(H, n = 50)
fit <- eigs(F, k = 1, which = "largest_algebraic")
fit$values - sqrt(5)
#> [1] -4.440892e-16
fit$certificate
#> solved by primme (hermitian operator, PRIMME with implicit restarting (highest priority))
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         construction     deterministic_bound   1
#> finite residual                pass         computation      identity              1
#> orthogonality                  pass         computation      identity              1
#> truncation bound               pass         theorem          deterministic_bound   1
#> isolation                      pass         theorem          deterministic_bound   1
#> enclosure                      pass         theorem          estimate              -
#> convergence                    pass         computation      identity              1
#> forward error                  pass         theorem          deterministic_bound   1
#> ------------------------------------------------------------------------------------------
#> overall                        pass   no deterministic bound on: enclosure
```
