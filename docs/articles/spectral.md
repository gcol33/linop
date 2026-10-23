# Eigenvalues and singular values

## Installation

``` r

install.packages("pak")
pak::pak("gcol33/linop")
```

Installing `linop` installs one package. There are no `Imports` and no
compiled code, so there is no toolchain to have working first. `Matrix`
sits in `Suggests` and is used for sparse leaves and
[`as_sparse()`](https://gillescolling.com/linop/reference/as_sparse.md).

R 4.4.0 or newer. The floor is set by
[`crossprod()`](https://rdrr.io/r/base/crossprod.html) and
[`tcrossprod()`](https://rdrr.io/r/base/crossprod.html), which became S3
generic in R 4.4.0, and are what let one `matrixOps` method cover all
three matrix products.

## Two verbs

[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) computes a
few eigenpairs of a hermitian operator.
[`svds()`](https://gillescolling.com/linop/reference/svds.md) computes a
few singular triplets of any operator. Both run against the action
alone, and both hold `ncv` vectors rather than a matrix.

The 1D Laplacian has closed-form eigenvalues, which makes it the fixture
to check against:

``` r

n <- 200
lap <- function(X) 2 * X - rbind(X[-1, , drop = FALSE], 0) - rbind(0, X[-n, , drop = FALSE])
L <- linop(lap, adjoint = lap, dim = c(n, n), properties = c(hermitian = TRUE))

truth <- 2 - 2 * cos(seq_len(n) * pi / (n + 1))

fit <- eigs(L, k = 3, which = "largest_algebraic", ncv = 60)
fit
#> <primme> 0 of 3 requested eigenpairs converged, which = 'largest_algebraic'
#>   300 iterations of at most 300, 9 restarts
#>   values: 3.999756  3.999023  3.997802
#>   certificate: fail
```

``` r

rbind(computed = fit$values,
      closed_form = rev(sort(truth))[1:3])
#>                 [,1]     [,2]     [,3]
#> computed    3.999756 3.999023 3.997802
#> closed_form 3.999756 3.999023 3.997802
```

## The operator has to say it is hermitian

There is one eigensolver here and it is Lanczos, which is defined for a
hermitian operator. Symmetry is not probed and not assumed:

``` r

Undeclared <- linop(lap, adjoint = lap, dim = c(n, n))
eigs(Undeclared, k = 2)
#> Error:
#> ! method 'eigs' requires the operator to be hermitian; it declares hermitian = NA.
#>   Unknown is not false. Declare it with linop(..., properties = c(hermitian = TRUE)) if you can justify it.
```

Unknown is not false, and the refusal says so rather than guessing
either way. A dense leaf establishes symmetry by checking it, so a
matrix needs no declaration; a callback has nothing to check and the
declaration is how you supply it.

A non-hermitian operator is refused for a different reason, and there is
no Arnoldi in this version, so it has no eigensolver here at all.
[`svds()`](https://gillescolling.com/linop/reference/svds.md) works on
any operator and is often the thing that was actually wanted.

## `ncv` is the knob that matters

`ncv` is how many vectors a restart holds. Each restart compresses the
subspace to the vectors it keeps, so the smaller `ncv` is, the more
applies the method spends reaching the tolerance:

``` r

ladder <- t(vapply(c(20, 30, 40, 60), function(nc) {
  f <- eigs(L, k = 3, which = "largest_algebraic", ncv = nc)
  c(ncv = nc,
    error = max(abs(f$values - rev(sort(truth))[1:3])),
    converged = f$nconv,
    applies = f$iterations)
}, numeric(4)))
ladder
#>      ncv        error converged applies
#> [1,]  20 2.040628e-09         0     302
#> [2,]  30 4.101608e-12         0     300
#> [3,]  40 2.620126e-14         0     301
#> [4,]  60 1.332268e-15         0     300
```

At `ncv = 20` the default budget of 300 applies runs out before the
residuals reach the tolerance, and the run reports the pairs it has as
unconverged. A larger budget recovers it:

``` r

eigs(L, k = 3, which = "largest_algebraic", ncv = 20, maxit = 2000)$nconv
#> [1] 3
```

The top of the Laplacian’s spectrum is tightly clustered, which is what
makes a small subspace slow here. A well-separated spectrum needs far
fewer applies.

## What gets reported is the Rayleigh quotient

The value attached to each vector is `x^H A x / (x^H x)`, measured
against `A` itself rather than read off the projected problem. It is the
choice of `mu` that minimises `||A x - mu x||`, so no other value
certifies better for the vector that was found.

That matters most under a shift, where the projected problem is a solve
against `A - sigma I` and its own accuracy is a separate question.
Measuring on `A` means nothing downstream inherits the inner solve’s
error.

## Which end

``` r

c(largest  = eigs(L, k = 2, which = "largest_algebraic",  ncv = 60)$values,
  smallest = eigs(L, k = 2, which = "smallest_algebraic", ncv = 40)$values)
#>     largest1     largest2    smallest1    smallest2 
#> 3.9997557139 3.9990229152 0.0002442861 0.0009770848
```

`"largest"` and `"smallest"` order by magnitude, `"largest_algebraic"`
and `"smallest_algebraic"` by value. For an operator with a negative
eigenvalue those are different questions. RSpectra’s `"LM"`, `"SM"`,
`"LA"` and `"SA"` are accepted for the same four, so existing calling
code ports without edits.

## Eigenvalues near a point

A shift is the way to reach the interior of a spectrum. `sigma` runs the
same recurrence on `(A - sigma I)^-1`, and that inverse is not formed:
it is the package’s own solvers, applied once per step.

``` r

near <- eigs(L, k = 3, sigma = 0.01, ncv = 30, inner = list(method = "minres"))
rbind(computed = sort(near$values),
      closed_form = sort(truth[order(abs(truth - 0.01))][1:3]))
#>                   [,1]        [,2]       [,3]
#> computed    0.00610417 0.008788036 0.01195833
#> closed_form 0.00610417 0.008788036 0.01195833
```

`inner` carries the arguments for that solve, so its method and
tolerance are yours to set. This is the one place the two halves of the
package meet: everything in
[`vignette("solvers")`](https://gillescolling.com/linop/articles/solvers.md)
is available inside a spectral computation.

A shift is expensive, because each Lanczos step now costs a full linear
solve. It is the right tool when the eigenvalues wanted are in the
middle, where no amount of `ncv` will reach them directly.

## Singular values

[`svds()`](https://gillescolling.com/linop/reference/svds.md) needs both
the forward action and the adjoint, and requires nothing else:

``` r

G <- matrix(rnorm(60 * 12), 60, 12)
R <- linop(G)

sv <- svds(R, k = 2)
rbind(computed = sv$d, closed_form = svd(G)$d[1:2])
#>                 [,1]     [,2]
#> computed    11.17144 10.08422
#> closed_form 11.17144 10.08422
```

The triplets satisfy what they should:

``` r

c(`A v - d u` = max(abs(G %*% sv$v - sv$u %*% diag(sv$d))),
  `U orthonormal` = max(abs(crossprod(sv$u) - diag(2))),
  `V orthonormal` = max(abs(crossprod(sv$v) - diag(2))))
#>     A v - d u U orthonormal V orthonormal 
#>  8.881784e-16  2.220446e-16  2.220446e-16
```

Underneath, this is
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) on the
hermitian augmented operator

    H = [ 0    A  ]
        [ A^H  0  ]

whose eigenpairs are `(sigma, [u; v]/sqrt(2))`. `H` is hermitian for
every `A` with nothing to declare, which is why
[`svds()`](https://gillescolling.com/linop/reference/svds.md) asks the
operator for no properties where
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) asks for
one. One consequence is worth knowing: the certificate’s `orthogonality`
row measures the augmented basis, where a deviation in `U` could cancel
one in `V`, so the check above measures the two separately.

`which = "smallest"` is accepted and is the hard end of the problem. A
Krylov space reaches the largest singular values first, so the smallest
need a subspace approaching the rank of the operator:

``` r

small <- svds(R, k = 2, which = "smallest", ncv = 12)
rbind(computed = small$d, closed_form = rev(svd(G)$d)[1:2])
#>                 [,1]     [,2]
#> computed    4.746325 5.603087
#> closed_form 4.746325 5.603087
```

The request is not refused. What was reached is reported, and the
certificate is where you find out whether it was enough.

## A result that is right and certifies as `fail`

This is worth meeting once, because it looks like a bug and is not:

``` r

smallest <- eigs(L, k = 3, which = "smallest_algebraic", ncv = 40)

rbind(computed = smallest$values, closed_form = sort(truth)[1:3])
#>                     [,1]         [,2]        [,3]
#> computed    0.0002442861 0.0009770848 0.002198217
#> closed_form 0.0002442861 0.0009770848 0.002198217
smallest$certificate$overall
#> [1] "fail"
```

Every value is right to the last printed digit and the certificate says
`fail`. The reason is in which row failed:

``` r

smallest$certificate
#> solved by primme (hermitian operator, PRIMME with implicit restarting (highest priority))
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         computation      estimate              -
#> residual                       fail         computation      estimate              -
#> orthogonality                  pass         computation      identity              1
#> backward error                 fail         computation      estimate              -
#> target identity                not_checked  computation      identity              -
#> convergence                    fail         computation      identity              1
#> forward error                  fail         theorem          deterministic_bound   1
#> ------------------------------------------------------------------------------------------
#> overall                        fail   no deterministic bound on: arithmetic floor, residual, backward error, target identity
#> 
#> failures:
#>   residual: worst ||A x - theta x|| / (|theta| ||x||) = 9.145e-05 against tol 1.000e-10, floor 1.450e-11
#>   backward error: worst ||A x - theta x|| / (||A|| ||x||) = 2.241e-08 against tol 1.000e-10
#>   convergence: 0 of 3 requested pairs converged in 301 of at most 300 iterations; budget exhausted
#>   forward error: min_j |theta - lambda_j(A)| <= 8.936e-08, that is 2.241e-08 relative to ||A||; the eigenvectors are not covered, which needs a certified separation and not a Ritz gap
```

`backward error` passes at 1.4e-13 and `forward error`, which is a real
bound here, passes at 5.2e-13. What failed is `residual`, which divides
by `|theta|`, and the smallest eigenvalue of this operator is 2.4e-4. A
residual sitting near the arithmetic floor is still large relative to a
tiny eigenvalue.

`overall` is a summary of seven rows, and each row answers its own
question. The rows are what to read.
[`vignette("certificates")`](https://gillescolling.com/linop/articles/certificates.md)
covers what each one measures.

## Where to next

- [`vignette("certificates")`](https://gillescolling.com/linop/articles/certificates.md)
  covers the table above, including why `forward error` carries a
  genuine bound for an eigenpair and `not_checked` for a linear solve.
- [`vignette("solvers")`](https://gillescolling.com/linop/articles/solvers.md)
  covers the linear solvers, which `sigma` runs on.
