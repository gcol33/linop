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
few eigenpairs of a square operator.
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

fit <- eigs(L, k = 3, which = "largest_algebraic", method = "lanczos", ncv = 60)
fit
#> <lanczos> 3 of 3 requested eigenpairs converged, which = 'largest_algebraic'
#>   450 iterations of at most 4000, 12 restarts
#>   values: 3.999756  3.999023  3.997802
#>   certificate: qualified
```

``` r

rbind(computed = fit$values,
      closed_form = rev(sort(truth))[1:3])
#>                 [,1]     [,2]     [,3]
#> computed    3.999756 3.999023 3.997802
#> closed_form 3.999756 3.999023 3.997802
```

## Hermitian or not

[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) has two
methods. Lanczos is defined for a hermitian operator, and its projected
problem is a real symmetric tridiagonal matrix. Arnoldi with
Krylov-Schur restarting runs on any square operator. Symmetry is not
probed and not assumed, so a callback that declares nothing gets
Arnoldi:

``` r

Undeclared <- linop(lap, adjoint = lap, dim = c(n, n))
plain <- eigs(Undeclared, k = 3, which = "largest_algebraic", ncv = 60)
plain$method
#> [1] "arnoldi"
plain$values - rev(sort(truth))[1:3]
#> [1] 8.881784e-16 1.332268e-15 4.440892e-16
```

The values are right, and the dispatch says why Lanczos was passed over:

``` r

d <- plain$certificate$dispatch$candidates
d[d$backend %in% c("lanczos", "arnoldi"), c("backend", "status")]
#>   backend   status
#> 2 lanczos rejected
#> 3 arnoldi selected
d$reason[d$backend == "lanczos"]
#> [1] "requires the operator to be hermitian; it declares hermitian = NA, and unknown is not false"
```

Unknown is not false, and the dispatch says so rather than guessing
either way. A dense leaf establishes symmetry by checking it, so a
matrix needs no declaration; a callback has nothing to check and the
declaration is how you supply it.

A declaration is your word, and `method = "auto"` does not act on it
alone. It takes Lanczos when symmetry is computed, constructed or
proven, as [`solve()`](https://rdrr.io/r/base/solve.html) takes MINRES
only then, and otherwise runs Arnoldi and says why. Naming the method is
how you rely on your own declaration, which is why the calls on `L` in
this article say `method = "lanczos"`:

``` r

eigs(L, k = 3, which = "largest_algebraic", ncv = 60)$method
#> [1] "arnoldi"
d <- eigs(L, k = 3, which = "largest_algebraic", ncv = 60)$certificate$dispatch$candidates
d$reason[d$backend == "lanczos"]
#> [1] "auto acts on hermitian only at computed, constructed or proven evidence, and the operator's is user_declaration; name the method to rely on it"
```

Whichever method runs, the declaration buys a bound on the eigenvalues,
the `forward error` row of the certificate, and that row records the
declaration it rests on:

``` r

row <- function(fit, check) fit$certificate$checks$status[fit$certificate$checks$check == check]
c(declared = row(fit, "forward error"), undeclared = row(plain, "forward error"))
#>      declared    undeclared 
#>        "pass" "not_checked"
```

## A non-normal operator

Convection-diffusion on a grid,
`2 u_j - (1 + mu) u_{j-1} - (1 - mu) u_{j+1}`, is not symmetric for
`mu != 0`, and its eigenvalues are still known in closed form,
`2 - 2 sqrt(1 - mu^2) cos(j pi / (n + 1))`. Its eigenvectors are the
sine vectors scaled by `rho^j`, `rho = sqrt((1 + mu)/(1 - mu))`, so they
are far from orthogonal once `rho^n` is large.

``` r

convdiff <- function(n, mu) {
  force(n); force(mu)
  linop(function(X) 2 * X - (1 + mu) * rbind(0, X[-n, , drop = FALSE]) -
                    (1 - mu) * rbind(X[-1, , drop = FALSE], 0), dim = c(n, n))
}
cd_truth <- function(n, mu) 2 - 2 * sqrt(1 - mu^2) * cos(seq_len(n) * pi / (n + 1))

mild <- eigs(convdiff(80, 0.04), k = 3, which = "largest_real")
mild$values - sort(cd_truth(80, 0.04), decreasing = TRUE)[1:3]
#> [1] -7.993606e-15  3.552714e-14 -2.975398e-13
```

At `mu = 0.04` and `n = 80`, `rho^n` is 24 and the values agree with the
closed form to rounding. At `mu = 0.5` and `n = 100` it is `3^50`, and
the same computation goes differently:

``` r

hard <- eigs(convdiff(100, 0.5), k = 2, which = "largest_real")
hard$values - sort(cd_truth(100, 0.5), decreasing = TRUE)[1:2]
#> [1] 0.03611324+0.01361576i 0.03862589-0.01361577i
hard$certificate
#> solved by arnoldi (operator not established hermitian, Arnoldi with Krylov-Schur restarting (highest priority eligible; primme requires the operator to be hermitian; it declares hermitian = NA, and unknown is not false; lanczos requires the operator to be hermitian; it declares hermitian = NA, and unknown is not false))
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         computation      estimate              -
#> residual                       pass         computation      identity              1
#> orthogonality                  not_checked  computation      identity              -
#> backward error                 pass         computation      estimate              -
#> target identity                not_checked  computation      identity              -
#> convergence                    pass         computation      identity              1
#> forward error                  not_checked  computation      identity              -
#> ------------------------------------------------------------------------------------------
#> overall                        qualified   no deterministic bound on: arithmetic floor, orthogonality, backward error, target identity, forward error
```

The backward error passes: each value is an exact eigenvalue of an
operator within `tol ||A||` of this one. The eigenvalues of this
operator move by far more than that under such a perturbation, which is
what an eigenvector condition number of `3^50` means, so the computed
values are wrong in the first digit, and they come back as a complex
pair although every eigenvalue is real. The certificate does not claim
otherwise. Its `forward error` row places each value in the operator’s
pseudospectrum and stops there, because the bound from there to an
eigenvalue needs the eigenvector condition number, and that is not
computed. The `orthogonality` row is not checked either, since these
eigenvectors are not orthogonal.

A real operator’s non-real eigenvalues come in conjugate pairs, and
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) computes
the two members together, so asking for `k` values that end inside a
pair still returns both members converged. `which` takes
`"largest_imaginary"` and `"smallest_imaginary"` for the magnitude of
the imaginary part, and a complex `sigma` for the values nearest a point
off the real axis.

[`svds()`](https://gillescolling.com/linop/reference/svds.md) works on
any operator and is sometimes the thing that was actually wanted:
singular values are well conditioned whatever the eigenvectors are.

## `ncv` is the knob that matters

`ncv` is how many vectors a restart holds. Each restart compresses the
subspace to the vectors it keeps, so the smaller `ncv` is, the more
applies the method spends reaching the tolerance:

``` r

ladder <- t(vapply(c(20, 30, 40, 60), function(nc) {
  f <- eigs(L, k = 3, which = "largest_algebraic", method = "lanczos", ncv = nc)
  c(ncv = nc,
    error = max(abs(f$values - rev(sort(truth))[1:3])),
    converged = f$nconv,
    applies = f$iterations)
}, numeric(4)))
ladder
#>      ncv        error converged applies
#> [1,]  20 8.881784e-16         3     618
#> [2,]  30 4.440892e-16         3     493
#> [3,]  40 4.440892e-16         3     460
#> [4,]  60 8.881784e-16         3     450
```

The top of the Laplacian’s spectrum is tightly clustered, which is what
makes a small subspace slow here. A well-separated spectrum needs far
fewer applies.

The budget `maxit` counts applies, and its default is `20 n`: twenty
times what forming the operator column by column would cost. It grows
with `n` because a crowded end of the spectrum gets more crowded as `n`
grows, and at a fixed `ncv` the applies grow with it: the reference
Lanczos takes 1.6 `n` for the three largest here and 6 `n` at
`n = 2000`. At the defaults the run converges well inside it:

``` r

dflt <- eigs(L, k = 3, which = "largest_algebraic", method = "lanczos")
c(converged = dflt$nconv, applies = dflt$iterations, budget = dflt$maxit)
#> converged   applies    budget 
#>         3       549      4000
```

A run that stops improving ends before the budget does and says so, so a
larger default costs time only on a run that is still converging. The
far end of a spectrum spread over many decades is the case it does not
cover; a shift near that end, `sigma`, is the tool there.

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

c(largest  = eigs(L, k = 2, which = "largest_algebraic", method = "lanczos", ncv = 60)$values,
  smallest = eigs(L, k = 2, which = "smallest_algebraic", method = "lanczos", ncv = 40)$values)
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

near <- eigs(L, k = 3, sigma = 0.01, method = "lanczos", ncv = 30,
             inner = list(method = "minres"))
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

## A second operator on the right

A finite-element discretisation puts a mass matrix on the right-hand
side, `K u = lambda M u`. For linear elements on the unit interval the
stiffness is the Laplacian above and the mass is `tridiag(1, 4, 1) / 6`;
the two share the sine eigenvectors, so the pencil’s eigenvalues are
known in closed form, `6 (1 - cos t) / (2 + cos t)` with
`t = j pi / (n + 1)`:

``` r

mass <- function(X) (4 * X + rbind(X[-1, , drop = FALSE], 0) +
                     rbind(0, X[-n, , drop = FALSE])) / 6
M <- linop(mass, adjoint = mass, dim = c(n, n),
           properties = c(hermitian = TRUE, positive_definite = TRUE))
t_j <- seq_len(n) * pi / (n + 1)
pencil_truth <- sort(6 * (1 - cos(t_j)) / (2 + cos(t_j)), decreasing = TRUE)

gen <- eigs(L, k = 3, which = "largest_algebraic", method = "lanczos", B = M, ncv = 60)
gen$values - pencil_truth[1:3]
#> [1]  1.776357e-15 -3.552714e-15  0.000000e+00
```

`B` is the argument. With `A` hermitian and `B` hermitian positive
definite, Lanczos runs on `B^-1 A` in the inner product `x^H B x`, where
that operator is self-adjoint, so the projected problem stays real and
symmetric. Each step solves with `B`, by the package’s own solvers, and
`inner` carries that solve’s arguments as it does under a shift. The
vectors come back orthonormal in the same inner product:

``` r

max(abs(crossprod(gen$vectors, mass(gen$vectors)) - diag(3)))
#> [1] 1.110223e-15
```

Without the two declarations on `B` the same call runs Arnoldi on
`B^-1 A` and says why in the dispatch. The certificate measures each
pair as `x^H A x / x^H B x` and perturbs both operators in its backward
error, `||A x - theta B x|| / ((||A|| + |theta| ||B||) ||x||)`, since a
pencil’s data is the pair. `sigma` with `B` runs on
`(A - sigma B)^-1 B`.

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

## A result that is right and certifies as `qualified`

This is worth meeting once, because it looks like a hedge and is not:

``` r

smallest <- eigs(L, k = 3, which = "smallest_algebraic", method = "lanczos", ncv = 40)

rbind(computed = smallest$values, closed_form = sort(truth)[1:3])
#>                     [,1]         [,2]        [,3]
#> computed    0.0002442861 0.0009770848 0.002198217
#> closed_form 0.0002442861 0.0009770848 0.002198217
smallest$certificate$overall
#> [1] "qualified"
```

Every value is right to the last printed digit and the certificate says
`qualified`. The reason is in which row is qualified:

``` r

smallest$certificate
#> solved by lanczos (hermitian operator, Lanczos with full reorthogonalisation (requested by name))
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         computation      estimate              -
#> residual                       pass         computation      identity              1
#> orthogonality                  pass         computation      identity              1
#> backward error                 pass         computation      estimate              -
#> target identity                not_checked  computation      identity              -
#> convergence                    pass         computation      identity              1
#> forward error                  pass         theorem          deterministic_bound   1
#> ------------------------------------------------------------------------------------------
#> overall                        qualified   no deterministic bound on: arithmetic floor, backward error, target identity
```

`backward error` and `forward error`, which is a rigorous bound here,
pass with room to spare. The row that does not is `residual`, which
divides by `|theta|`, and the smallest eigenvalue of this operator is
2.4e-4. A residual near the arithmetic floor is still large relative to
a tiny eigenvalue, so the pair meets the tolerance only once the floor
is added to it, and `qualified` is the status that says so.

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
