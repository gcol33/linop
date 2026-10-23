# Solving a system

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

## One verb

[`solve()`](https://rdrr.io/r/base/solve.html) is R’s generic for
`a x = b`, and a `linop` has a method for it. There is no second verb to
learn and no separate object to construct.

``` r

n <- 200
L <- linop(
  function(X) 2 * X - rbind(X[-1, , drop = FALSE], 0) - rbind(0, X[-n, , drop = FALSE]),
  dim = c(n, n)
)
xstar <- sin(seq_len(n) * pi / (n + 1))
b <- as.numeric(L %*% xstar)

x <- solve(L, b)
max(abs(x - xstar))
#> [1] 1.061568e-12
```

`L` is the 1D Laplacian written as its action. Nothing was factorized,
because for an operator given as a closure there is nothing to
factorize. What came back is a plain numeric vector, so it indexes and
multiplies like any other, and it carries a record of how it was
produced:

``` r

attr(x, "certificate")
#> solved by gmres (no capability is established, and gmres requires none)
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         computation      estimate              -
#> residual                       pass         computation      identity              1
#> backward error                 pass         computation      estimate              -
#> convergence                    pass         computation      identity              1
#> forward error                  not_checked  computation      identity              -
#> ------------------------------------------------------------------------------------------
#> overall                        qualified   no deterministic bound on: arithmetic floor, backward error, forward error
```

The first line is the one to read here. The rest of that table is the
subject of
[`vignette("certificates")`](https://gillescolling.com/linop/articles/certificates.md).

## What `auto` chose, and why

`method = "auto"` picks from what the operator has established about
itself and from its shape, and records both the choice and the reason.
Three routes are reachable.

A diagonal operator of full rank is solved exactly, with no iteration at
all:

``` r

D <- linop_scaling(seq(1, 4, length.out = 50))
attr(solve(D, rep(1, 50)), "certificate")$dispatch$reason
#> [1] "the operator is diagonal and of full rank, so the solve is exact"
```

A dense leaf checks its own symmetry when it is built, on data it
already holds, so a symmetric matrix reaches MINRES:

``` r

M <- crossprod(matrix(rnorm(60 * 60), 60, 60)) + diag(60)
Ad <- linop(M)
bd <- as.numeric(Ad %*% rep(1, 60))
attr(solve(Ad, bd), "certificate")$dispatch$reason
#> [1] "the operator is hermitian, and definiteness is not established"
```

The operator above is positive definite, and `auto` still chose MINRES
rather than CG. Symmetry was computed; definiteness was not, and `auto`
will not assume it. There is a CG branch above the MINRES one, and no
constructor in the package currently reaches it, because the routes that
establish definiteness at a source `auto` accepts are diagonal and are
taken by the direct route first. CG is reached by naming it.

That is the general rule, and it applies to declarations too. `L` at the
top of this article was built with no properties at all. Declaring them
changes nothing about which method `auto` picks:

``` r

Ldecl <- linop(
  function(X) 2 * X - rbind(X[-1, , drop = FALSE], 0) - rbind(0, X[-n, , drop = FALSE]),
  dim = c(n, n),
  properties = c(hermitian = TRUE, positive_definite = TRUE)
)
attr(solve(Ldecl, b), "certificate")$dispatch$reason
#> [1] "no capability is established, and gmres requires none"
```

Both capabilities are `TRUE` and `auto` acted on neither, because a bare
declaration is not the kind of argument `auto` will choose a method on.
The capability and the evidence behind it are two separate questions,
which is what
[`vignette("capabilities")`](https://gillescolling.com/linop/articles/capabilities.md)
is about. Here it decides a dispatch.

GMRES is the fallback, and it is what makes `auto` total rather than
partial: it requires nothing of the operator, and it runs on the forward
action alone. `L` has no adjoint, which is why every route above it was
unavailable in the first place.

## Naming a method

Naming a method is the caller making the assertion, so no evidence
minimum applies:

``` r

attr(solve(Ldecl, b, method = "cg"), "certificate")$dispatch
#> $requested
#> [1] "cg"
#> 
#> $chosen
#> [1] "cg"
#> 
#> $reason
#> [1] "named by the caller, so no evidence minimum was applied"
```

The declaration on `Ldecl` is true, CG is the right method for it, and
saying so is how you get it. The two questions stay separate: `auto` is
the package choosing on your behalf, and a named method is your own
claim about your own operator.

A declaration that is false is not silently tolerated. CG forms
`p^H A p` every step, and a value that is not positive far from
convergence means the operator is not what it said it was, which is an
error naming the capability rather than a wrong answer.

## The seven methods

| method   | requires                    | adjoint | shape  |
|:---------|:----------------------------|:--------|:-------|
| cg       | hermitian positive definite | no      | square |
| minres   | hermitian                   | no      | square |
| gmres    | nothing                     | no      | square |
| fgmres   | nothing                     | no      | square |
| bicgstab | nothing                     | no      | square |
| lsqr     | nothing                     | yes     | any    |
| lsmr     | nothing                     | yes     | any    |

GMRES and BiCGSTAB both require nothing and both run without an adjoint.
They are not interchangeable. GMRES minimises the residual over the
whole Krylov space, so its residual cannot rise and it has no breakdown;
BiCGSTAB is cheaper per step and can break down, in which case it
returns a `fail` certificate naming the breakdown rather than throwing.
`auto` prefers GMRES for the first reason.

FGMRES is GMRES with the one difference that it accepts a preconditioner
that changes between applications, such as an inner iterative solve.

## Least squares is a different request

[`solve()`](https://rdrr.io/r/base/solve.html) with `method = "auto"`
refuses a rectangular operator, and says what to name instead:

``` r

G <- matrix(rnorm(120 * 40), 120, 40)
R <- linop(G)
bl <- rnorm(120)
solve(R, bl)
#> Error:
#> ! solve() needs a square operator; this one is 120 x 40.
#>   A rectangular system is a least-squares problem, min ||b - a x||, which is a
#>   different request rather than the same one on a different shape. Name the method
#>   to make it: solve(a, b, method = "lsqr") or method = "lsmr".
```

This is deliberate. `min ||b - A x||` is a different mathematical
question from `A x = b`, rather than the same one on a different shape,
so the caller states which they meant. Both methods are reachable
through the same verb:

``` r

xl <- solve(R, bl, method = "lsqr")
max(abs(xl - qr.solve(G, bl)))
#> [1] 3.203147e-09
```

LSQR and LSMR solve the same problem and minimise different things over
the same bidiagonalisation. LSQR minimises `||r||`. LSMR minimises
`||A^H r||`, which is the quantity that goes to zero at a least-squares
solution and the numerator of the backward error the certificate
reports:

``` r

xm <- solve(R, bl, method = "lsmr")
c(lsqr = attr(xl, "certificate")$values$backward_error,
  lsmr = attr(xm, "certificate")$values$backward_error)
#>         lsqr         lsmr 
#> 4.008744e-09 3.731806e-09
```

## Preconditioners

A preconditioner is built from one direction, `M^-1`, and carries what
it declares about itself:

``` r

set.seed(2)
d <- 10^runif(200, 0, 5)
Adense <- diag(d)
Adense[cbind(1:199, 2:200)] <- -1
Adense[cbind(2:200, 1:199)] <- -1
A <- linop(Adense)
xs <- rep(1, 200)
bb <- as.numeric(A %*% xs)

P <- preconditioner(function(R) R / d, hermitian = TRUE,
                    positive_definite = TRUE, side = "left", dim = c(200, 200))
P
#> <preconditioner>
#>   side:              left
#>   fixed:             TRUE
#>   hermitian:         TRUE
#>   positive_definite: TRUE
#>   M^-H:              equal to M^-1, by the hermitian declaration
```

The diagonal of that operator spans five orders of magnitude. At a fixed
budget of 40 iterations, with everything else equal:

``` r

plain  <- solve(A, bb, method = "minres", maxit = 40, tol = 1e-10)
scaled <- solve(A, bb, method = "minres", maxit = 40, tol = 1e-10,
                preconditioner = P)

c(plain  = max(abs(plain  - xs)),
  scaled = max(abs(scaled - xs)))
#>        plain       scaled 
#> 1.007307e+00 1.294592e-07
c(plain  = attr(plain,  "certificate")$overall,
  scaled = attr(scaled, "certificate")$overall)
#>       plain      scaled 
#>      "fail" "qualified"
```

Each method states which properties it needs of a preconditioner and
which sides it is defined for, and enforces both at the point the
preconditioner is accepted:

| method   | requires                            | sides              |
|:---------|:------------------------------------|:-------------------|
| cg       | fixed, hermitian, positive definite | left, right, split |
| minres   | fixed, hermitian, positive definite | left, right, split |
| gmres    | fixed                               | left, right, split |
| fgmres   | nothing                             | right              |
| bicgstab | fixed                               | left, right, split |
| lsqr     | fixed                               | right              |
| lsmr     | fixed                               | right              |

A refusal names what was missing:

``` r

flexible <- preconditioner(function(R) R / d, fixed = FALSE,
                           hermitian = TRUE, positive_definite = TRUE,
                           side = "left", dim = c(200, 200))
solve(A, bb, method = "minres", preconditioner = flexible)
#> Error:
#> ! method 'minres' requires a preconditioner with fixed = TRUE; this one has fixed = FALSE.
#>   A flexible preconditioner changes between applications; use fgmres, which takes a right-side preconditioner.
```

For CG the three sides coincide, so its row is unrestricted because the
distinction does not arise. For GMRES all three are implemented and they
produce different iterates, so its row is unrestricted for the opposite
reason:

``` r

set.seed(1)
Gm <- matrix(rnorm(80 * 80), 80, 80)
diag(Gm) <- diag(Gm) + 12
Ag <- linop(Gm)
bg <- as.numeric(Ag %*% rep(1, 80))
dg <- diag(Gm)

vapply(c("left", "right", "split"), function(sd) {
  Ps <- preconditioner(function(R) R / dg, hermitian = TRUE,
                       positive_definite = TRUE, side = sd, dim = c(80, 80))
  attr(solve(Ag, bg, method = "gmres", maxit = 6, tol = 1e-12,
             preconditioner = Ps), "certificate")$values$residual
}, numeric(1))
#>       left      right      split 
#> 0.04805537 0.04799045 0.04800726
```

Right preconditioning minimises the residual that gets reported. Left
and split minimise a preconditioned one, so the outer loop converts the
target going in and measures in your units coming out.

The narrow rows are the two least-squares methods, and the reason is the
shape. `M^-1` acts on the domain of the operator; left and split act on
the codomain, which for a rectangular operator is a different space.

## A solve that does not converge

An exhausted budget is a result, not an exception. It comes back as a
`fail` certificate that says what ran out:

``` r

attr(plain, "certificate")
#> solved by minres (named by the caller, so no evidence minimum was applied)
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         computation      identity              1
#> residual                       fail         computation      estimate              -
#> backward error                 fail         computation      estimate              -
#> convergence                    fail         computation      identity              1
#> forward error                  not_checked  computation      identity              -
#> ------------------------------------------------------------------------------------------
#> overall                        fail   no deterministic bound on: residual, backward error, forward error
#> 
#> failures:
#>   residual: worst ||b - A x|| / ||b|| = 7.135e-04 against tol 1.000e-10, floor 3.360e-15
#>   backward error: worst ||r|| / (||A|| ||x|| + ||b||) = 1.886e-04 against tol 1.000e-10
#>   convergence: 40 of at most 40 iterations; budget exhausted
```

The iterate is finite and usable, the number reached is reported, and
nothing was thrown. A contradicted declaration is the case that does
throw, because there the operator is not what it claimed and no iterate
is worth returning.

`details = TRUE` returns the whole object rather than the solution,
including the convergence history:

``` r

fit <- solve(A, bb, method = "minres", maxit = 40, tol = 1e-10, details = TRUE)
names(fit)
#>  [1] "x"           "certificate" "method"      "iterations"  "restarts"   
#>  [6] "converged"   "residual"    "history"     "requested"   "reason"
fit$iterations
#> [1] 40
fit$converged
#> [1] FALSE
```

## Where to next

- [`vignette("certificates")`](https://gillescolling.com/linop/articles/certificates.md)
  covers the table printed above: what each row measures, the arithmetic
  floor underneath it, and which rows carry a bound rather than a
  measurement.
- [`vignette("spectral")`](https://gillescolling.com/linop/articles/spectral.md)
  covers [`eigs()`](https://gillescolling.com/linop/reference/eigs.md)
  and [`svds()`](https://gillescolling.com/linop/reference/svds.md),
  which are built on the solvers here when a shift is given.
- [`vignette("capabilities")`](https://gillescolling.com/linop/articles/capabilities.md)
  covers what an operator establishes about itself, which is what `auto`
  reads.
