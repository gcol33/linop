# Reading a certificate

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

## One object, five shapes

Every computation in this package that could be wrong comes back with a
certificate. It is the same class each time, and it takes five shapes,
because the questions differ:

| produced_by | subject |
|:---|:---|
| verify() | does the operator behave like a linear map |
| solve(), square | how close is x to solving A x = b |
| solve(), least squares | how close is x to minimising \|\|b - A x\|\| |
| eigs(), svds() | how close is (theta, x) to an eigenpair |
| eigs() on a sequence space | how close is theta to the spectrum of the operator that was truncated |

This article is about the last four.
[`verify()`](https://gillescolling.com/linop/reference/verify.md) is
covered in
[`vignette("adapters")`](https://gillescolling.com/linop/articles/adapters.md),
where it belongs to the job of writing an operator.

## Getting hold of one

Each verb puts its certificate where its own return value has room. A
solve returns the solution a matrix solve would return and carries the
certificate in an attribute; the spectral verbs return a result object
with a field;
[`verify()`](https://gillescolling.com/linop/reference/verify.md)
returns the certificate itself.
[`certificate()`](https://gillescolling.com/linop/reference/certificate.md)
reads it off any of them.

``` r

n <- 200
lap <- function(X) 2 * X - rbind(X[-1, , drop = FALSE], 0) - rbind(0, X[-n, , drop = FALSE])
L <- linop(lap, adjoint = lap, dim = c(n, n), properties = c(hermitian = TRUE))
b <- as.numeric(L %*% sin(seq_len(n) * pi / (n + 1)))

cert <- certificate(solve(L, b, method = "minres", tol = 1e-10))
cert
#> solved by minres (named by the caller, so no evidence minimum was applied)
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

A certificate belongs to a result rather than to an operator. An
operator on its own has nothing to certify: the one claim it can make
alone is that it satisfies the contract, and that is
[`verify()`](https://gillescolling.com/linop/reference/verify.md), which
returns a certificate of its own.

The attribute on a solve survives arithmetic and is dropped by indexing,
[`as.numeric()`](https://rdrr.io/r/base/numeric.html),
[`c()`](https://rdrr.io/r/base/c.html) and the reductions. So the cost
of [`solve()`](https://rdrr.io/r/base/solve.html) returning a plain
matrix is not that the certificate is easily lost; it is that `2 * x`
still carries one describing `x`.

## Status and evidence are separate questions

Every row of that table carries a status and, separately, the kind of
argument behind it.

The status column answers whether the check passed. `pass`, `fail`,
`not_checked`, and `qualified`, which means it passed for a reason
weaker than the clean one.

The three columns after it are the evidence, and they are three
independent fields rather than a ranking. `source` is what kind of
argument produced the number: a `computation` on data at hand, a
`theorem`, a `construction`, a `user_declaration`. `guarantee` is what
the argument gives you: an `identity` holds exactly, a
`deterministic_bound` holds always, an `estimate` holds usually.
`confidence` is a probability where one applies and `-` where none does.

A `pass` behind an `estimate` and a `pass` behind an `identity` are both
passes and they are not the same claim. The summary line names the rows
with neither kind of deterministic guarantee, so that distinction
survives being skimmed:

``` r

cert$without_deterministic_bound
#> [1] "arithmetic floor" "backward error"   "forward error"
```

## The rows of a linear solve

| row | measures |
|:---|:---|
| arithmetic floor | the level below which no residual is meaningful |
| residual | \|\|b - A x\|\| / \|\|b\|\| |
| backward error | \|\|r\|\| / (\|\|A\|\| \|\|x\|\| + \|\|b\|\|), the Rigal-Gaches quantity |
| convergence | whether the iteration stopped because it was done |
| forward error | how far x is from the true solution |

The backward error is the row to read. It answers a question that has an
exact answer: it is the size of the smallest perturbation to `A` and `b`
for which your computed `x` is the exact solution. A backward error of
1e-14 means you solved a problem within 1e-14 of the one you asked
about.

`||A||` in that denominator is a lower bound rather than an estimate in
either direction, deliberately. Shrinking the denominator can only make
the reported backward error larger, so the number can overstate the
error and never understate it.

## The arithmetic floor

Ask for an exact solve and no floating-point iteration can deliver it:

``` r

tight <- solve(L, b, method = "minres", tol = 0, maxit = 400)
bare  <- solve(L, b, method = "minres", tol = 0, maxit = 400,
               control = list(floor_const = 0))

identical(as.numeric(tight), as.numeric(bare))
#> [1] TRUE
```

The same iteration ran twice and produced the same vector, bit for bit.
Only the certificates differ:

``` r

c(with_floor = attr(tight, "certificate")$overall,
  without    = attr(bare,  "certificate")$overall)
#>  with_floor     without 
#> "qualified"      "fail"
```

``` r

attr(bare, "certificate")
#> solved by minres (named by the caller, so no evidence minimum was applied)
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         computation      estimate              -
#> residual                       fail         computation      estimate              -
#> backward error                 fail         computation      estimate              -
#> convergence                    fail         computation      identity              1
#> forward error                  not_checked  computation      identity              -
#> ------------------------------------------------------------------------------------------
#> overall                        fail   no deterministic bound on: arithmetic floor, residual, backward error, forward error
#> 
#> failures:
#>   residual: worst ||b - A x|| / ||b|| = 9.961e-13 against tol 0.000e+00, floor 0.000e+00
#>   backward error: worst ||r|| / (||A|| ||x|| + ||b||) = 6.084e-17 against tol 0.000e+00
#>   convergence: 400 of at most 400 iterations; budget exhausted
```

A backward error of 1.8e-16 is reported as a failure, because it was
tested against a tolerance of zero with nothing underneath it. The floor
is what stops a fully converged solve from certifying as `fail`, and it
is `c * eps * (||A|| ||x|| + ||b||)`, which is Higham’s bound on the
computed residual. That is the same quantity as the Rigal-Gaches
denominator, so once the backward error is written relatively the floor
is a plain `c * eps`.

A row that meets its tolerance only because of the floor says
`qualified` and carries an `estimate`, never a clean `identity`. `c`
counts terms in an inner product a matrix-free operator never exposes,
so it is an argument with a default rather than a constant.

## The same table, measured differently

A least-squares solve has no exact answer. `b - A x` does not go to
zero, so testing it against a tolerance would report a converged
solution as a failure. What goes to zero is `A^H r`.

Both are backward errors, so one row carries both readings, and which
one applies is decided by measurement rather than by the method that was
called:

``` r

G <- matrix(rnorm(120 * 40), 120, 40)
R <- linop(G)

incompatible <- solve(R, rnorm(120), method = "lsmr", tol = 1e-10)
compatible   <- solve(R, as.numeric(G %*% rnorm(40)), method = "lsmr", tol = 1e-10)

status_of <- function(x, row) {
  ch <- attr(x, "certificate")$checks
  ch$status[ch$check == row]
}
c(incompatible = status_of(incompatible, "residual"),
  compatible   = status_of(compatible,   "residual"))
#>  incompatible    compatible 
#> "not_checked"        "pass"
```

For the incompatible right-hand side the residual reads `not_checked`,
because it has a nonzero minimum and comparing that to `tol` would mean
nothing. For the compatible one the same code found a residual that did
go to zero and checked it.

``` r

attr(incompatible, "certificate")
#> solved by lsmr (named by the caller, so no evidence minimum was applied)
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         computation      identity              1
#> residual                       not_checked  computation      estimate              -
#> backward error                 pass         computation      identity              1
#> convergence                    pass         computation      identity              1
#> forward error                  not_checked  computation      identity              -
#> ------------------------------------------------------------------------------------------
#> overall                        qualified   no deterministic bound on: residual, forward error
```

The perturbation behind the least-squares reading is Stewart’s,
`dA = -(r r^H A) / ||r||^2`, which is exhibited rather than argued: it
is an actual matrix of that norm making `x` an exact minimiser. That is
what makes `||A^H r|| / (||A|| ||r||)` a backward error rather than a
heuristic, and it is also Paige and Saunders’ stopping rule.

## A third shape, for an eigenpair

For a solve, `b` is data and `A` is data, and the Rigal-Gaches
denominator has a term for each. For an eigenpair, `theta` and `x` are
both outputs and `A` is the only datum, so the denominator loses a term
and the table gains rows:

``` r

fit <- eigs(L, k = 3, which = "largest_algebraic", ncv = 60)
fit$certificate
#> solved by primme (hermitian operator, PRIMME with implicit restarting (highest priority))
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

`orthogonality` is new: the computed vectors should be orthonormal, and
this measures whether they are. `target identity` is new and is always
`not_checked`, which is the honest entry. A small residual proves the
pair is an approximate eigenpair. That it is the *largest* one needs an
inertia count, an enclosure or a separation bound, and matrix-free there
is generally none available.

## The one row that is sometimes a bound

`forward error` is `not_checked` on every linear solve in this article,
and it is not laziness. Nothing constrains it. At a condition number of
1e10:

``` r

m <- 40
Q <- qr.Q(qr(matrix(rnorm(m * m), m, m)))
M <- Q %*% diag(10^seq(0, -10, length.out = m)) %*% t(Q)
M <- (M + t(M)) / 2
Aill <- linop(M)
xtrue <- rnorm(m)
bill <- as.numeric(Aill %*% xtrue)

t(vapply(c("minres", "gmres", "bicgstab"), function(mth) {
  z <- solve(Aill, bill, method = mth, tol = 1e-10, maxit = 400)
  c(backward = attr(z, "certificate")$values$backward_error,
    forward  = max(abs(z - xtrue)) / max(abs(xtrue)))
}, numeric(2)))
#>              backward   forward
#> minres   6.433681e-11 0.3207462
#> gmres    4.258965e-11 0.2493074
#> bicgstab 2.537527e-09 0.6736162
```

Backward errors around 1e-11 sit beside true forward errors around 1e-1.
Both columns are correct. `kappa^2 * eps` is 2.2e4 here, and no backward
error constrains the forward one at that conditioning, which is exactly
why the row refuses to report a number it cannot stand behind. Do not
read that spread as one method being more accurate than another.

For an eigenpair the same row does carry a bound, and it is the first in
the package. `A + E` with `E = -r x^H - x r^H` is hermitian and has
`theta` as an exact eigenvalue, so Weyl’s inequality gives
`min_j |theta - lambda_j(A)| <= ||r|| / ||x||`:

``` r

fwd <- fit$certificate$checks
fwd <- fwd[fwd$check == "forward error", ]
fwd[, c("status", "source", "guarantee")]
#>   status  source           guarantee
#> 7   pass theorem deterministic_bound
cat(fwd$detail)
#> min_j |theta - lambda_j(A)| <= 3.969e-10, that is 9.924e-11 relative to ||A||; the eigenvectors are not covered, which needs a certified separation and not a Ritz gap
```

`source` is `theorem` and `guarantee` is `deterministic_bound`, the only
row of this shape where either appears. The bound is about *some*
eigenvalue of `A`; which one is the `target identity` question and stays
`not_checked`.

## Evidence can carry a dependency

That bound holds because `A` is hermitian. If the only argument for that
is your say-so, the bound inherits it, and the certificate records the
inheritance rather than hiding it:

``` r

Md <- crossprod(matrix(rnorm(20 * 20), 20, 20)) + diag(20)

checked  <- eigs(linop(Md), k = 2, which = "largest_algebraic")
declared <- eigs(linop(function(X) Md %*% X, adjoint = function(X) Md %*% X,
                       dim = c(20, 20), properties = c(hermitian = TRUE)),
                 k = 2, which = "largest_algebraic")

checked$certificate$evidence[["forward error"]]
#> theorem / deterministic_bound / conf 1 <- [computation / identity / conf 1]
declared$certificate$evidence[["forward error"]]
#> theorem / deterministic_bound / conf 1 <- [user_declaration / identity / conf 1]
```

Both say `theorem / deterministic_bound`. What differs is in the
brackets, and it is what a requirement tests:

``` r

strict <- requirement(sources = c("construction", "computation", "theorem"),
                      guarantees = c("identity", "deterministic_bound"))

c(checked  = evidence_satisfies(checked$certificate$evidence[["forward error"]],  strict),
  declared = evidence_satisfies(declared$certificate$evidence[["forward error"]], strict))
#>  checked declared 
#>     TRUE    FALSE
```

A theorem applied to a declaration is worth what the declaration is
worth. This is the same rule that stops a sum of two declared
positive-definite operators from satisfying CG’s requirement, described
in
[`vignette("capabilities")`](https://gillescolling.com/linop/articles/capabilities.md),
reaching the certificate.

## A fifth shape, about an operator that is not there

Every shape so far certifies a computation against the operator the
computation ran on. A finite section certifies a statement about the
operator it truncated, which no computation touched:

``` r

H <- linop_jacobi(diagonal = 1)
dim(H)
#> [1] Inf Inf

fit <- eigs(H, k = 1, which = "largest_algebraic")
certificate(fit)
#> solved by primme on a finite section (operator on l^2(Z); width 41 -> 49, stopped because the tail is below half the tolerance)
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

The width was not given either.
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) on an
operator with no matrix truncates it, and the header records which
widths it went through and why it stopped;
[`vignette("sequence-spaces")`](https://gillescolling.com/linop/articles/sequence-spaces.md)
covers how that width is chosen.

`residual` and `backward error` are gone, because a truncated eigenvalue
answers no equation with data in it. What replaces them is a split, and
the split is the point:

``` r

v <- fit$certificate$values
c(finite_residual = v$residual, truncation = v$truncation, bound = v$bound)
#> finite_residual      truncation           bound 
#>    1.653582e-10    5.313083e-11    1.736869e-10
```

The residual of the extended eigenvector against the infinite operator
is those two terms in quadrature, exactly. The first is what the
eigensolve left behind and answers to `tol` and `ncv`; the second is
what the truncation costs and answers to `n`. Reported as one number you
would know the size of the error without knowing which knob moves it:

``` r

narrow <- eigs(finite_section(H, n = 12), k = 1, which = "largest_algebraic")
nv <- narrow$certificate$values
c(finite_residual = nv$residual, truncation = nv$truncation)
#> finite_residual      truncation 
#>    7.913351e-11    1.815373e-03
```

At `n = 12` the eigensolve is four orders better than it needs to be and
the section is six orders short, and the convergence row says which:

``` r

ch <- narrow$certificate$checks
cat(ch$detail[ch$check == "convergence"])
#> 0 of 1 requested pairs converged in 33 of at most 500 iterations; the eigensolve met its target and the section did not; the truncation term is 1.815e-03 at n = 12
```

`isolation` is the row with no counterpart anywhere else. The essential
spectrum of this class is known exactly rather than estimated, so a
value can be tested against it, and a value that has not cleared it is a
discretisation of continuous spectrum rather than an eigenvalue:

``` r

free <- eigs(finite_section(linop_jacobi(diagonal = 0), n = 100), k = 1,
             which = "largest_algebraic")
free$values
#> [1] 1.999758
free$certificate$checks[, c("check", "status")]
#>              check      status
#> 1 arithmetic floor        pass
#> 2  finite residual        pass
#> 3    orthogonality        pass
#> 4 truncation bound        fail
#> 5        isolation        fail
#> 6        enclosure not_checked
#> 7      convergence        fail
#> 8    forward error        fail
```

The operator with no potential has no eigenvalues at all. The
eigensolver does not notice: it converges cleanly and hands back a
number just inside the band. Only the isolation row separates *no
eigenvalue here* from *a hard eigenvalue here*, and it is what makes the
whole certificate refuse.

The `forward error` row here rests on nothing declared. The class is
self-adjoint by construction and the residual identity is arithmetic, so
unlike the eigenpair version it carries no dependency at all:

``` r

fit$certificate$evidence[["forward error"]]
#> theorem / deterministic_bound / conf 1
```

## Reading one in code

The printed table is for people. The object is a list, so a package can
test against it:

``` r

names(cert)
#> [1] "subject"                     "checks"                     
#> [3] "overall"                     "without_deterministic_bound"
#> [5] "probes"                      "values"                     
#> [7] "evidence"                    "dispatch"
cert$checks[, c("check", "status")]
#>              check      status
#> 1 arithmetic floor        pass
#> 2         residual        pass
#> 3   backward error        pass
#> 4      convergence        pass
#> 5    forward error not_checked
cert$values$backward_error
#> [1] 4.895612e-15
```

``` r

stopifnot(cert$overall %in% c("pass", "qualified"))
```

`cert$dispatch` records which method ran and why, `cert$values` holds
the measured quantities, and `cert$evidence` holds the evidence objects
for rows that carry a dependency.

## Where to next

- [`vignette("solvers")`](https://gillescolling.com/linop/articles/solvers.md)
  covers the methods that produce the first two shapes.

- [`vignette("spectral")`](https://gillescolling.com/linop/articles/spectral.md)
  covers [`eigs()`](https://gillescolling.com/linop/reference/eigs.md)
  and [`svds()`](https://gillescolling.com/linop/reference/svds.md),
  including a result whose values are exact and whose certificate reads
  `qualified`.

- [`vignette("capabilities")`](https://gillescolling.com/linop/articles/capabilities.md)
  covers evidence away from certificates, where the same three fields
  decide what a solver will accept.

- [`vignette("adapters")`](https://gillescolling.com/linop/articles/adapters.md)
  covers
  [`verify()`](https://gillescolling.com/linop/reference/verify.md), the
  fourth shape.

- [`vignette("sequence-spaces")`](https://gillescolling.com/linop/articles/sequence-spaces.md)
  covers the operators the fifth shape is about.
