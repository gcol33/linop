# Operators on a sequence space

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

## An operator with no matrix and no size

A `linop` is defined by its action rather than its entries, and a
dimension of `Inf` is what that promise looks like taken all the way.
[`linop_jacobi()`](https://gillescolling.com/linop/reference/linop_jacobi.md)
builds a self-adjoint operator on `l^2(Z)`, the square-summable
sequences indexed by the integers:

``` r

H <- linop_jacobi(diagonal = 1)
H
#> <linop> Inf x Inf, double
#>   hermitian, symmetric, real
#> 
#> jacobi on l^2(Z), sigma_ess [-2, 2]
dim(H)
#> [1] Inf Inf
```

It is tridiagonal in the standard basis:

`(H u)[j] = u[j-1] + u[j+1] + V[j] u[j]`

with `V` here the potential `1` at the origin and zero elsewhere. Both
coefficient sequences are real and *eventually constant*: they take the
values you supply on a finite window and their limiting values
everywhere else.

Nothing numeric runs on it, and the refusals say so by name:

``` r

as.matrix(H)
#> Error:
#> ! as.matrix() needs an operator with finite dimensions; this one is Inf x Inf.
#>   An operator on a sequence space has no matrix and no block to act on.
#>   Truncate it with finite_section() and call as.matrix() on the result.
```

`Inf` is the extent of the index set here, not a flag. The ordinary
shape arithmetic keeps working on it with no special case, so the
algebra composes these operators as it composes any other:

``` r

dim(t(H))
#> [1] Inf Inf
dim(2 * H + H)
#> [1] Inf Inf
```

## Truncating it

`finite_section(H, n)` is `P_n H P_n` on the indices `-n` to `n`. What
comes back is an ordinary operator of dimension `2n + 1`:

``` r

F <- finite_section(H, n = 20)
F
#> <linop> 41 x 41, double
#>   hermitian, symmetric, real
#> 
#> section |j| <= 20
#> `- jacobi on l^2(Z), sigma_ess [-2, 2]
dim(F)
#> [1] 41 41
```

The printed tree shows the infinite operator underneath the finite one,
because that is what the section holds. It is not an annotation: the
child is the operator itself, and everything the certificate says about
it is read off that child.

``` r

explain(F)
#> operator 41 x 41, dtype double, estimated 205 flops per column
#> 2 nodes, depth 1
#> 
#>  depth    node nrow ncol  dtype cost
#>      0 section   41   41 double  205
#>      1  jacobi  Inf  Inf double  Inf
```

A section is an ordinary `linop` in every other respect. It has a
matrix, the solvers run on it, and it composes:

``` r

round(as.matrix(finite_section(H, n = 2)), 3)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    1    0    0    0
#> [2,]    1    0    1    0    0
#> [3,]    0    1    1    1    0
#> [4,]    0    0    1    0    1
#> [5,]    0    0    0    1    0
x <- solve(F, rep(1, 41), method = "minres")
max(abs(as.numeric(F %*% x) - 1))
#> [1] 4.318546e-12
```

`n` has to be large enough that the window sits strictly inside it.
Outside the window the coefficients are at their limits and the
eigenvalue equation is the free recurrence, which is what makes
everything below a closed form:

``` r

finite_section(linop_jacobi(diagonal = c(1, 2, 3)), n = 2)
#> Error:
#> ! n = 2 cuts through the window; it has to exceed 2.
#>   Outside the window the coefficients are at their limits and the eigenvalue
#>   equation is the free recurrence, which is what makes the truncation error a
#>   closed form. A section that keeps part of the perturbation outside itself
#>   has no such statement behind it.
```

## Two things that are known exactly

The reason this class comes first is that two quantities are exact
rather than estimated.

The perturbation is finitely supported, hence finite rank, hence
compact, so by Weyl’s theorem the essential spectrum is that of the free
operator. The free operator has symbol `a + 2 b cos(theta)`, so:

`sigma_ess(H) = [a - 2|b|, a + 2|b|]`

which for this operator is `[-2, 2]`, from the coefficients and with
nothing computed.

The second is the decay. Outside the window an eigenfunction satisfies
the free recurrence exactly, whose solutions are `z^j` for the two roots
of `z^2 - mu z + 1` with `mu = (lambda - a)/b`. The roots are
reciprocal, so for a value outside the band exactly one is inside the
unit disc, and square summability forces it on each side:

``` r

decay_rate(H, sqrt(5))
#> [1] 0.618034
(sqrt(5) - 1) / 2
#> [1] 0.618034
```

The eigenfunction is *exactly* geometric outside the window at that
rate, not asymptotically so. Inside the band both roots sit on the unit
circle, nothing decays, and there is nothing to report:

``` r

decay_rate(H, c(0, 1.9, 2, 3))
#> [1]       NA       NA       NA 0.381966
```

## Solving it

[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) on the
section certifies a statement about `H`:

``` r

fit <- eigs(F, k = 1, which = "largest_algebraic")
fit$values
#> [1] 2.236068
```

For a single site at the origin the answer is closed form,
`lambda = sign(v) sqrt(v^2 + 4)`, so the error is available:

``` r

c(computed = fit$values, exact = sqrt(5), error = abs(fit$values - sqrt(5)))
#>     computed        exact        error 
#> 2.236068e+00 2.236068e+00 1.493015e-09
```

and the certificate bounded it without knowing it:

``` r

v <- certificate(fit)$values
c(bound = v$bound, true_error = abs(fit$values - sqrt(5)))
#>        bound   true_error 
#> 3.863956e-05 1.493015e-09
```

[`certificate()`](https://gillescolling.com/linop/reference/certificate.md)
is the reader. Every verb in the package returns one somewhere in its
result, and this reaches it wherever that is.

That bound is two terms in quadrature and a floor. The truncation term
is two entries of the computed eigenvector and nothing else, and it
falls at the closed- form rate:

``` r

rho <- decay_rate(H, sqrt(5))
t(vapply(c(10, 20, 30, 40), function(n) {
  cv <- eigs(finite_section(H, n = n), k = 1,
             which = "largest_algebraic")$certificate$values
  c(n = n, truncation = cv$truncation, over_rho_n = cv$truncation / rho^n)
}, numeric(3)))
#>       n   truncation over_rho_n
#> [1,] 10 4.754377e-03  0.5847498
#> [2,] 20 3.863956e-05  0.5845006
#> [3,] 30 3.141608e-07  0.5844955
#> [4,] 40 2.566185e-09  0.5872101
```

The last column is constant, which says the rate is `rho` with a
computable constant rather than merely `O(rho^n)`. Push far enough and
it stops being constant, because the tail has fallen below what the
finite eigensolve resolves and the term is reading noise rather than
truncation:

``` r

deep <- eigs(finite_section(H, n = 80), k = 1, which = "largest_algebraic")
c(measured = deep$certificate$values$truncation, analytic = 0.5845 * rho^80)
#>     measured     analytic 
#> 1.392069e-11 1.116280e-17
```

The bound still holds, because it is an identity about the vector
actually computed rather than a statement about the exact eigenvector.
Where the term has gone all the way down to the arithmetic floor the row
says so:

``` r

sharp <- eigs(finite_section(linop_jacobi(diagonal = 4), n = 40), k = 1,
              which = "largest_algebraic")
ch <- sharp$certificate$checks
cat(ch$detail[ch$check == "truncation bound"])
#> worst |b| sqrt(u_{-n}^2 + u_n^2) / ||u|| = 1.352e-11 against 6.000e-10 at n = 40; decay rate 0.236068
```

## Letting it choose the width

`n = 20` above was a guess.
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) on the
operator itself does not need one:

``` r

chosen <- eigs(H, k = 1, which = "largest_algebraic")
chosen
#> <primme on a finite section> 1 of 1 requested eigenpairs converged, which = 'largest_algebraic'
#>   55 iterations of at most 1980, 6 restarts
#>   section n = 49, from 2 widths
#>   values: 2.236068
#>   certificate: pass
```

No width is available before the solve: the tail falls at
`decay_rate(H, lambda)`, and `lambda` is what the solve produces. So the
width is measured rather than derived, and one measurement is enough.
Outside the window the eigenvector is exactly geometric, so the
truncation term at one width gives it at every other,

`eta(n') = eta(n) rho^(n' - n)`,

and the width that brings it to a target is a division:

`n' = n + log(target / eta(n)) / log(rho)`

The widths it went through are on the result, and there are two of them:
one to measure and one to confirm.

``` r

chosen$widths
#> [1] 41 49
```

That holds at every rate. A potential nearer the band edge decays more
slowly and wants a wider section, and the number of widths stays the
same:

``` r

t(vapply(c(0.2, 0.3, 0.5, 1, 2), function(v) {
  Hv <- linop_jacobi(diagonal = v)
  f <- eigs(Hv, k = 1, which = "largest_algebraic")
  c(v = v, rho = decay_rate(Hv, f$values), n = f$n, widths = length(f$widths))
}, numeric(4)))
#>        v       rho   n widths
#> [1,] 0.2 0.9049876 219      2
#> [2,] 0.3 0.8611874 150      2
#> [3,] 0.5 0.7807764  93      2
#> [4,] 1.0 0.6180340  49      2
#> [5,] 2.0 0.4142136  41      1
```

The last row took one, because the first section was already wide enough
at `v = 2`. The rows above it are what the closed form buys: from a
first width of 41 a doubling scheme reaches 82 and then 164, and
`v = 0.3` wants 150.

What stopped the search is on the certificate, where every other
dispatch in the package is recorded:

``` r

certificate(chosen)
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

Three things stop it. The tail can reach the target, which is the line
above. It can stop falling with `n`, which is the plateau: past the
level at which the finite eigensolve stores the tail at all, `eta`
measures the eigenvector’s own error and no width improves it. And
`n_max` is a budget, reported as the width it would have taken:

``` r

capped <- eigs(linop_jacobi(diagonal = 0.3), k = 1, which = "largest_algebraic",
               section = list(n_max = 60))
capped$certificate$dispatch$reason
#> [1] "operator on l^2(Z); width 41 -> 60, stopped because n_max = 60"
ch <- certificate(capped)$checks
ch[ch$check %in% c("finite residual", "truncation bound"), c("check", "status")]
#>              check status
#> 2  finite residual   pass
#> 4 truncation bound   fail
```

The eigensolve was fine and the section was short, and the two rows say
so separately, because they answer to different knobs.

## The subspace is the knob to watch

A section’s spectrum is a discretisation of the essential spectrum
together with whatever separated from it. Enlarging `n` densifies that
cluster rather than adding isolated eigenvalues, so a value close to the
band edge has to be told apart from an ever tighter crowd. A restarted
run separates it in more rounds the smaller its subspace, and
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) widens its
default subspace for a section because of it:

``` r

W <- finite_section(linop_jacobi(diagonal = 0.1), n = 400)
narrow <- eigs(W, k = 1, which = "largest_algebraic", ncv = 21, maxit = 2000)
wide   <- eigs(W, k = 1, which = "largest_algebraic", maxit = 2000)

c(narrow = narrow$values, wide = wide$values, exact = sqrt(0.1^2 + 4))
#>   narrow     wide    exact 
#> 2.002498 2.002498 2.002498
```

Both runs reach the value. The wider subspace costs more storage and
fewer applies of the operator:

``` r

c(narrow = narrow$iterations, wide = wide$iterations)
#> narrow   wide 
#>    404    404
```

See
[`vignette("spectral")`](https://gillescolling.com/linop/articles/spectral.md)
for `ncv` on operators generally.

## When there is nothing to find

A potential of zero leaves the free operator, which has no eigenvalues
at all. Every finite section still produces numbers, and they look like
eigenvalues:

``` r

free <- eigs(finite_section(linop_jacobi(diagonal = 0), n = 40), k = 1,
             which = "largest_algebraic")
free$values
#> [1] 1.998532
free$certificate$values$residual
#> [1] 1.89025e-10
```

The value sits just inside `[-2, 2]` and approaches the edge from below
as `n` grows. The eigensolve converged cleanly and has nothing to
report. What does is the isolation row, which tests the value against
the essential spectrum:

``` r

free$certificate
#> solved by primme (hermitian operator, PRIMME with implicit restarting (highest priority))
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         construction     deterministic_bound   1
#> finite residual                pass         computation      identity              1
#> orthogonality                  pass         computation      identity              1
#> truncation bound               fail         theorem          deterministic_bound   1
#> isolation                      fail         theorem          deterministic_bound   1
#> enclosure                      not_checked  theorem          estimate              -
#> convergence                    fail         computation      identity              1
#> forward error                  fail         theorem          deterministic_bound   1
#> ------------------------------------------------------------------------------------------
#> overall                        fail   no deterministic bound on: enclosure
#> 
#> failures:
#>   truncation bound: worst |b| sqrt(u_{-n}^2 + u_n^2) / ||u|| = 8.460e-03 against 2.000e-10 at n = 40; decay rate undefined inside the band
#>   isolation: sigma_ess(H) = [-2, 2] by Weyl; the closest value clears it by -1.468e-03 against the error bound 8.460e-03; 1 of 1 values have not separated from the continuum
#>   convergence: 0 of 1 requested pairs converged in 131 of at most 1620 iterations; the eigensolve met its target and the section did not; the truncation term is 8.460e-03 at n = 40
#>   forward error: min_j |q - lambda_j(H)| <= 8.460e-03 = sqrt(1.890e-10^2 + 8.460e-03^2) + 1.776e-15, against 2.000e-10; the eigenvectors of H are not covered
```

This is the case the exact essential spectrum exists for. A value that
has not cleared the band is a discretisation of continuous spectrum, and
the certificate refuses rather than issuing a number that would read
like an eigenvalue.

The width search meets the same wall. With nothing outside the band
there is no rate to predict with, so it doubles until its budget and
reports that:

``` r

hunt <- eigs(linop_jacobi(diagonal = 0), k = 1, which = "largest_algebraic",
             section = list(n_max = 200))
hunt$widths
#> [1]  41  82 164 200
hunt$certificate$dispatch$reason
#> [1] "operator on l^2(Z); width 41 -> 82 -> 164 -> 200, stopped because n_max = 200"
```

## Where to next

- [`vignette("certificates")`](https://gillescolling.com/linop/articles/certificates.md)
  covers the fifth certificate shape beside the other four.

- [`vignette("spectral")`](https://gillescolling.com/linop/articles/spectral.md)
  covers [`eigs()`](https://gillescolling.com/linop/reference/eigs.md)
  and [`svds()`](https://gillescolling.com/linop/reference/svds.md) on
  ordinary operators.

- [`vignette("expressions")`](https://gillescolling.com/linop/articles/expressions.md)
  covers the algebra these operators compose under.
