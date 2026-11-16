# An operator on a function space

The one constructor for operators on functions, held as operators on
their coefficients in the space's basis. `kind = "fourier"` is the
circle: a function is its Fourier coefficients, a sequence on the
integers, and the operator is a constant-coefficient differential
operator plus multiplication by a trigonometric polynomial,
`H u = sum_m derivative[m + 1] (d/dx)^m u + V u`, with
`V(x) = sum_j multiplier[j + d + 1] exp(i j x)`.

## Usage

``` r
linop_space(
  kind = "fourier",
  derivative = c(0, 0, -1),
  multiplier = 0,
  boundary = NULL,
  tol = .Machine$double.eps,
  nmax = 65536L
)
```

## Arguments

- kind:

  The function space: `"fourier"`, the circle; `"chebyshev"`, the
  interval `[-1, 1]`; `"hermite"`, the line; `"laguerre"`, the half
  line.

- derivative:

  Coefficients `c_0, c_1, ..., c_p` of `sum c_m (d/dx)^m`. The default,
  `c(0, 0, -1)`, is `-d^2/dx^2`. For `"chebyshev"` it may also be a list
  of polynomial coefficients `a_0(x), ..., a_p(x)`, each given by its
  Chebyshev coefficients, for `sum a_m(x) u^(m)`. Each entry of that
  list may also be a function of `x`.

- multiplier:

  The coefficients of the multiplying function: Fourier coefficients
  from frequency `-d` to `d` (an odd number of them) on the circle,
  Chebyshev coefficients `a_0, a_1, ...` on the interval, monomial
  coefficients `v_0, v_1, ...` of `V(x)` on the line and half line. A
  single number is a constant. A function of `x` is sampled adaptively
  in the selected space. It must return one finite value per point.
  Scalar functions that error on vector input are evaluated pointwise
  after one detected fallback. Certificate construction may make one
  cached complex vector call to estimate an analytic tail on an ellipse
  or strip. The sampled maximum and analyticity are estimates. Real-only
  functions retain coefficient and off-grid estimates.

- boundary:

  For `"chebyshev"`, a list with `at` (each `-1` or `1`), `value`, and
  either `derivative` (one per condition) or `weights` (a list of weight
  vectors over `u, u', ...`, for mixed conditions), one condition per
  order of the operator. A term of a sum may omit it; the problem needs
  it. `value` may be a function of the endpoint, or a list with a
  function or scalar per endpoint. Functions are approximated on the
  interval first.

- tol:

  Relative tolerance for adaptive function coefficients. Defaults to
  machine epsilon. Chopping is an estimate recorded in the certificate.
  Solve and eigen requests may shorten stored multipliers when the
  estimated sup-norm error fits one tenth of the request and the
  arithmetic floor permits it.

- nmax:

  Maximum sample grid length for each function. An unresolved function
  is refused by argument name. Defaults to `65536`.

## Value

A `linop` of dimension `Inf x Inf`.

## Details

In coefficients that is a banded operator on `l^2(Z)` with an unbounded
diagonal, so it has infinite dimensions:
[`finite_section()`](https://gillescolling.com/linop/reference/finite_section.md)
truncates it to the frequencies `-n` to `n`, and
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) runs on it
directly and chooses the truncation. With a real symbol (the derivative
coefficients times `i^m` real) and a real `V` (the multiplier
conjugate-symmetric) it is self-adjoint by construction, and of order at
least one its spectrum is discrete, so the certificate's residual bound
is a distance to an eigenvalue of the operator on the circle itself.

`kind = "chebyshev"` is the interval `[-1, 1]`: a function is its
Chebyshev coefficients, the operator `sum_m c_m u^(m) + a(x) u` is held
in the ultraspherical method (Olver and Townsend), mapping Chebyshev
coefficients to those of the ultraspherical basis `C^(p)` where it is
banded, and `boundary` gives one condition `u^(d)(+-1) = value` per
order.
[`solve()`](https://gillescolling.com/linop/reference/solve.linop.md)
takes the Chebyshev coefficients of the right-hand side and returns
those of the solution, and
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) with
`which = "smallest"` the eigenvalues of smallest magnitude; both choose
the number of coefficients and certify the exact residual on the
untruncated system. The boundary conditions define the operator's
domain:
[`adjoint()`](https://gillescolling.com/linop/reference/adjoint.md)
returns the adjoint operator with the adjoint domain, by integration by
parts (mixed conditions where the boundary form requires them), and the
operator records whether it is self-adjoint in L2. For a self-adjoint
`-(k u')' + q u` with `k > 0` and homogeneous Dirichlet, Neumann or
Robin conditions, `eigs(..., method = "galerkin")` discretises the weak
form (natural conditions enter it as boundary terms, Dirichlet ones as a
basis that vanishes there) and returns upper bounds of the eigenvalues
(Courant-Fischer), lifted back to Chebyshev coefficients.

`kind = "hermite"` is the line, in the Hermite functions, where `x` and
`d/dx` are tridiagonal: `sum_m c_m (d/dx)^m + V(x)` with `V` a
polynomial given by its monomial coefficients. `kind = "laguerre"` is
the half line, in the Laguerre functions `exp(-x/2) L_n(x)`, where `x`
is tridiagonal and the Laguerre operator `-(x u')' + x u / 4` is
diagonal: `c` times that operator plus a polynomial `V(x)`, with
`derivative = c`. Both are banded operators on coefficient sequences,
truncated and certified as on the circle.

## References

Aurentz, J. L. and Trefethen, L. N. (2017). Chopping a Chebyshev series.
ACM Transactions on Mathematical Software 43(4).
[doi:10.1145/2998442](https://doi.org/10.1145/2998442) .

## Examples

``` r
## -u'' + V u with V = a^2 sin^2 x - a cos x has the ground state
## exp(a cos x) at eigenvalue 0
a <- 2
H <- linop_space("fourier", multiplier = c(-a^2 / 4, -a / 2, a^2 / 2, -a / 2, -a^2 / 4))
fit <- eigs(H, k = 1, which = "smallest_algebraic")
fit$values
#> [1] 3.650369e-17
fit$n
#> [1] 16

## u'' = -pi^2 sin(pi x) on [-1, 1] with u(-1) = u(1) = 0, whose solution is
## sin(pi x); the right-hand side's Chebyshev coefficients are Bessel values
k <- 0:40
f <- -pi^2 * ifelse(k %% 2 == 1, 2 * (-1)^((k - 1) / 2) * besselJ(pi, k), 0)
L <- linop_space("chebyshev", derivative = c(0, 0, 1),
                 boundary = list(at = c(-1, 1), derivative = c(0, 0), value = c(0, 0)))
u <- solve(L, f)
x <- 0.3
sum(u * cos((seq_along(u) - 1) * acos(x))) - sin(pi * x)
#> [1] -4.440892e-16
```
