# An operator on a function space

The one constructor for operators on functions, held as operators on
their coefficients in the space's basis. \`kind = "fourier"\` is the
circle: a function is its Fourier coefficients, a sequence on the
integers, and the operator is a constant-coefficient differential
operator plus multiplication by a trigonometric polynomial,

## Usage

``` r
linop_space(
  kind = "fourier",
  derivative = c(0, 0, -1),
  multiplier = 0,
  boundary = NULL
)
```

## Arguments

- kind:

  The function space: \`"fourier"\`, the circle; \`"chebyshev"\`, the
  interval \`\[-1, 1\]\`; \`"hermite"\`, the line; \`"laguerre"\`, the
  half line.

- derivative:

  Coefficients \`c_0, c_1, ..., c_p\` of \`sum c_m (d/dx)^m\`. The
  default, \`c(0, 0, -1)\`, is \`-d^2/dx^2\`.

- multiplier:

  The coefficients of the multiplying function: Fourier coefficients
  from frequency \`-d\` to \`d\` (an odd number of them) on the circle,
  Chebyshev coefficients \`a_0, a_1, ...\` on the interval, monomial
  coefficients \`v_0, v_1, ...\` of \`V(x)\` on the line and half line.
  A single number is a constant.

- boundary:

  For \`"chebyshev"\`, a list with \`at\` (each \`-1\` or \`1\`),
  \`value\`, and either \`derivative\` (one per condition) or
  \`weights\` (a list of weight vectors over \`u, u', ...\`, for mixed
  conditions), one condition per order of the operator. A term of a sum
  may omit it; the problem needs it.

## Value

A \`linop\` of dimension \`Inf x Inf\`.

## Details

\`H u = sum_m derivative\[m + 1\] (d/dx)^m u + V u\`, with \`V(x) =
sum_j multiplier\[j + d + 1\] exp(i j x)\`.

In coefficients that is a banded operator on \`l^2(Z)\` with an
unbounded diagonal, so it has infinite dimensions: \[finite_section()\]
truncates it to the frequencies \`-n\` to \`n\`, and \[eigs()\] runs on
it directly and chooses the truncation. With a real symbol (the
derivative coefficients times \`i^m\` real) and a real \`V\` (the
multiplier conjugate-symmetric) it is self-adjoint by construction, and
of order at least one its spectrum is discrete, so the certificate's
residual bound is a distance to an eigenvalue of the operator on the
circle itself.

\`kind = "chebyshev"\` is the interval \`\[-1, 1\]\`: a function is its
Chebyshev coefficients, the operator \`sum_m c_m u^(m) + a(x) u\` is
held in the ultraspherical method (Olver and Townsend), mapping
Chebyshev coefficients to those of the ultraspherical basis \`C^(p)\`
where it is banded, and \`boundary\` gives one condition \`u^(d)(+-1) =
value\` per order. \[solve()\] takes the Chebyshev coefficients of the
right-hand side and returns those of the solution, and \[eigs()\] with
\`which = "smallest"\` the eigenvalues of smallest magnitude; both
choose the number of coefficients and certify the exact residual on the
untruncated system. The boundary conditions define the operator's
domain: \[adjoint()\] returns the adjoint operator with the adjoint
domain, by integration by parts (mixed conditions where the boundary
form requires them), and the operator records whether it is self-adjoint
in L2. For \`-c u” + (c_0 + a) u\` with homogeneous Dirichlet
conditions, \`eigs(..., method = "galerkin")\` discretises the weak form
on a basis that vanishes at both ends and returns upper bounds of the
eigenvalues (Courant-Fischer), lifted back to Chebyshev coefficients.

\`kind = "hermite"\` is the line, in the Hermite functions, where \`x\`
and \`d/dx\` are tridiagonal: \`sum_m c_m (d/dx)^m + V(x)\` with \`V\` a
polynomial given by its monomial coefficients. \`kind = "laguerre"\` is
the half line, in the Laguerre functions \`exp(-x/2) L_n(x)\`, where
\`x\` is tridiagonal and the Laguerre operator \`-(x u')' + x u / 4\` is
diagonal: \`c\` times that operator plus a polynomial \`V(x)\`, with
\`derivative = c\`. Both are banded operators on coefficient sequences,
truncated and certified as on the circle.

## Examples

``` r
## -u'' + V u with V = a^2 sin^2 x - a cos x has the ground state
## exp(a cos x) at eigenvalue 0
a <- 2
H <- linop_space("fourier", multiplier = c(-a^2 / 4, -a / 2, a^2 / 2, -a / 2, -a^2 / 4))
fit <- eigs(H, k = 1, which = "smallest_algebraic")
fit$values
#> [1] -4.516328e-17
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
