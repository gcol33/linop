## The conversion graph's arithmetic: function-space operators of one kind are
## held as coefficients, so + and scalar * merge them into one operator, whose
## order (hence range space) is the larger of the terms', and whose spectrum and
## solutions are those of the merged problem.

dirichlet <- function() list(at = c(-1, 1), derivative = c(0, 0), value = c(0, 0))

test_that("a sum of Hermite terms is the operator, with its spectrum", {
  kinetic <- linop_space("hermite", derivative = c(0, 0, -1), multiplier = 0)
  potential <- linop_space("hermite", derivative = 0, multiplier = c(0, 0, 1))
  H <- kinetic + 4 * potential
  expect_equal(H$node, "basis")
  expect_equal(H$args$v, c(0, 0, 4))
  f <- eigs(H, k = 3, which = "smallest_algebraic")
  expect_equal(f$values, 2 * (2 * (0:2) + 1), tolerance = 1e-10)
})

test_that("Fourier terms of different degrees align at frequency zero", {
  F1 <- linop_space("fourier")
  F2 <- linop_space("fourier", derivative = 0, multiplier = c(0.5, 0, 0.5))
  F3 <- linop_space("fourier", derivative = 0, multiplier = c(-1, 0, 0, 0, -1))
  S <- F1 + F2 + 2 * F3
  expect_equal(S$args$vhat, c(-2, 0.5, 0, 0.5, -2))
  expect_equal(S$args$coef, c(0, 0, -1))
})

test_that("a Chebyshev sum takes its range from the highest order and solves", {
  ## u'' + u' with u(+-1) = 0: the first-order term is built without conditions,
  ## and the sum carries the problem's
  second <- linop_space("chebyshev", derivative = c(0, 0, 1), boundary = dirichlet())
  first <- linop_space("chebyshev", derivative = c(0, 1))
  expect_equal(first$args$range, "ultraspherical(1)")
  L <- second + first
  expect_equal(L$args$range, "ultraspherical(2)")
  expect_equal(L$args$coef, c(0, 1, 1))
  ## u = sin(pi x) gives f = -pi^2 sin(pi x) + pi cos(pi x); compare values
  k <- 0:40
  s <- ifelse(k %% 2 == 1, 2 * (-1)^((k - 1) / 2) * besselJ(pi, k), 0)
  c0 <- ifelse(k %% 2 == 0, 2 * (-1)^(k / 2) * besselJ(pi, k), 0); c0[1] <- besselJ(pi, 0)
  u <- solve(L, -pi^2 * s + pi * c0)
  xs <- seq(-1, 1, length.out = 21)
  val <- vapply(xs, function(t) sum(u * cos((seq_along(u) - 1) * acos(t))), numeric(1))
  expect_lt(max(abs(val - sin(pi * xs))), 1e-12)
})

test_that("what a sum refuses", {
  dir2 <- list(at = c(-1, 1), derivative = c(0, 1), value = c(0, 0))
  A <- linop_space("chebyshev", derivative = c(0, 0, 1), boundary = dirichlet())
  B <- linop_space("chebyshev", derivative = c(0, 0, 1), boundary = dir2)
  expect_error(A + B, "different boundary conditions")
  expect_error(linop_space("hermite") + linop_space("laguerre", derivative = 1), "one basis")
  ## different kinds of space do not merge; the sum is an ordinary node with no
  ## truncation rule
  mixed <- linop_space("fourier") + linop_space("hermite")
  expect_equal(mixed$node, "sum")
})
