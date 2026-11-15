## Plan sections 8.4 to 8.6 for the Chebyshev operator: its domain, its adjoint
## with the adjoint domain, and Galerkin on its weak form. The adjoint is checked
## by the identity that defines it, int (L u) v = int u (L* v) for u in the domain
## of L and v in the domain of L*, on random polynomials in each domain.

dirichlet <- function() list(at = c(-1, 1), derivative = c(0, 0), value = c(0, 0))
gl <- linop:::gauss_legendre(80)
cheb_at <- function(cf, x) vapply(x, function(t) sum(cf * cos((seq_along(cf) - 1) * acos(t))), numeric(1))
cheb_d <- function(cf) {
  n <- length(cf)
  if (n <= 1L) return(0)
  d <- numeric(n + 1L)
  for (k in (n - 1L):1L) d[k] <- d[k + 2L] + 2 * k * cf[k + 1L]
  d[1L] <- d[1L] / 2
  d[seq_len(n - 1L)]
}
apply_strong <- function(P, cf) {
  out <- 0; dc <- cf
  for (m in seq_along(P$args$coef) - 1L) {
    if (m > 0L) dc <- cheb_d(dc)
    out <- out + P$args$coef[m + 1L] * cheb_at(dc, gl$x)
  }
  out + cheb_at(linop:::cheb_mult(P$args$a, cf), gl$x)
}
in_domain <- function(P, N, seed) {
  set.seed(seed)
  cf <- rnorm(N)
  bc <- P$args$boundary
  R <- t(vapply(seq_along(bc$at), function(i) linop:::bc_row(bc$at[i], bc$weights[[i]], N), numeric(N)))
  as.vector(cf - t(R) %*% solve(R %*% t(R), R %*% cf))
}
lagrange_gap <- function(P, seed) {
  Pa <- adjoint(P)
  u <- in_domain(P, 14, seed); v <- in_domain(Pa, 14, seed + 100)
  lhs <- sum(gl$w * apply_strong(P, u) * cheb_at(v, gl$x))
  rhs <- sum(gl$w * cheb_at(u, gl$x) * apply_strong(Pa, v))
  abs(lhs - rhs) / max(abs(lhs), 1)
}

test_that("the adjoint satisfies its defining identity, Robin conditions included", {
  ops <- list(
    linop_space("chebyshev", derivative = c(0, 0, -1), boundary = dirichlet()),
    linop_space("chebyshev", derivative = c(0, 1, 1),
                boundary = list(at = c(-1, 1), derivative = c(1, 0), value = c(0, 0))),
    linop_space("chebyshev", derivative = c(2, -1, 1), multiplier = c(0.5, 0.2, 0.3),
                boundary = list(at = c(-1, -1), derivative = c(0, 1), value = c(0, 0))),
    linop_space("chebyshev", derivative = c(0, 0, 0, 1),
                boundary = list(at = c(-1, 1, 1), derivative = c(0, 0, 1), value = c(0, 0, 0))))
  for (i in seq_along(ops)) for (seed in 1:5) {
    expect_lt(lagrange_gap(ops[[i]], seed), 1e-10, label = sprintf("operator %d seed %d", i, seed))
  }
  ## u'' + u' with u'(-1) = 0: the adjoint v'' - v' has the mixed condition
  ## -v(-1) + v'(-1) = 0 there
  Ma <- adjoint(ops[[2]])
  expect_equal(Ma$args$coef, c(0, -1, 1))
  w <- Ma$args$boundary$weights[[which(Ma$args$boundary$at == -1)]]
  expect_equal(w / w[1L], c(1, -1), tolerance = 1e-12)
})

test_that("self-adjointness in L2 is read off the coefficients and the conditions", {
  expect_true(linop_space("chebyshev", derivative = c(0, 0, -1), boundary = dirichlet())$args$self_adjoint)
  neumann <- list(at = c(-1, 1), derivative = c(1, 1), value = c(0, 0))
  expect_true(linop_space("chebyshev", derivative = c(1, 0, -1), multiplier = c(0, 0, 1),
                          boundary = neumann)$args$self_adjoint)
  expect_false(linop_space("chebyshev", derivative = c(0, 1, -1), boundary = dirichlet())$args$self_adjoint)
  ## two conditions at one end: an initial value problem is not self-adjoint
  ivp <- list(at = c(-1, -1), derivative = c(0, 1), value = c(0, 0))
  expect_false(linop_space("chebyshev", derivative = c(0, 0, -1), boundary = ivp)$args$self_adjoint)
  expect_true(is.na(linop_space("chebyshev", derivative = c(0, 0, -1))$args$self_adjoint))
})

test_that("the adjoint of the adjoint has the original domain", {
  M <- linop_space("chebyshev", derivative = c(0, 1, 1),
                   boundary = list(at = c(-1, 1), derivative = c(1, 0), value = c(0, 0)))
  M2 <- adjoint(adjoint(M))
  expect_equal(M2$args$coef, M$args$coef)
  for (e in c(-1, 1)) {
    expect_true(linop:::same_rowspace(linop:::cheb_endpoint_rows(M$args$boundary, e, 2L),
                                      linop:::cheb_endpoint_rows(M2$args$boundary, e, 2L)))
  }
  expect_error(adjoint(linop_space("chebyshev", derivative = c(0, 0, 1))), "needs its domain")
})

test_that("adjoints on the circle and the line are the coefficients'", {
  D <- linop_space("fourier", derivative = c(0, 1))
  expect_equal(adjoint(D)$args$coef, c(0, -1))
  V <- linop_space("fourier", derivative = c(0, 0, -1), multiplier = c(1i, 0, 1, 0, -1i))
  expect_equal(adjoint(V)$args$vhat, Conj(rev(V$args$vhat)))
  H <- linop_space("hermite", derivative = c(0, 1, -1), multiplier = c(0, 0, 1))
  expect_equal(adjoint(H)$args$coef, c(0, -1, -1))
})

test_that("Galerkin on the weak form gives upper bounds that agree with the ultraspherical method", {
  L <- linop_space("chebyshev", derivative = c(0, 0, -1), boundary = dirichlet())
  g <- eigs(L, k = 4, which = "smallest", method = "galerkin")
  expect_equal(g$values, ((1:4) * pi / 2)^2, tolerance = 1e-12)
  expect_equal(cert_status(g$certificate, "upper bound"), "pass")
  ## the lifted vector vanishes at both ends and is the sine mode
  u <- g$vectors[, 1L]
  expect_lt(abs(sum(u)), 1e-12)
  expect_lt(abs(sum(u * (-1)^(seq_along(u) - 1))), 1e-12)
  ## with a potential, against the ultraspherical pencil, an independent method
  W <- linop_space("chebyshev", derivative = c(0, 0, -1), multiplier = c(5, 0, 5),
                   boundary = dirichlet())
  gw <- eigs(W, k = 3, which = "smallest", method = "galerkin")
  uw <- eigs(W, k = 3, which = "smallest")
  expect_equal(gw$values, Re(uw$values), tolerance = 1e-10)
})

test_that("Galerkin refuses what its form does not cover", {
  first <- linop_space("chebyshev", derivative = c(0, 1, -1), boundary = dirichlet())
  expect_error(eigs(first, k = 1, which = "smallest", method = "galerkin"), "no first-order term")
  neu <- linop_space("chebyshev", derivative = c(0, 0, -1),
                     boundary = list(at = c(-1, 1), derivative = c(1, 1), value = c(0, 0)))
  expect_error(eigs(neu, k = 1, which = "smallest", method = "galerkin"), "homogeneous Dirichlet")
})

test_that("print shows the domain", {
  M <- linop_space("chebyshev", derivative = c(0, 1, 1),
                   boundary = list(at = c(-1, 1), derivative = c(1, 0), value = c(0, 0)))
  out <- capture.output(print(adjoint(M)))
  expect_true(any(grepl("u(-1) - u'(-1) = 0", out, fixed = TRUE)))
  expect_true(any(grepl("self-adjoint in L2",
                        capture.output(print(linop_space("chebyshev", derivative = c(0, 0, -1),
                                                         boundary = dirichlet()))))))
})
