## Phase 6 step 5: the line in Hermite functions and the half line in Laguerre
## functions, where x is tridiagonal and one differential operator is diagonal.
## Closed forms: -u'' + w^2 x^2 u has eigenvalues w (2n + 1); adding b x shifts
## them by -b^2/4; the Laguerre operator T alone has n + 1/2, and T + b x has
## (2n + 1) sqrt(1/4 + b).

smallest <- function(A, k = 3) eigs(A, k = k, which = "smallest_algebraic")

test_that("the harmonic oscillator, scaled and shifted, in Hermite functions", {
  f <- smallest(linop_space("hermite", multiplier = c(0, 0, 1)))
  expect_equal(f$values, c(1, 3, 5), tolerance = 1e-12)
  for (w in c(0.5, 2)) {
    f <- smallest(linop_space("hermite", multiplier = c(0, 0, w^2)))
    expect_equal(f$values, w * (2 * (0:2) + 1), tolerance = 1e-10, label = sprintf("w = %g", w))
    expect_equal(cert_status(f$certificate, "forward error"), "pass")
    expect_equal(cert_status(f$certificate, "truncation bound"), "pass")
  }
  f <- smallest(linop_space("hermite", multiplier = c(0, 1, 1)))
  expect_equal(f$values, 2 * (0:2) + 1 - 1 / 4, tolerance = 1e-10)
})

test_that("the Laguerre operator and its linear perturbation on the half line", {
  f <- smallest(linop_space("laguerre", derivative = 1))
  expect_equal(f$values, (0:2) + 0.5, tolerance = 1e-12)
  for (b in c(0.75, 2)) {
    f <- smallest(linop_space("laguerre", derivative = 1, multiplier = c(0, b)))
    expect_equal(f$values, (2 * (0:2) + 1) * sqrt(1 / 4 + b), tolerance = 1e-10,
                 label = sprintf("b = %g", b))
    expect_equal(cert_status(f$certificate, "forward error"), "pass")
  }
})

test_that("the quartic oscillator agrees with a wide dense section and certifies", {
  A <- linop_space("hermite", multiplier = c(0, 0, 0, 0, 1))
  f <- smallest(A, k = 1)
  ref <- min(eigen(as.matrix(finite_section(A, 300)), symmetric = TRUE, only.values = TRUE)$values)
  expect_equal(f$values, ref, tolerance = 1e-10)
  expect_equal(cert_status(f$certificate, "forward error"), "pass")
  expect_match(f$certificate$dispatch$reason, "l^2(N)", fixed = TRUE)
})

test_that("the section block is exact: the truncation of a wider block", {
  A <- linop_space("hermite", derivative = c(0, 0, -1), multiplier = c(0.5, 0, 1, 0, 0.1))
  small <- as.matrix(finite_section(A, 20))
  big <- as.matrix(finite_section(A, 60))
  expect_lt(max(abs(small - big[1:21, 1:21])), 1e-12)
  ## -u'' + x^2 u is exactly diagonal in Hermite functions
  ho <- as.matrix(finite_section(linop_space("hermite", multiplier = c(0, 0, 1)), 10))
  expect_lt(max(abs(ho - diag(2 * (0:10) + 1))), 1e-12)
})

test_that("the edge term is the residual of the zero-extended vector on a wider section", {
  A <- linop_space("hermite", multiplier = c(0, 0, 0, 0, 1))
  for (n in c(12L, 24L)) {
    S <- as.matrix(finite_section(A, n))
    e <- eigen(S, symmetric = TRUE); q <- min(e$values); u <- e$vectors[, which.min(e$values)]
    W <- as.matrix(finite_section(A, n + 10L))
    ut <- c(u, numeric(10))
    full <- sqrt(sum((W %*% ut - q * ut)^2))
    res <- sqrt(sum((S %*% u - q * u)^2))
    edge <- linop:::basis_edge(A, matrix(u, ncol = 1L))
    expect_equal(full, sqrt(res^2 + edge^2), tolerance = 1e-10)
  }
})

test_that("what is not confining or not self-adjoint is said so", {
  free <- linop_space("hermite")
  expect_false(free$args$discrete)
  odd <- linop_space("hermite", derivative = c(0, 1))
  expect_false(isTRUE(capv(odd, "hermitian")))
  fit <- eigs(finite_section(odd, 20), k = 1, which = "smallest")
  expect_equal(cert_status(fit$certificate, "forward error"), "not_checked")
  expect_error(linop_space("laguerre", derivative = c(0, 1)), "single coefficient")
  expect_error(finite_section(linop_space("hermite", multiplier = c(0, 0, 0, 1)), 2), "narrower than the band")
})
