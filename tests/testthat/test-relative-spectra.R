relative_fixture <- function(v = 0, coupling = 0.4, shift = 0, scale = 1, site = 0L) {
  b <- c(1, coupling)
  linop_jacobi(shift + scale * v, from = site,
               offdiagonal = scale * b[((site - 1L):site %% 2L) + 1L],
               diagonal_limit = shift, offdiagonal_limit = scale * b)
}

## The diagonal Green function of the period-two operator gives
## (lambda^2-(1+b)^2)(lambda^2-(1-b)^2) = v^2 lambda^2.
## Its smaller squared root is the gap eigenvalue, with sign opposite to v.
relative_truth <- function(v, coupling = 0.4) {
  s <- (1 + coupling)^2 + (1 - coupling)^2 + v^2
  p <- (1 - coupling^2)^2
  -sign(v) * sqrt(2 * p / (s + sqrt(s^2 - 4 * p)))
}

test_that("second order retains the cut rows when it squares the infinite operator", {
  H <- relative_fixture()
  F <- finite_section(H, 20)
  J <- as.matrix(F)
  wide <- as.matrix(finite_section(H, 21))
  rectangular <- wide[, 2:42]
  K <- crossprod(rectangular)
  cuts <- diag(K - J %*% J)
  expect_equal(cuts[c(1, 41)], c(0.4, 1)^2, tolerance = 1e-14)
  expect_equal(cuts[2:40], rep(0, 39), tolerance = 1e-14)
  fit <- eigs(F, k = 1, sigma = 0.1, method = "second_order")
  z <- fit$roots
  x <- fit$vectors
  Qx <- (K - 2 * z * J + diag(z^2, 41)) %*% x
  expect_lt(sqrt(sum(Mod(Qx)^2)), 1e-11)
  expect_gt(abs(Im(z)), 0.6)
  ordinary <- eigen(J, symmetric = TRUE, only.values = TRUE)$values
  expect_lt(min(abs(ordinary)), 1e-13)
  expect_gt(min(Mod(fit$all_roots)), 0.6)
})

test_that("the periodic operator has no relative enclosure isolated inside its empty gap", {
  H <- relative_fixture()
  for (n in c(12L, 20L, 41L, 60L)) {
    fit <- eigs(finite_section(H, n), k = 1, sigma = 0.1, method = "second_order")
    z <- fit$all_roots
    lo <- Re(z) - abs(Im(z))
    hi <- Re(z) + abs(Im(z))
    expect_false(any(lo > -0.6 & hi < 0.6))
    expect_false(any(fit$enclosures[, 1L] > -0.6 & fit$enclosures[, 2L] < 0.6))
    expect_equal(cert_status(certificate(fit), "isolation"), "fail")
    expect_equal(cert_status(certificate(fit), "enclosure"), "fail")
    expect_equal(fit$nconv, 0L)
    expect_equal(cert_status(certificate(fit), "gap completeness"), "not_checked")
  }
  fit <- eigs(H, k = 1, sigma = 0.1, method = "second_order")
  expect_lte(length(fit$widths), 3L)
  expect_match(certificate(fit)$dispatch$reason, "stopped shrinking")
  expect_equal(certificate(fit)$overall, "fail")
})

test_that("the defect enclosure recovers independent truth and shrinks with section width", {
  H <- relative_fixture(0.3)
  dense <- eigen(as.matrix(finite_section(H, 600)), symmetric = TRUE,
                 only.values = TRUE)$values
  truth <- relative_truth(0.3)
  reference <- dense[which.min(abs(dense - truth))]
  expect_equal(reference, truth, tolerance = 1e-12)
  expect_equal(truth, -0.583987233542, tolerance = 1e-12)
  widths <- c(20L, 40L, 80L, 120L)
  radii <- numeric(length(widths))
  for (i in seq_along(widths)) {
    fit <- eigs(finite_section(H, widths[i]), k = 1, sigma = -0.58,
                 method = "second_order", tol = 1e-3)
    radii[i] <- certificate(fit)$values$bound
    expect_lte(fit$enclosures[1, 1], reference)
    expect_gte(fit$enclosures[1, 2], reference)
    expect_lte(Re(fit$roots) - abs(Im(fit$roots)), reference)
    expect_gte(Re(fit$roots) + abs(Im(fit$roots)), reference)
    relative <- certificate(fit)$values$relative_enclosures
    expect_lte(relative[1, 1], reference)
    expect_gte(relative[1, 2], reference)
    expect_equal(certificate(fit)$evidence[["second-order inclusion"]]$source, "theorem")
    expect_equal(certificate(fit)$evidence[["second-order inclusion"]]$guarantee,
                 "deterministic_bound")
    expect_equal(cert_status(certificate(fit), "isolation"), "pass")
    if (i > 1L) {
      expect_equal(cert_status(certificate(fit), "enclosure"), "pass")
      expect_equal(cert_status(certificate(fit), "convergence"), "pass")
    }
  }
  expect_true(all(diff(radii) < 0))
  expect_lt(radii[4L] / radii[1L], 1e-3)
  adaptive <- eigs(H, k = 1, sigma = -0.58, method = "second_order", tol = 1e-5)
  expect_equal(cert_status(certificate(adaptive), "enclosure"), "pass")
  expect_equal(cert_status(certificate(adaptive), "convergence"), "pass")
  expect_lte(adaptive$enclosures[1, 1], reference)
  expect_gte(adaptive$enclosures[1, 2], reference)
})

test_that("second-order certificates cover translated and scaled closed-form truth over 20 seeds", {
  for (seed in 1:20) {
    set.seed(seed + 4000L)
    b <- runif(1, 0.3, 0.5)
    v <- runif(1, 0.35, 0.75) * sample(c(-1, 1), 1)
    shift <- runif(1, -0.2, 0.2)
    scale <- runif(1, 0.8, 1.3)
    site <- sample(-3:3, 1)
    truth <- shift + scale * relative_truth(v, b)
    H <- relative_fixture(v, b, shift, scale, site)
    fit <- eigs(finite_section(H, 80), k = 1, sigma = truth,
                 method = "second_order", tol = 1e-3, seed = seed)
    expect_lte(fit$enclosures[1, 1], truth)
    expect_gte(fit$enclosures[1, 2], truth)
    expect_equal(cert_status(certificate(fit), "isolation"), "pass")
    expect_equal(cert_status(certificate(fit), "enclosure"), "pass")
    ev <- certificate(fit)$evidence$enclosure
    expect_true(linop:::evidence_satisfies(ev,
                  linop:::requirement(guarantees = c("identity", "deterministic_bound"))))
    expect_equal(ev$source, "theorem")
    expect_equal(ev$guarantee, "deterministic_bound")
  }
})

test_that("barriers separate genuine and spurious gap points over 20 seeds and three gammas", {
  genuine_ratio <- spurious_ratio <- numeric(0)
  for (seed in 1:20) {
    set.seed(seed + 5000L)
    b <- runif(1, 0.3, 0.5)
    v <- runif(1, 0.4, 0.75) * sample(c(-1, 1), 1)
    shift <- runif(1, -0.2, 0.2)
    scale <- runif(1, 0.8, 1.3)
    site <- sample(-3:3, 1)
    truth <- shift + scale * relative_truth(v, b)
    H <- relative_fixture(v, b, shift, scale, site)
    free <- relative_fixture(0, b, shift, scale, site)
    for (gamma in c(0.002, 0.01, 0.05) * scale) {
      opts <- list(gamma = gamma, n = 50L)
      good <- eigs(finite_section(H, 100), k = 1, sigma = truth,
                    method = "dissipative", barrier = opts, tol = 1e-3, seed = seed)
      bad <- eigs(finite_section(free, 100), k = 1, sigma = shift,
                   method = "dissipative", barrier = opts, tol = 1e-3, seed = seed)
      g <- certificate(good)$values
      f <- certificate(bad)$values
      genuine_ratio <- c(genuine_ratio, g$barrier_ratio)
      spurious_ratio <- c(spurious_ratio, f$barrier_ratio)
      expect_equal(g$barrier_ratio, g$projection_mass, tolerance = 1e-10)
      expect_lte(good$enclosures[1, 1], truth)
      expect_gte(good$enclosures[1, 2], truth)
      expect_equal(cert_status(certificate(good), "barrier separation"), "pass")
      expect_equal(cert_status(certificate(good), "enclosure"), "pass")
      expect_equal(cert_status(certificate(good), "isolation"), "pass")
      expect_equal(cert_status(certificate(bad), "barrier separation"), "fail")
      expect_equal(cert_status(certificate(bad), "isolation"), "fail")
      expect_equal(bad$nconv, 0L)
      expect_equal(certificate(good)$checks$guarantee[
                     certificate(good)$checks$check == "barrier separation"], "estimate")
      expect_true("barrier separation" %in% certificate(good)$without_deterministic_bound)
    }
  }
  expect_gt(min(genuine_ratio), 0.99)
  expect_lt(max(abs(spurious_ratio)), 1e-8)
})

test_that("gamma and projection width change the perturbation at the same section budget", {
  F <- finite_section(relative_fixture(0.3), 100)
  small <- eigs(F, k = 1, sigma = -0.58, method = "dissipative",
                barrier = list(gamma = 0.002, n = 40), tol = 1e-3)
  large <- eigs(F, k = 1, sigma = -0.58, method = "dissipative",
                barrier = list(gamma = 0.01, n = 40), tol = 1e-3)
  expect_gt(Im(large$roots) / Im(small$roots), 4.9)
  expect_lt(Im(large$roots) / Im(small$roots), 5.1)
  phase <- sum(Conj(small$vectors) * large$vectors)
  expect_gt(max(Mod(small$vectors - large$vectors * Conj(phase) / Mod(phase))), 1e-6)
  narrow <- eigs(F, k = 1, sigma = -0.58, method = "dissipative",
                 barrier = list(gamma = 0.01, n = 8), tol = 1e-3)
  expect_gt(certificate(large)$values$barrier_ratio - certificate(narrow)$values$barrier_ratio, 0.1)
  expect_lt(certificate(large)$values$bound / certificate(narrow)$values$bound, 0.1)
  expect_equal(cert_status(certificate(large), "enclosure"), "pass")
  expect_equal(cert_status(certificate(narrow), "enclosure"), "fail")
})

test_that("a finite barrier diagnostic alone cannot certify a false positive", {
  H <- relative_fixture(coupling = 0.9)
  fit <- eigs(finite_section(H, 30), k = 1, sigma = 0,
               method = "dissipative", barrier = list(gamma = 0.01, n = 29))
  expect_gt(certificate(fit)$values$barrier_ratio, 0.5)
  expect_equal(cert_status(certificate(fit), "barrier separation"), "pass")
  expect_equal(cert_status(certificate(fit), "isolation"), "fail")
  expect_equal(fit$nconv, 0L)
  expect_equal(certificate(fit)$overall, "fail")
})

test_that("gap methods reject unsupported or ineffective controls", {
  H <- relative_fixture(0.3)
  run <- function(...) eigs(H, k = 1, sigma = -0.58, method = "second_order", ...)
  expect_error(run(maxit = 1), "dense section")
  expect_error(run(ncv = 5), "dense section")
  expect_error(run(barrier = list(gamma = 0.1)), "barrier applies")
  expect_error(run(section = list(n_start = 2.5)), "whole")
  expect_error(run(section = list(n_start = 40, n_max = 20)), "first width")
  expect_error(eigs(H, 1, method = "second_order"), "requires a real sigma")
  expect_error(eigs(H, 1, method = "second_order", sigma = 1), "validated gap")
  expect_error(eigs(linop_jacobi(1), 1, method = "second_order", sigma = 0), "validated gap")
  expect_error(eigs(linop(diag(3)), 1, method = "second_order", sigma = 0),
               "defined on the 'jacobi' and 'section' node types")
  expect_error(eigs(H, 1, sigma = -0.58, method = "dissipative",
                     barrier = list(gamma = 0)), "gamma")
  expect_error(eigs(finite_section(H, 20), 1, sigma = -0.58, method = "dissipative",
                     barrier = list(n = 20)), "below section width")
})

test_that("coalescing conjugate roots do not count the same eigenvector twice", {
  H <- relative_fixture(0.3)
  fit <- eigs(finite_section(H, 220), 2, sigma = -0.58,
               method = "second_order", tol = 1e-7)
  truth <- relative_truth(0.3)
  expect_equal(sum(abs(fit$values - truth) < 1e-7), 1L)
  expect_equal(fit$nconv, 1L)
  expect_equal(cert_status(certificate(fit), "convergence"), "fail")
  expect_lte(fit$enclosures[1, 1], truth)
  expect_gte(fit$enclosures[1, 2], truth)
})

test_that("an empty finite candidate set leaves gap completeness unchecked", {
  H <- linop_jacobi(1, diagonal_limit = c(1, -1))
  for (method in c("second_order", "dissipative")) {
    fit <- eigs(H, 1, sigma = 0, method = method, section = list(n_start = 12))
    expect_length(fit$values, 0)
    expect_equal(fit$nconv, 0L)
    expect_equal(cert_status(certificate(fit), "enclosure"), "not_checked")
    expect_equal(cert_status(certificate(fit), "gap completeness"), "not_checked")
    expect_equal(cert_status(certificate(fit), "convergence"), "fail")
  }
})
