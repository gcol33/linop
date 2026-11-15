## The RSpectra engine behind method = "rspectra": recovery against closed forms
## and against eigen(), the refusals the registry makes for it, the budget mapped
## from applies to restarts, and the certificate built on its vectors.

skip_if_not_installed("RSpectra")

nearest_error <- function(values, truth) {
  max(vapply(values, function(v) min(Mod(truth - v)), numeric(1)))
}

test_that("a non-normal operator: convection-diffusion at both ends", {
  n <- 80; mu <- 0.04
  A <- convdiff_1d(n, mu)
  truth <- convdiff_1d_eigenvalues(n, mu)
  fit <- eigs(A, k = 4, which = "largest_real", method = "rspectra")
  expect_equal(fit$method, "rspectra")
  expect_equal(fit$nconv, 4L)
  expect_true(is.double(fit$values))
  expect_equal(fit$values, sort(truth, decreasing = TRUE)[1:4], tolerance = 1e-9)
  R <- linop:::linop_apply(A, fit$vectors, "N") - fit$vectors * rep(fit$values, each = n)
  expect_lt(max(sqrt(colSums(R^2))), 1e-8)
  sm <- eigs(A, k = 3, which = "smallest_real", method = "rspectra")
  expect_equal(sm$values, sort(truth)[1:3], tolerance = 1e-9)
})

test_that("a hermitian operator goes to the symmetric driver at every end", {
  n <- 80
  A <- laplacian_1d(n)
  ev <- laplacian_1d_eigenvalues(n)
  want <- list(largest_algebraic = rev(ev)[1:3], smallest_algebraic = ev[1:3],
               largest = rev(ev)[1:3])
  for (w in names(want)) {
    fit <- eigs(A, k = 3, which = w, method = "rspectra")
    expect_equal(fit$nconv, 3L, info = w)
    expect_lt(max(abs(fit$values - want[[w]])), 1e-9, label = w)
    expect_equal(cert_status(fit$certificate, "backward error"), "pass")
    expect_equal(cert_status(fit$certificate, "forward error"), "pass")
  }
  ## the values are Rayleigh quotients measured on A, not the engine's own
  fit <- eigs(A, k = 3, which = "largest_algebraic", method = "rspectra")
  X <- fit$vectors
  expect_equal(fit$values, colSums(X * linop:::linop_apply(A, X, "N")) / colSums(X^2),
               tolerance = 1e-15)
})

test_that("a real operator's complex eigenvalues, against eigen() over 20 seeds", {
  for (s in 1:20) {
    n <- 60
    M <- rmat(n, n, seed = s) / sqrt(n)
    ev <- eigen(M, only.values = TRUE)$values
    fit <- eigs(linop(M), k = 4, which = "largest", method = "rspectra", seed = s)
    info <- sprintf("seed %d", s)
    expect_lt(nearest_error(fit$values, ev), 1e-9, label = info)
    expect_equal(cert_status(fit$certificate, "backward error"), "pass", info = info)
    expect_false(fit$certificate$overall == "fail", info = info)
    expect_equal(fit$certificate$checks$status[fit$certificate$checks$check == "forward error"],
                 "not_checked", info = info)
  }
})

test_that("auto never chooses it, and the dispatch lists it as eligible", {
  A <- convdiff_1d(60, 0.04)
  for (p in list(NULL, c(primme = -Inf))) {
    old <- options(linop.backend_priority = p)
    fit <- eigs(A, k = 2, which = "largest_real")
    options(old)
    cand <- fit$certificate$dispatch$candidates
    expect_equal(fit$method, "arnoldi")
    expect_equal(cand$status[cand$backend == "rspectra"], "eligible")
  }
})

test_that("a complex operator and a complex shift are refused before any call", {
  ## RSpectra coerces a complex callback to double and returns the spectrum of
  ## Re(A); the refusal is linop's, on the dtype it tracks.
  set.seed(3)
  Z <- matrix(complex(real = rnorm(400), imaginary = rnorm(400)), 20)
  H <- linop(Z + Conj(t(Z)))
  expect_error(eigs(H, k = 2, method = "rspectra"), "does not take a complex operator")
  expect_error(svds(linop(Z), k = 2, method = "rspectra"), "does not take a complex operator")
  expect_error(eigs(convdiff_1d(40, 0.04), k = 2, sigma = 1 + 0.5i, method = "rspectra"),
               "complex shift")
})

test_that("shift-invert runs on eigs()'s own (A - sigma I)^-1", {
  n <- 80; mu <- 0.04
  A <- convdiff_1d(n, mu)
  truth <- convdiff_1d_eigenvalues(n, mu)
  fit <- eigs(A, k = 3, sigma = 1.3, method = "rspectra")
  expect_equal(fit$method, "rspectra shift-invert")
  expect_equal(sort(fit$values), sort(truth[order(abs(truth - 1.3))[1:3]]), tolerance = 1e-10)
  expect_equal(cert_status(fit$certificate, "backward error"), "pass")
})

test_that("the budget in applies maps onto restarts and is never overspent", {
  A <- convdiff_1d(200, 0.04)
  for (maxit in c(23L, 30L, 60L, 100L, 250L)) {
    fit <- eigs(A, k = 3, which = "largest_real", method = "rspectra", maxit = maxit)
    expect_lte(fit$iterations, maxit)
  }
  ## below the first factorisation nothing can stop the run, so it is refused
  expect_error(eigs(A, k = 3, which = "largest_real", method = "rspectra", maxit = 10L),
               "first factorisation spends 23 applies, more than maxit = 10")
  ## a run that converges nothing returns nothing and a fail certificate
  fit <- eigs(A, k = 3, which = "largest_real", method = "rspectra", maxit = 30L)
  expect_length(fit$values, 0L)
  expect_equal(fit$certificate$overall, "fail")
  expect_match(fit$certificate$checks$detail[fit$certificate$checks$check == "convergence"],
               "budget exhausted")
})

test_that("singular triplets, the largest only, and the refusals that follow", {
  sigma <- c(20, 9, 5, seq(2, 0.4, length.out = 17))
  f <- lsq_prescribed(30, 20, sigma, seed = 5L)
  A <- linop(f$A)
  fit <- svds(A, k = 3, method = "rspectra")
  expect_equal(fit$method, "rspectra")
  expect_equal(fit$nconv, 3L)
  expect_lt(max(abs(fit$d - sigma[1:3])), 1e-9)
  expect_equal(fit$d, colSums(fit$u * (f$A %*% fit$v)), tolerance = 1e-15)
  expect_error(svds(A, k = 2, which = "smallest", method = "rspectra"),
               "does not select which = 'SM'")
  expect_error(svds(A, k = 2, v0 = rep(1, 20), method = "rspectra"), "no starting vector")
  for (maxit in c(60L, 100L)) {
    expect_lte(svds(A, k = 3, method = "rspectra", maxit = maxit)$iterations, maxit)
  }
})
