## Section 7.2, eigs() on an operator that is not established hermitian: Arnoldi
## with Krylov-Schur restarting. Recovery against spectra known by construction,
## the certificate's rows for a non-normal operator, and the Schur kernel the
## restart rests on.

## A real matrix whose spectrum is the given conjugate pairs, by construction:
## 2 x 2 rotation blocks under a random similarity.
real_with_pairs <- function(lam, seed) {
  p <- length(lam)
  D <- matrix(0, 2 * p, 2 * p)
  for (i in seq_len(p)) {
    j <- 2L * i - 1L
    D[j:(j + 1L), j:(j + 1L)] <- matrix(c(Re(lam[i]), -Im(lam[i]), Im(lam[i]), Re(lam[i])), 2L)
  }
  S <- rmat(2 * p, 2 * p, seed = seed)
  S %*% D %*% solve(S, diag(2 * p))
}

## A complex non-normal matrix in Schur form under a random unitary: its
## eigenvalues are the diagonal of the triangular factor.
complex_schur_fixture <- function(ev, offdiag, seed) {
  n <- length(ev)
  set.seed(seed)
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n) + 1i * matrix(stats::rnorm(n * n), n)))
  U <- diag(ev)
  up <- upper.tri(U)
  U[up] <- offdiag * complex(real = stats::rnorm(sum(up)), imaginary = stats::rnorm(sum(up)))
  Q %*% U %*% Conj(t(Q))
}

nearest_error <- function(values, truth) {
  max(vapply(values, function(v) min(Mod(truth - v)), numeric(1)))
}

## ----------------------------------------------------------- recovery --

test_that("a non-normal operator with a real spectrum: convection-diffusion", {
  ## rho^n = ((1 + mu)/(1 - mu))^(n/2) = 24, inside the fixture's budget.
  n <- 80; mu <- 0.04
  A <- convdiff_1d(n, mu)
  truth <- convdiff_1d_eigenvalues(n, mu)
  fit <- eigs(A, k = 4, which = "largest_real")
  expect_equal(fit$method, "arnoldi")
  expect_equal(fit$nconv, 4L)
  expect_true(is.double(fit$values))
  expect_true(is.double(fit$vectors))
  expect_equal(fit$values, sort(truth, decreasing = TRUE)[1:4], tolerance = 1e-10)
  ## A S - S Lambda, never the sorted spectrum alone
  R <- linop:::linop_apply(A, fit$vectors, "N") -
    fit$vectors * rep(fit$values, each = n)
  expect_lt(max(sqrt(colSums(R^2))), 1e-9)

  sm <- eigs(A, k = 3, which = "smallest_real")
  expect_equal(sm$values, sort(truth)[1:3], tolerance = 1e-10)
})

test_that("a real operator's conjugate pairs come back together and conjugate", {
  lam <- complex(real = seq(1, 3, length.out = 20), imaginary = seq(0.5, 2, length.out = 20))
  M <- real_with_pairs(lam, seed = 3L)
  truth <- c(lam, Conj(lam))
  for (w in c("largest", "largest_real", "largest_imaginary", "smallest_real")) {
    fit <- eigs(linop(M), k = 4, which = w)
    expect_equal(fit$nconv, 4L, info = w)
    expect_lt(nearest_error(fit$values, truth), 1e-9)
    ## two pairs: each value's conjugate is among the four, to rounding
    gap <- vapply(fit$values, function(v) min(Mod(fit$values - Conj(v))), numeric(1))
    expect_lt(max(gap), 1e-9)
    expect_true(all(Im(fit$values) != 0), info = w)
  }
  ## the top of the spectrum by magnitude and by imaginary part is 3 +- 2i
  top <- eigs(linop(M), k = 2, which = "largest_imaginary")$values
  expect_equal(sort(Im(top)), c(-2, 2), tolerance = 1e-9)
  expect_equal(Re(top), c(3, 3), tolerance = 1e-9)
})

test_that("a split pair is completed when k ends inside it", {
  lam <- complex(real = c(5, 2, 1), imaginary = c(1, 0.5, 0.25))
  M <- real_with_pairs(lam, seed = 8L)
  fit <- eigs(linop(M), k = 1, which = "largest")
  expect_length(fit$values, 1L)
  expect_equal(Mod(fit$values), Mod(5 + 1i), tolerance = 1e-10)
  expect_equal(fit$nconv, 1L)
})

test_that("a complex non-normal operator", {
  set.seed(40)
  ev <- complex(real = stats::rnorm(70), imaginary = stats::rnorm(70))
  ev[1:3] <- c(5 + 1i, 4.5 - 2i, 4 + 3i)
  M <- complex_schur_fixture(ev, offdiag = 0.3, seed = 4L)
  fit <- eigs(linop(M), k = 3, which = "largest_real")
  expect_equal(fit$values, ev[1:3], tolerance = 1e-10)
  expect_true(is.complex(fit$vectors))
})

test_that("an operator declaring nothing runs Arnoldi, and the dispatch says why", {
  n <- 60
  L <- laplacian_1d(n)
  undeclared <- linop(function(X) linop:::linop_apply(L, X, "N"), dim = c(n, n))
  fit <- eigs(undeclared, k = 3, which = "largest_algebraic")
  expect_equal(fit$method, "arnoldi")
  expect_equal(fit$values, rev(sort(laplacian_1d_eigenvalues(n)))[1:3], tolerance = 1e-10)
  expect_true(is.double(fit$values))
  cand <- fit$certificate$dispatch$candidates
  expect_equal(cand$status[cand$backend == "lanczos"], "rejected")
  expect_match(cand$reason[cand$backend == "lanczos"], "unknown is not false")
  expect_false(grepl("\n", cand$reason[cand$backend == "lanczos"]))
  expect_error(eigs(undeclared, k = 2, method = "lanczos"),
               "requires the operator to be hermitian; it declares hermitian = NA")
  expect_match(fit$certificate$dispatch$reason, "not established hermitian")
})

test_that("on a hermitian operator Arnoldi agrees with Lanczos", {
  A <- laplacian_1d(80)
  la <- eigs(A, k = 4, which = "largest_algebraic", method = "lanczos")
  ar <- eigs(A, k = 4, which = "largest_algebraic", method = "arnoldi")
  expect_equal(ar$method, "arnoldi")
  expect_true(is.double(ar$values))
  expect_equal(ar$values, la$values, tolerance = 1e-12)
})

test_that("shift-invert with a real shift and with a complex one", {
  n <- 80; mu <- 0.04
  A <- convdiff_1d(n, mu)
  truth <- convdiff_1d_eigenvalues(n, mu)
  fit <- eigs(A, k = 3, sigma = 1.3)
  expect_equal(fit$method, "arnoldi shift-invert")
  expect_equal(sort(fit$values), sort(truth[order(abs(truth - 1.3))[1:3]]), tolerance = 1e-10)

  lam <- complex(real = seq(1, 3, length.out = 20), imaginary = seq(0.5, 2, length.out = 20))
  M <- real_with_pairs(lam, seed = 3L)
  truth <- c(lam, Conj(lam))
  sigma <- 2 + 1i
  fz <- eigs(linop(M), k = 2, sigma = sigma)
  want <- truth[order(Mod(truth - sigma))[1:2]]
  expect_lt(nearest_error(fz$values, want), 1e-10)
  ## a complex shift on a hermitian operator passes over Lanczos to Arnoldi
  h <- eigs(linop(as.matrix(laplacian_1d(40))), k = 2, sigma = 1 + 0.1i)
  expect_equal(h$method, "arnoldi shift-invert")
  expect_match(h$certificate$dispatch$candidates$reason[
    h$certificate$dispatch$candidates$backend == "lanczos"], "complex shift")
})

test_that("an exhausted budget comes back as a fail certificate", {
  A <- convdiff_1d(80, 0.04)
  fit <- eigs(A, k = 4, which = "largest_real", maxit = 30L)
  expect_equal(fit$certificate$overall, "fail")
  expect_lt(fit$nconv, 4L)
  expect_match(fit$certificate$checks$detail[fit$certificate$checks$check == "convergence"],
               "budget exhausted")
})

test_that("an invariant subspace reached from the start vector is a restart, not an end", {
  ## block diagonal of non-normal 2 x 2 blocks; a start vector inside the first
  ## block spans an invariant subspace after two steps
  blocks <- list(matrix(c(4, 0, 1, 3), 2), matrix(c(2, 0, 5, 1.5), 2),
                 matrix(c(-1, 0, 2, -0.5), 2))
  M <- matrix(0, 6, 6)
  for (i in 1:3) M[(2 * i - 1):(2 * i), (2 * i - 1):(2 * i)] <- blocks[[i]]
  fit <- eigs(linop(M), k = 3, which = "largest_real", v0 = c(1, 1, 0, 0, 0, 0), ncv = 5L)
  expect_equal(fit$nconv, 3L)
  ## the eigenvalue 2 couples by 5 across a gap of 0.5 to its neighbour 1.5, so
  ## a backward error at tol moves it by about ten times that
  expect_equal(fit$values, c(4, 3, 2), tolerance = 1e-8)
})

## ---------------------------------------------------------- certificate --

test_that("the certificate makes no claim a non-normal operator does not support", {
  A <- convdiff_1d(80, 0.04)
  cert <- eigs(A, k = 3, which = "largest_real")$certificate
  st <- function(row) cert$checks$status[cert$checks$check == row]
  expect_equal(st("backward error"), "pass")
  expect_equal(st("residual"), "pass")
  expect_equal(st("orthogonality"), "not_checked")
  expect_equal(st("forward error"), "not_checked")
  expect_match(cert$checks$detail[cert$checks$check == "forward error"], "pseudospectrum")

  ## the same method on an operator established hermitian keeps Weyl's bound
  h <- eigs(laplacian_1d(60), k = 2, which = "largest_algebraic", method = "arnoldi")$certificate
  expect_equal(h$checks$status[h$checks$check == "forward error"], "pass")
  expect_equal(h$checks$guarantee[h$checks$check == "forward error"], "deterministic_bound")
  expect_equal(h$checks$status[h$checks$check == "orthogonality"], "pass")
})

test_that("a backward error of 1e-10 says nothing about an ill-conditioned eigenvalue", {
  ## rho^n = 3^50 here: the computed values are exact for an operator within
  ## 1e-10 of A and wrong in the first digit, and the certificate claims the
  ## first and not the second.
  n <- 100; mu <- 0.5
  fit <- eigs(convdiff_1d(n, mu), k = 2, which = "largest_real")
  cert <- fit$certificate
  expect_equal(cert$checks$status[cert$checks$check == "backward error"], "pass")
  expect_equal(cert$checks$status[cert$checks$check == "forward error"], "not_checked")
  expect_gt(nearest_error(fit$values, convdiff_1d_eigenvalues(n, mu)), 1e-3)
})

test_that("over 20 seeds the k returned are the leading k and none certifies as fail", {
  for (s in 1:20) {
    n <- 60
    M <- rmat(n, n, seed = s) / sqrt(n)
    ev <- eigen(M, only.values = TRUE)$values
    for (w in c("LM", "LR")) {
      fit <- eigs(linop(M), k = 4, which = w, seed = s)
      key <- switch(w, LM = Mod, LR = Re)
      expect_false(fit$certificate$overall == "fail", info = sprintf("seed %d %s", s, w))
      expect_lt(nearest_error(fit$values, ev), 1e-9)
      expect_equal(sort(key(fit$values), decreasing = TRUE),
                   sort(key(ev), decreasing = TRUE)[1:4], tolerance = 1e-9,
                   info = sprintf("seed %d %s", s, w))
    }
  }
})

## --------------------------------------------------------- Schur kernel --

test_that("the Schur form is a unitary similarity to a triangular matrix", {
  for (cplx in c(FALSE, TRUE)) {
    for (m in c(1L, 2L, 7L, 30L)) {
      M <- rmat(m, m, seed = m)
      if (cplx) M <- M + 1i * rmat(m, m, seed = m + 100L)
      s <- linop:::schur_complex(M)
      expect_lt(max(Mod(s$U %*% s$T %*% Conj(t(s$U)) - M)), 1e-13 * max(1, max(Mod(M))) * m)
      expect_lt(max(Mod(Conj(t(s$U)) %*% s$U - diag(m))), 1e-13 * m)
      expect_true(all(s$T[lower.tri(s$T)] == 0))
      ev <- eigen(M, only.values = TRUE)$values
      expect_lt(nearest_error(diag(s$T), ev), 1e-11)
    }
  }
})

test_that("reordering brings the chosen eigenvalues to the front and keeps the similarity", {
  M <- rmat(20, 20, seed = 7L)
  s <- linop:::schur_complex(M)
  d <- diag(s$T)
  front <- order(Mod(d))[1:4]
  r <- linop:::schur_reorder(s$T, s$U, front)
  expect_equal(diag(r$T)[1:4], d[front])
  expect_identical(r$pos[1:4], front)
  expect_lt(max(Mod(r$U %*% r$T %*% Conj(t(r$U)) - M)), 1e-12)
  expect_true(all(r$T[lower.tri(r$T)] == 0))
  ## the leading columns span the invariant subspace of those eigenvalues
  Q <- r$U[, 1:4]
  expect_lt(max(Mod(M %*% Q - Q %*% r$T[1:4, 1:4])), 1e-12)
})
