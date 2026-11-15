# Extracted from test-eigs-arnoldi.R:140

# prequel ----------------------------------------------------------------------
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

# test -------------------------------------------------------------------------
n <- 80
mu <- 0.04
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
h <- eigs(laplacian_1d(40), k = 2, sigma = 1 + 0.1i)
expect_equal(h$method, "arnoldi shift-invert")
expect_match(h$certificate$dispatch$candidates$reason[
    h$certificate$dispatch$candidates$backend == "lanczos"], "complex shift")
