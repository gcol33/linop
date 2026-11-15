# Extracted from test-eigs-generalized.R:110

# prequel ----------------------------------------------------------------------
fem_mass <- function(n) {
  M <- diag(4, n)
  M[cbind(1:(n - 1), 2:n)] <- 1
  M[cbind(2:n, 1:(n - 1))] <- 1
  M / 6
}
fem_eigenvalues <- function(n) {
  t <- seq_len(n) * pi / (n + 1)
  sort(6 * (1 - cos(t)) / (2 + cos(t)))
}
pencil_residual <- function(A, B, X, values) {
  R <- linop:::linop_apply(A, X, "N") -
    linop:::linop_apply(B, X, "N") * rep(values, each = nrow(X))
  max(sqrt(colSums(Mod(R)^2)) / sqrt(colSums(Mod(X)^2)))
}
nearest_error <- function(values, truth) {
  max(vapply(values, function(v) min(Mod(truth - v)), numeric(1)))
}

# test -------------------------------------------------------------------------
n <- 60
Md <- fem_mass(n)
M <- linop(function(X) Md %*% X, dim = c(n, n))
fit <- eigs(laplacian_1d(n), k = 3, which = "largest_real", B = M)
expect_equal(fit$method, "arnoldi")
expect_equal(fit$values, rev(fem_eigenvalues(n))[1:3], tolerance = 1e-10)
cand <- fit$certificate$dispatch$candidates
expect_match(cand$reason[cand$backend == "lanczos"], "the B inner product requires B to be hermitian")
