# Extracted from test-eigs-generalized.R:42

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
n <- 100
K <- laplacian_1d(n)
M <- linop(fem_mass(n), properties = c(positive_definite = TRUE))
truth <- fem_eigenvalues(n)
top <- eigs(K, k = 3, which = "largest_algebraic", B = M)
expect_equal(top$method, "lanczos")
expect_equal(top$values, rev(truth)[1:3], tolerance = 1e-12)
expect_lt(pencil_residual(K, M, top$vectors, top$values), 1e-9)
G <- crossprod(top$vectors, fem_mass(n) %*% top$vectors)
expect_lt(max(abs(G - diag(3))), 1e-12)
