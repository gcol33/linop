# Extracted from test-eigs.R:90

# prequel ----------------------------------------------------------------------
pair_backward <- function(A, values, vectors, anorm) {
  R <- linop_apply(A, vectors, "N") - vectors * rep(values, each = nrow(vectors))
  max(sqrt(colSums(Mod(R)^2)) / sqrt(colSums(Mod(vectors)^2))) / anorm
}
nearest_gap <- function(values, truth) {
  vapply(values, function(v) min(abs(v - truth)), numeric(1))
}

# test -------------------------------------------------------------------------
n <- 30
lambda <- c(50, 20, 8, seq(2, 0.5, length.out = n - 3L))
M <- hpd_prescribed(n, lambda, seed = 3L)
A <- linop(M, properties = c(hermitian = TRUE, positive_definite = TRUE))
fit <- eigs(A, k = 3, which = "largest_algebraic")
expect_equal(fit$nconv, 3L)
expect_equal(fit$values, lambda[1:3], tolerance = 1e-9)
