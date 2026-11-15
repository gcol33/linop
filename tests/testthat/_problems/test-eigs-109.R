# Extracted from test-eigs.R:109

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
lambda <- c(10, 10, 10, seq(2, 0.5, length.out = n - 3L))
M <- spd_prescribed(n, lambda, seed = 11L)
A <- as_spd_linop(M)
fit <- eigs(A, k = 3, which = "largest_algebraic")
expect_equal(fit$nconv, 3L)
expect_equal(fit$values, rep(10, 3), tolerance = 1e-9)
expect_lt(max(abs(crossprod(fit$vectors) - diag(3))), 1e-10)
