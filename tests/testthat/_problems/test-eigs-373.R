# Extracted from test-eigs.R:373

# prequel ----------------------------------------------------------------------
pair_backward <- function(A, values, vectors, anorm) {
  R <- linop_apply(A, vectors, "N") - vectors * rep(values, each = nrow(vectors))
  max(sqrt(colSums(Mod(R)^2)) / sqrt(colSums(Mod(vectors)^2))) / anorm
}
nearest_gap <- function(values, truth) {
  vapply(values, function(v) min(abs(v - truth)), numeric(1))
}

# test -------------------------------------------------------------------------
n <- 500
A <- laplacian_1d(n)
for (k in c(1L, 3L)) {
    fit <- eigs(A, k = k, which = "smallest_algebraic")
    expect_equal(fit$nconv, k)
    expect_equal(fit$values, laplacian_1d_eigenvalues(n)[seq_len(k)], tolerance = 1e-12)
    expect_equal(cert_status(fit$certificate, "convergence"), "pass")
  }
