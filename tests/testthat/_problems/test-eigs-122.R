# Extracted from test-eigs.R:122

# prequel ----------------------------------------------------------------------
pair_backward <- function(A, values, vectors, anorm) {
  R <- linop_apply(A, vectors, "N") - vectors * rep(values, each = nrow(vectors))
  max(sqrt(colSums(Mod(R)^2)) / sqrt(colSums(Mod(vectors)^2))) / anorm
}
nearest_gap <- function(values, truth) {
  vapply(values, function(v) min(abs(v - truth)), numeric(1))
}

# test -------------------------------------------------------------------------
n <- 60
A <- laplacian_1d(n)
truth <- laplacian_1d_eigenvalues(n)
sigma <- 2
want <- truth[order(abs(truth - sigma))][1:3]
fs <- eigs(A, k = 3, sigma = sigma)
expect_equal(fs$nconv, 3L)
expect_equal(sort(fs$values), sort(want), tolerance = 1e-9)
