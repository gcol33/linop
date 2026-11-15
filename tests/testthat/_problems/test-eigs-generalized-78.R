# Extracted from test-eigs-generalized.R:78

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
for (s in 1:20) {
    n <- 40
    dt <- if (s %% 2) "double" else "complex"
    Ad <- hpd_prescribed(n, seq(-3, 5, length.out = n), seed = s)
    if (dt == "double") Ad <- Re(Ad)
    Bd <- hpd_prescribed(n, exp(seq(log(10), 0, length.out = n)), seed = 100L + s)
    if (dt == "double") Bd <- Re(Bd)
    Ad <- (Ad + Conj(t(Ad))) / 2; Bd <- (Bd + Conj(t(Bd))) / 2
    ref <- sort(Re(eigen(solve(Bd, Ad), only.values = TRUE)$values), decreasing = TRUE)
    A <- linop(Ad)
    B <- linop(Bd, properties = c(positive_definite = TRUE))
    fit <- eigs(A, k = 3, which = "largest_algebraic", B = B, seed = s)
    info <- sprintf("seed %d %s", s, dt)
    expect_equal(fit$method, "lanczos", info = info)
    expect_equal(fit$values, ref[1:3], tolerance = 1e-10, info = info)
    expect_false(fit$certificate$overall == "fail", info = info)
    expect_equal(cert_status(fit$certificate, "orthogonality"), "pass", info = info)
    expect_equal(cert_status(fit$certificate, "forward error"), "pass", info = info)
  }
