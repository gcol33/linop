# Extracted from test-backends.R:114

# prequel ----------------------------------------------------------------------
with_priority <- function(priority, code) {
  old <- options(linop.backend_priority = priority)
  on.exit(options(old), add = TRUE)
  force(code)
}
with_default_priority <- function(code) with_priority(NULL, code)
with_backend <- function(name, family, code, ...) {
  key <- linop:::register_backend(name, family, ...)
  on.exit(rm(list = key, envir = linop:::.backends), add = TRUE)
  force(code)
}
triplet_backward <- function(A, d, u, v, anorm) {
  ru <- linop:::linop_apply(A, v, "N") - u * rep(d, each = nrow(u))
  rv <- linop:::linop_apply(A, u, "C") - v * rep(d, each = nrow(v))
  max(sqrt(colSums(Mod(ru)^2) + colSums(Mod(rv)^2)) / sqrt(2)) / anorm
}
rq_on <- function(A, X) {
  AX <- linop:::linop_apply(A, X, "N")
  Re(colSums(Conj(X) * AX)) / Re(colSums(Conj(X) * X))
}

# test -------------------------------------------------------------------------
skip_if_not_installed("linop.primme")
A <- laplacian_1d(30)
for (tol in c(0, 1e-17)) {
    fit <- with_default_priority(eigs(A, k = 2, tol = tol, maxit = 40L))
    expect_equal(fit$method, "lanczos")
    cand <- fit$certificate$dispatch$candidates
    expect_match(cand$reason[cand$backend == "primme"], "below machine epsilon")
    expect_error(eigs(A, k = 2, tol = tol, method = "primme"), "below machine epsilon")
  }
