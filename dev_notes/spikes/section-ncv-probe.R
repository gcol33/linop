## Why a finite section gets a wider default subspace than every other operator.
##
## Run from the repository root:  Rscript dev_notes/spikes/section-ncv-probe.R
##
## A section's spectrum is a discretisation of the essential spectrum together
## with whatever separated from it. Enlarging n densifies that cluster instead of
## adding isolated eigenvalues, so the question is whether a fixed subspace stays
## adequate as the section grows, and whether it has to grow with n if it does not.

suppressMessages(devtools::load_all(".", quiet = TRUE))

top <- function(H, n, ...) {
  eigs(finite_section(H, n = n), k = 1, which = "largest_algebraic", ...)
}

row <- function(v, n, ncv, maxit = 2000L) {
  truth <- sqrt(v^2 + 4)
  fit <- top(linop_jacobi(diagonal = v), n, ncv = ncv, maxit = maxit)
  cv <- fit$certificate$values
  data.frame(v = v, n = n, ncv = ncv, err = abs(fit$values - truth),
             residual = cv$residual, truncation = cv$truncation,
             steps = fit$iterations, nconv = fit$nconv)
}

cat("== the default subspace against the section size ==\n")
cat("lambda = sqrt(v^2 + 4); a small v puts it just outside the band [-2, 2].\n\n")
out <- do.call(rbind, c(
  lapply(c(20, 40, 80, 140), function(n) row(0.3, n, 21L, maxit = 300L)),
  lapply(c(20, 40, 80, 140), function(n) row(0.3, n, 40L, maxit = 300L))))
print(out, row.names = FALSE, digits = 4)

cat("\n== does a constant subspace survive a growing section ==\n\n")
big <- do.call(rbind, unlist(lapply(c(0.3, 0.1), function(v) {
  lapply(c(80, 200, 400, 800), function(n) row(v, n, 40L))
}), recursive = FALSE))
print(big, row.names = FALSE, digits = 4)

cat("\n== several bound states, where k rather than the cluster sets the width ==\n")
cat("a well of depth 2 and width 5 carries more than one eigenvalue above the band.\n\n")
W <- linop_jacobi(diagonal = rep(2, 5))
ref <- sort(eigen(as.matrix(finite_section(W, n = 400)),
                  symmetric = TRUE, only.values = TRUE)$values,
            decreasing = TRUE)[1:3]
multi <- do.call(rbind, lapply(c(21L, 40L, 80L), function(ncv) {
  fit <- eigs(finite_section(W, n = 120), k = 3, which = "largest_algebraic",
              ncv = ncv, maxit = 2000L)
  data.frame(ncv = ncv, nconv = fit$nconv, steps = fit$iterations,
             worst_err = max(abs(sort(fit$values, decreasing = TRUE) - ref)),
             worst_res = max(fit$certificate$values$residual))
}))
print(multi, row.names = FALSE, digits = 4)
