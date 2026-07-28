## BiCGSTAB against the published recurrence, and the numbers the test file's
## thresholds are set from.
##
##   Rscript dev_notes/spikes/bicgstab_reference.R

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

bicgstab <- linop:::bicgstab_solve

## van der Vorst 1992, single column, no preconditioner and no stopping test, so
## agreement is a statement about the recurrence rather than about the
## bookkeeping around it.
reference_bicgstab <- function(M, b, steps) {
  n <- ncol(M)
  x <- rep(0, n)
  if (is.complex(M) || is.complex(b)) storage.mode(x) <- "complex"
  r <- as.vector(b)
  rhat <- r
  p <- r
  rho <- sum(Conj(rhat) * r)
  for (j in seq_len(steps)) {
    v <- as.vector(M %*% p)
    alpha <- rho / sum(Conj(rhat) * v)
    s <- r - alpha * v
    tt <- as.vector(M %*% s)
    omega <- sum(Conj(tt) * s) / sum(Conj(tt) * tt)
    x <- x + alpha * p + omega * s
    r <- s - omega * tt
    rho_new <- sum(Conj(rhat) * r)
    beta <- (rho_new / rho) * (alpha / omega)
    p <- r + beta * (p - omega * v)
    rho <- rho_new
  }
  matrix(x, n, 1L)
}

fixtures <- list(
  "convdiff mu=.3" = as.matrix(convdiff_1d(40, 0.3)),
  "convdiff mu=.7" = as.matrix(convdiff_1d(40, 0.7)),
  "laplacian" = as.matrix(laplacian_1d(40)),
  "kms rho=.7" = kms_matrix(40, 0.7))
set.seed(9)
fixtures[["complex dense"]] <- zmat(30, 30) + diag(8, 30)

cat("=== step-for-step agreement with the published recurrence ===\n")
cat(sprintf("  %-16s %s\n", "fixture",
            paste(sprintf("%9s", paste0("k=", c(1, 2, 4, 8, 12, 16, 24))), collapse = "")))
for (nm in names(fixtures)) {
  M <- fixtures[[nm]]
  n <- nrow(M)
  out <- character(0)
  for (k in c(1, 2, 4, 8, 12, 16, 24)) {
    worst <- 0
    for (seed in 1:5) {
      set.seed(300L + seed)
      b <- if (is.complex(M))
        matrix(complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)), n, 1)
      else matrix(stats::rnorm(n), n, 1)
      mine <- bicgstab(linop(M), b, tol = 0, maxit = k)
      ref <- reference_bicgstab(M, b, k)
      worst <- max(worst, max(Mod(mine$x - ref)) / max(Mod(ref)))
    }
    out <- c(out, sprintf("%9.1e", worst))
  }
  cat(sprintf("  %-16s %s\n", nm, paste(out, collapse = "")))
}

cat("\n=== the breakdown guard is silent where it is not needed ===\n")
for (nm in names(fixtures)) {
  M <- fixtures[[nm]]
  n <- nrow(M)
  same <- TRUE
  for (seed in 1:5) {
    set.seed(400L + seed)
    b <- if (is.complex(M))
      matrix(complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)), n, 1)
    else matrix(stats::rnorm(n), n, 1)
    a <- bicgstab(linop(M), b, tol = 1e-11, maxit = 1000L)
    z <- bicgstab(linop(M), b, tol = 1e-11, maxit = 1000L, breakdown_tol = 0)
    same <- same && identical(a$x, z$x) && a$iterations == z$iterations
  }
  cat(sprintf("  %-16s bitwise identical with the guard off: %s\n", nm, same))
}

cat("\n=== how far the true residual rises, over seeds and fixtures ===\n")
for (nm in names(fixtures)) {
  M <- fixtures[[nm]]
  n <- nrow(M)
  A <- linop(M)
  rises <- integer(0); worst <- 0
  for (seed in 1:8) {
    set.seed(500L + seed)
    b <- if (is.complex(M))
      matrix(complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)), n, 1)
    else matrix(stats::rnorm(n), n, 1)
    v <- vapply(seq_len(20), function(k) {
      x <- bicgstab(A, b, tol = 0, maxit = k)$x
      sqrt(sum(Mod(b - M %*% x)^2))
    }, numeric(1))
    rises <- c(rises, sum(diff(v) > 0))
    worst <- max(worst, max(c(0, diff(v) / utils::head(v, -1))))
  }
  cat(sprintf("  %-16s rises per seed %d to %d, worst +%.0f%%\n", nm,
              min(rises), max(rises), 100 * worst))
}
