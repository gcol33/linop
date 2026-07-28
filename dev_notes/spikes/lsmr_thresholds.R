## The numbers the LSMR test file's thresholds are set from.
##
##   Rscript dev_notes/spikes/lsmr_thresholds.R

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")
source("dev_notes/spikes/lsmr_reference_fn.R")

lsmr <- linop:::lsmr_solve
lsqr <- linop:::lsqr_solve

cat("=== reference agreement at 4 steps, across conditioning ===\n")
for (kappa in c(1e2, 1e4, 1e6, 1e8, 1e10, 1e12)) {
  worst <- 0
  for (seed in 1:5) {
    nn <- 20
    sg <- exp(seq(log(1), log(1 / kappa), length.out = nn))
    f <- lsq_prescribed(60, nn, sg, seed = seed)
    set.seed(600L + seed)
    b <- matrix(stats::rnorm(60), 60, 1)
    mine <- lsmr(linop(f$A), b, tol = 0, maxit = 4L)$x
    ref <- reference_lsmr(f$A, b, 4L)$x
    worst <- max(worst, max(Mod(mine - ref)) / max(Mod(ref)))
  }
  cat(sprintf("  kappa %.0e   worst relative gap over 5 seeds %.3e\n", kappa, worst))
}

cat("\n=== exhaustion: a Krylov space that closes early ===\n")
## b = A v for a single right singular vector, so alpha_2 = 0 at step 1.
f <- lsq_prescribed(12, 4, c(1, 2, 3, 4), seed = 1)
v1 <- f$V[, 1L, drop = FALSE]
b <- f$A %*% v1
fit <- lsmr(linop(f$A), b, tol = 1e-12)
cat(sprintf("  iterations %d, rel err %.3e, converged %s\n", fit$iterations,
            max(Mod(fit$x - v1)) / max(Mod(v1)), fit$converged))
fitq <- lsqr(linop(f$A), b, tol = 1e-12)
cat(sprintf("  lsqr for comparison: iterations %d, rel err %.3e\n", fitq$iterations,
            max(Mod(fitq$x - v1)) / max(Mod(v1))))

cat("\n=== monotonicity of the certificate's numerator, more seeds ===\n")
for (kappa in c(1e2, 1e4, 1e6)) {
  lsmr_rises <- integer(0); lsqr_rises <- integer(0)
  worst_lsmr <- 0
  for (seed in 1:12) {
    nn <- 16
    sg <- exp(seq(log(1), log(1 / kappa), length.out = nn))
    f <- lsq_prescribed(50, nn, sg, seed = seed)
    M <- f$A
    set.seed(700L + seed)
    b <- matrix(stats::rnorm(50), 50, 1)
    A <- linop(M)
    at <- function(solver, k) {
      x <- solver(A, b, tol = 0, maxit = k)$x
      sqrt(sum(Mod(crossprod(Conj(M), b - M %*% x))^2))
    }
    ks <- seq_len(20)
    a <- vapply(ks, function(k) at(lsmr, k), numeric(1))
    q <- vapply(ks, function(k) at(lsqr, k), numeric(1))
    rel <- diff(a) / utils::head(a, -1)
    lsmr_rises <- c(lsmr_rises, sum(rel > 0))
    worst_lsmr <- max(worst_lsmr, max(c(0, rel)))
    lsqr_rises <- c(lsqr_rises, sum(diff(q) > 0))
  }
  cat(sprintf("  kappa %.0e  lsmr rises per seed: max %d, worst relative +%.3e | lsqr rises: min %d, max %d\n",
              kappa, max(lsmr_rises), worst_lsmr, min(lsqr_rises), max(lsqr_rises)))
}

cat("\n=== shared budget on the certificate's own quantity, 12 seeds ===\n")
for (kappa in c(1e4, 1e6, 1e8)) {
  for (budget in c(20L, 40L, 80L)) {
    wins <- 0; ratios <- numeric(0)
    for (seed in 1:12) {
      nn <- 20
      sg <- exp(seq(log(1), log(1 / kappa), length.out = nn))
      f <- lsq_prescribed(60, nn, sg, seed = seed)
      set.seed(1200L + seed)
      b <- matrix(stats::rnorm(60), 60, 1)
      A <- linop(f$A)
      a <- lsmr(A, b, tol = 0, maxit = budget)$certificate$values$backward_error
      q <- lsqr(A, b, tol = 0, maxit = budget)$certificate$values$backward_error
      if (a < q) wins <- wins + 1
      ratios <- c(ratios, q / a)
    }
    cat(sprintf("  kappa %.0e budget %3d   lsmr smaller on %2d of 12, ratio lsqr/lsmr min %.2f median %.2f\n",
                kappa, budget, wins, min(ratios), stats::median(ratios)))
  }
}

cat("\n=== the floor window ===\n")
for (seed in 1:5) {
  f <- lsq_prescribed(40, 12, seq(1, 4, length.out = 12), seed = 27 + seed)
  A <- linop(f$A)
  set.seed(29L + seed)
  b <- matrix(stats::rnorm(40), 40, 1)
  achieved <- lsmr(A, b, tol = 0, maxit = 200L)$certificate$values$backward_error
  cat(sprintf("  seed %d achieved %.3e against c eps %.3e\n", seed, achieved,
              linop:::SOLVE_FLOOR_CONST * .Machine$double.eps))
}
