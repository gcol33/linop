## What the two methods actually certify on an ill-conditioned least-squares
## problem, and whether the difference in forward error is a difference in what
## was achieved or in what the problem allows.
##
##   Rscript dev_notes/spikes/lsmr_illcond.R

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

lsmr <- linop:::lsmr_solve
lsqr <- linop:::lsqr_solve
cert_status <- linop:::cert_status

report <- function(label, fit, M, b, truth) {
  x <- fit$x
  r <- b - M %*% x
  atr <- crossprod(Conj(M), r)
  sv <- svd(M, nu = 0L, nv = 0L)$d
  true_bw <- sqrt(sum(Mod(atr)^2)) / (sv[1L] * sqrt(sum(Mod(r)^2)))
  cat(sprintf("    %-5s %4d it  conv %-5s  cert bw %.2e  true bw %.2e  fwd err %.2e  %s\n",
              label, fit$iterations, fit$converged,
              fit$certificate$values$backward_error, true_bw,
              max(Mod(x - truth)) / max(Mod(truth)),
              cert_status(fit$certificate, "convergence")))
}

cat("=== incompatible right-hand side, growing condition number ===\n")
cat("    the forward-error bound for least squares carries kappa^2 on the\n")
cat("    incompatible part, so kappa^2 eps is the scale to read against\n\n")
for (kappa in c(1e4, 1e6, 1e8, 1e10, 1e12)) {
  nn <- 20
  sg <- exp(seq(log(1), log(1 / kappa), length.out = nn))
  cat(sprintf("  kappa %.0e   kappa^2 eps = %.2e\n", kappa, kappa^2 * .Machine$double.eps))
  for (seed in 1:3) {
    f <- lsq_prescribed(60, nn, sg, seed = seed)
    M <- f$A
    set.seed(1000L + seed)
    b <- matrix(stats::rnorm(60), 60, 1)
    truth <- lsq_prescribed_solve(f, b)
    cat(sprintf("   seed %d\n", seed))
    report("lsmr", lsmr(linop(M), b, tol = 1e-12, maxit = 5000L), M, b, truth)
    report("lsqr", lsqr(linop(M), b, tol = 1e-12, maxit = 5000L), M, b, truth)
  }
}

cat("\n=== compatible right-hand side, same operators ===\n")
cat("    with b in the range of A the kappa^2 term drops out\n\n")
for (kappa in c(1e6, 1e8, 1e10)) {
  nn <- 20
  sg <- exp(seq(log(1), log(1 / kappa), length.out = nn))
  cat(sprintf("  kappa %.0e   kappa eps = %.2e\n", kappa, kappa * .Machine$double.eps))
  for (seed in 1:3) {
    f <- lsq_prescribed(60, nn, sg, seed = seed)
    M <- f$A
    set.seed(1100L + seed)
    truth <- matrix(stats::rnorm(nn), nn, 1)
    b <- M %*% truth
    cat(sprintf("   seed %d\n", seed))
    report("lsmr", lsmr(linop(M), b, tol = 1e-12, maxit = 5000L), M, b, truth)
    report("lsqr", lsqr(linop(M), b, tol = 1e-12, maxit = 5000L), M, b, truth)
  }
}

cat("\n=== what a shared budget buys, on the certificate's own quantity ===\n")
cat("    both stopped at the same iteration count, so the comparison is of\n")
cat("    what each reached rather than of what each spent\n\n")
for (kappa in c(1e6, 1e8, 1e10)) {
  nn <- 20
  sg <- exp(seq(log(1), log(1 / kappa), length.out = nn))
  for (budget in c(30L, 60L, 120L)) {
    win_lsmr <- 0; win_lsqr <- 0
    for (seed in 1:8) {
      f <- lsq_prescribed(60, nn, sg, seed = seed)
      M <- f$A
      set.seed(1200L + seed)
      b <- matrix(stats::rnorm(60), 60, 1)
      a <- lsmr(linop(M), b, tol = 0, maxit = budget)$certificate$values$backward_error
      c_ <- lsqr(linop(M), b, tol = 0, maxit = budget)$certificate$values$backward_error
      if (a < c_) win_lsmr <- win_lsmr + 1 else win_lsqr <- win_lsqr + 1
    }
    cat(sprintf("  kappa %.0e  budget %3d   lsmr better on %d of 8, lsqr on %d\n",
                kappa, budget, win_lsmr, win_lsqr))
  }
}
