## LSMR against the published recurrence, and against LSQR on the quantity the
## certificate reports.
##
##   Rscript dev_notes/spikes/lsmr_reference.R

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

lsmr <- linop:::lsmr_solve
lsqr <- linop:::lsqr_solve

source("dev_notes/spikes/lsmr_reference_fn.R")

cat("=== 1. step-for-step agreement with the published recurrence ===\n")
for (steps in c(1, 2, 4, 6, 8, 10, 12, 15, 20, 25)) {
  worst <- 0
  for (seed in 1:5) {
    f <- lsq_prescribed(45, 15, exp(seq(log(0.5), log(20), length.out = 15)), seed = seed)
    set.seed(600L + seed)
    b <- matrix(stats::rnorm(45), 45, 1)
    mine <- lsmr(linop(f$A), b, tol = 0, maxit = steps)
    ref <- reference_lsmr(f$A, b, steps)$x
    worst <- max(worst, max(Mod(mine$x - ref)) / max(Mod(ref)))
  }
  cat(sprintf("  steps %2d   worst relative gap %.3e\n", steps, worst))
}

cat("\n=== 2. the recurrence estimates against the true quantities ===\n")
f <- lsq_prescribed(45, 15, exp(seq(log(0.5), log(20), length.out = 15)), seed = 1)
M <- f$A
set.seed(601)
b <- matrix(stats::rnorm(45), 45, 1)
ref <- reference_lsmr(M, b, 30)
cat("  step   est ||r||    true ||r||   est ||A^H r||  true ||A^H r||\n")
for (k in c(1, 2, 4, 8, 12, 15, 18, 20, 25, 30)) {
  xk <- reference_lsmr(M, b, k)$x
  r <- b - M %*% xk
  cat(sprintf("  %4d   %.5e  %.5e  %.5e   %.5e\n", k,
              ref$normr[k], sqrt(sum(r^2)),
              ref$normar[k], sqrt(sum((crossprod(M, r))^2))))
}

cat("\n=== 3. monotonicity of ||A^H r||, LSMR against LSQR ===\n")
cat("  the certificate's least-squares backward error has ||A^H r|| in the numerator\n")
for (seed in 1:6) {
  f <- lsq_prescribed(50, 16, exp(seq(log(0.2), log(30), length.out = 16)), seed = seed)
  M <- f$A
  set.seed(700L + seed)
  b <- matrix(stats::rnorm(50), 50, 1)
  A <- linop(M)

  at <- function(solver, k) {
    x <- solver(A, b, tol = 0, maxit = k)$x
    sqrt(sum(Mod(crossprod(Conj(M), b - M %*% x))^2))
  }
  ks <- seq_len(24)
  a_lsmr <- vapply(ks, function(k) at(lsmr, k), numeric(1))
  a_lsqr <- vapply(ks, function(k) at(lsqr, k), numeric(1))
  rise <- function(v) sum(diff(v) > 0)
  worst_rise <- function(v) { d <- diff(v) / utils::head(v, -1); max(c(0, d)) }
  cat(sprintf("  seed %d  lsmr: %2d rises, worst +%6.2f%%   lsqr: %2d rises, worst +%8.2f%%\n",
              seed, rise(a_lsmr), 100 * worst_rise(a_lsmr),
              rise(a_lsqr), 100 * worst_rise(a_lsqr)))
}

cat("\n=== 4. recovery and cost, against the closed form ===\n")
for (seed in 1:5) {
  f <- lsq_prescribed(50, 16, exp(seq(log(0.2), log(30), length.out = 16)), seed = seed)
  M <- f$A
  set.seed(800L + seed)
  b <- matrix(stats::rnorm(50), 50, 1)
  A <- linop(M)
  truth <- lsq_prescribed_solve(f, b)
  a <- lsmr(A, b, tol = 1e-11, maxit = 2000L)
  c_ <- lsqr(A, b, tol = 1e-11, maxit = 2000L)
  cat(sprintf("  seed %d  lsmr %3d it, rel err %.2e, conv %-5s | lsqr %3d it, rel err %.2e, conv %s\n",
              seed, a$iterations, max(Mod(a$x - truth)) / max(Mod(truth)), a$converged,
              c_$iterations, max(Mod(c_$x - truth)) / max(Mod(truth)), c_$converged))
}

cat("\n=== 5. finite termination in floating point ===\n")
f <- lsq_prescribed(45, 15, exp(seq(log(0.5), log(20), length.out = 15)), seed = 2)
M <- f$A
set.seed(900)
b <- matrix(stats::rnorm(45), 45, 1)
truth <- lsq_prescribed_solve(f, b)
for (k in c(10, 15, 18, 20, 22, 25, 30)) {
  x <- reference_lsmr(M, b, k)$x
  cat(sprintf("  step %2d   relative error %.3e\n", k,
              max(Mod(x - truth)) / max(Mod(truth))))
}

cat("\n=== 6. complex, rank deficient, and lockstep ===\n")
fz <- lsq_prescribed(40, 15, seq(1, 7, length.out = 15), seed = 4, dtype = "complex")
set.seed(13)
bz <- matrix(complex(real = stats::rnorm(40), imaginary = stats::rnorm(40)), 40, 1)
fitz <- lsmr(linop(fz$A), bz, tol = 1e-12)
cat(sprintf("  complex: conv %s, rel err %.3e\n", fitz$converged,
            max(Mod(fitz$x - lsq_prescribed_solve(fz, bz))) /
              max(Mod(lsq_prescribed_solve(fz, bz)))))

n <- 40
D <- diff_1d(n)
set.seed(14)
bd <- matrix(stats::rnorm(n - 1L), n - 1L, 1)
fitd <- lsmr(D, bd, tol = 1e-12)
cat(sprintf("  rank deficient: conv %s, rel err %.3e, sum(x) %.3e\n", fitd$converged,
            max(Mod(fitd$x - diff_1d_min_norm_solve(n, bd))) /
              max(Mod(diff_1d_min_norm_solve(n, bd))), sum(fitd$x)))

f <- lsq_prescribed(50, 16, exp(seq(log(0.3), log(30), length.out = 16)), seed = 9)
A <- linop(f$A)
set.seed(21)
BB <- matrix(stats::rnorm(50 * 4), 50, 4)
blk <- lsmr(A, BB, tol = 1e-11)
ok <- TRUE
for (j in 1:4) {
  one <- lsmr(A, BB[, j, drop = FALSE], tol = 1e-11)
  ok <- ok && identical(blk$x[, j], one$x[, 1L])
}
cat(sprintf("  lockstep bitwise: %s\n", ok))

cat("\n=== 7. the floor on the least-squares line ===\n")
f <- lsq_prescribed(40, 12, seq(1, 4, length.out = 12), seed = 28)
A <- linop(f$A)
set.seed(29)
b <- matrix(stats::rnorm(40), 40, 1)
achieved <- lsmr(A, b, tol = 0, maxit = 200L)$certificate$values$backward_error
cat(sprintf("  achieved backward error %.3e, c eps = %.3e\n",
            achieved, linop:::SOLVE_FLOOR_CONST * .Machine$double.eps))
wf <- lsmr(A, b, tol = achieved / 2, maxit = 200L)
wo <- lsmr(A, b, tol = achieved / 2, maxit = 200L, floor_const = 0)
cat(sprintf("  identical iterates %s | with floor %s / %s | without %s / %s\n",
            identical(wf$x, wo$x),
            linop:::cert_status(wf$certificate, "backward error"),
            linop:::cert_status(wf$certificate, "convergence"),
            linop:::cert_status(wo$certificate, "backward error"),
            linop:::cert_status(wo$certificate, "convergence")))

cat("\n=== 8. ill-conditioned, against a dense reference ===\n")
for (kappa in c(1e4, 1e8, 1e10)) {
  nn <- 20
  sg <- exp(seq(log(1), log(1 / kappa), length.out = nn))
  f <- lsq_prescribed(60, nn, sg, seed = 5)
  M <- f$A
  set.seed(1000)
  b <- matrix(stats::rnorm(60), 60, 1)
  truth <- lsq_prescribed_solve(f, b)
  a <- lsmr(linop(M), b, tol = 1e-12, maxit = 5000L)
  c_ <- lsqr(linop(M), b, tol = 1e-12, maxit = 5000L)
  cat(sprintf("  kappa %.0e  lsmr %4d it rel err %.2e | lsqr %4d it rel err %.2e\n",
              kappa, a$iterations, max(Mod(a$x - truth)) / max(Mod(truth)),
              c_$iterations, max(Mod(c_$x - truth)) / max(Mod(truth))))
}
