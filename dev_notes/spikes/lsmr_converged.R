## This implementation run to convergence on the same exported fixtures, so the
## forward error at high condition number can be attributed to the method rather
## than to the code.
##
##   Rscript dev_notes/spikes/lsmr_converged.R truth    # before the python side
##   Rscript dev_notes/spikes/lsmr_converged.R compare  # after it

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

dir <- Sys.getenv("LSMR_XDIR")
cases <- list(
  list(name = "kappa1e2",  kappa = 1e2,  m = 60, n = 20, seed = 1),
  list(name = "kappa1e6",  kappa = 1e6,  m = 60, n = 20, seed = 2),
  list(name = "kappa1e10", kappa = 1e10, m = 60, n = 20, seed = 3),
  list(name = "tall",      kappa = 1e4,  m = 200, n = 30, seed = 4))

build <- function(cs) {
  sg <- exp(seq(log(1), log(1 / cs$kappa), length.out = cs$n))
  f <- lsq_prescribed(cs$m, cs$n, sg, seed = cs$seed)
  set.seed(5000L + cs$seed)
  b <- matrix(stats::rnorm(cs$m), cs$m, 1)
  list(f = f, M = f$A, b = b, truth = lsq_prescribed_solve(f, b))
}

action <- commandArgs(trailingOnly = TRUE)[1L]

if (identical(action, "truth")) {
  for (cs in cases) {
    d <- build(cs)
    utils::write.table(d$truth, file.path(dir, paste0(cs$name, "-truth.txt")),
                       row.names = FALSE, col.names = FALSE)
  }
  cat("truth written\n")
}

if (identical(action, "compare")) {
  for (cs in cases) {
    d <- build(cs)
    A <- linop(d$M)
    sv1 <- svd(d$M, nu = 0L, nv = 0L)$d[1L]
    line <- function(label, fit) {
      r <- d$b - d$M %*% fit$x
      bw <- sqrt(sum(Mod(crossprod(Conj(d$M), r))^2)) / (sv1 * sqrt(sum(Mod(r)^2)))
      sprintf("%s %4d it  fwd %.2e  bw %.2e", label, fit$iterations,
              max(Mod(fit$x - d$truth)) / max(Mod(d$truth)), bw)
    }
    cat(sprintf("  %-10s  %s  |  %s\n", cs$name,
                line("lsmr", linop:::lsmr_solve(A, d$b, tol = 1e-12, maxit = 5000L)),
                line("lsqr", linop:::lsqr_solve(A, d$b, tol = 1e-12, maxit = 5000L))))

    ## and against the converged SciPy iterate directly
    for (m in c("lsmr", "lsqr")) {
      fs <- file.path(dir, sprintf("%s-scipy-conv-%s.txt", cs$name, m))
      if (!file.exists(fs)) next
      theirs <- as.matrix(utils::read.table(fs))
      mine <- if (m == "lsmr") linop:::lsmr_solve(A, d$b, tol = 1e-12, maxit = 5000L)$x
              else linop:::lsqr_solve(A, d$b, tol = 1e-12, maxit = 5000L)$x
      cat(sprintf("               %s against scipy's own converged iterate: %.3e\n",
                  m, max(Mod(mine - theirs)) / max(Mod(theirs))))
    }
  }
}
