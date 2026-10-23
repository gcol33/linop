## Applies eigs() and svds() need to converge at the default subspace, for the
## default budget EIGEN_MAXIT_PER_DIM in R/eigen-common.R.
##
## Run from the repository root:
##   Rscript dev_notes/spikes/eigs-budget-sweep.R out.csv
##
## Every fixture runs at the default ncv with a budget of 60000, far above what
## any converging run here needs, so `applies` is what convergence costs and
## `nconv < k` with applies below the cap is a run that stopped improving. The
## reference engine and PRIMME (when linop.primme is installed) both run.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-linop.R")
out_file <- commandArgs(TRUE)[1]
if (is.na(out_file)) out_file <- "eigs-budget-sweep.csv"
CAP <- 60000L

sym_from_spectrum <- function(lam, seed) {
  set.seed(seed)
  Q <- qr.Q(qr(matrix(rnorm(length(lam)^2), length(lam))))
  M <- Q %*% (lam * t(Q))
  (M + t(M)) / 2
}

problems <- list()
for (n in c(50, 100, 200, 500, 1000, 2000)) {
  problems[[sprintf("lap1d_%d", n)]] <- list(A = laplacian_1d(n), n = n, kind = "eigs")
}
for (n in c(200, 500)) {
  set.seed(n); G <- matrix(rnorm(n * n), n); S <- (G + t(G)) / sqrt(2 * n)
  problems[[sprintf("goe_%d", n)]] <- list(A = linop(S), n = n, kind = "eigs")
  problems[[sprintf("spd_spread_%d", n)]] <-
    list(A = linop(sym_from_spectrum(exp(seq(0, log(1e4), length.out = n)), n)),
         n = n, kind = "eigs")
  problems[[sprintf("spd_cluster_%d", n)]] <-
    list(A = linop(sym_from_spectrum(c(1 + (1:(n - 20)) / n, 10 + 1e-3 * (1:20)), n + 1)),
         n = n, kind = "eigs")
}
for (n in c(100, 300, 1000)) {
  problems[[sprintf("diff1d_%d", n)]] <- list(A = diff_1d(n), n = n, kind = "svds")
}
for (n in c(200, 500)) {
  problems[[sprintf("randrect_%d", n)]] <- list(A = linop(rmat(n + 50, n, seed = n)),
                                                n = n, kind = "svds")
}

has_primme <- requireNamespace("linop.primme", quietly = TRUE)
rows <- list()
for (nm in names(problems)) {
  p <- problems[[nm]]
  whiches <- if (p$kind == "eigs") c("largest_algebraic", "smallest_algebraic")
             else c("largest", "smallest")
  methods <- if (p$kind == "eigs") c("lanczos", if (has_primme) "primme")
             else c("golub-kahan", if (has_primme) "primme")
  for (k in c(1L, 3L, 5L, 10L)) for (w in whiches) for (m in methods) {
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      if (p$kind == "eigs") eigs(p$A, k = k, which = w, method = m, maxit = CAP)
      else svds(p$A, k = k, which = w, method = m, maxit = CAP),
      error = function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    r <- if (inherits(fit, "error")) {
      data.frame(problem = nm, n = p$n, k = k, which = w, method = m,
                 ncv = NA, applies = NA, nconv = NA, restarts = NA,
                 stop = conditionMessage(fit), secs = el)
    } else {
      data.frame(problem = nm, n = p$n, k = k, which = w, method = m,
                 ncv = if (m == "primme") NA else default_ncv(p$A, k),
                 applies = fit$iterations, nconv = fit$nconv, restarts = fit$restarts,
                 stop = "", secs = el)
    }
    print(r)
    rows[[length(rows) + 1L]] <- r
    write.csv(do.call(rbind, rows), out_file, row.names = FALSE)
  }
}
