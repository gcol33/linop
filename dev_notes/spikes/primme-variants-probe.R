## Where PRIMME's time goes on the converging fixture, and what the adapter's
## choice of preset and block costs against the alternatives.
##
## The well of primme-restart-probe.R (n = 1e5 and 1e6, four lowest bound
## states, tol 1e-10), solved by the engine directly, so the times are PRIMME's
## and the operator's with no certificate. Five configurations: GD+k and JDQMR,
## each with a block of 1 and of k = 4, and PRIMME's own default, DYNAMIC,
## which switches between the two on measured timings and is run twice to show
## that it does not repeat itself. PRIMME's counters split the elapsed time into
## the applies, the orthogonalisation and the dense Rayleigh-Ritz updates.
##
## From the repository root: Rscript dev_notes/spikes/primme-variants-probe.R

if (!file.exists("DESCRIPTION")) stop("run this from the repository root")
OUT_DIR <- file.path("dev_notes", "spikes", "results")
CSV <- file.path(OUT_DIR, "primme-variants-probe.csv")

well <- function(n) {
  v <- numeric(n)
  mid <- n %/% 2
  v[(mid - 14):(mid + 15)] <- 1
  function(X) {
    X <- as.matrix(X)
    Y <- 2 * X
    Y[-n, ] <- Y[-n, , drop = FALSE] - X[-1L, , drop = FALSE]
    Y[-1L, ] <- Y[-1L, , drop = FALSE] - X[-n, , drop = FALSE]
    Y - v * X
  }
}

configs <- list(
  list(label = "GD+k, block 1", preset = "PRIMME_DEFAULT_MIN_MATVECS", block = 1L),
  list(label = "GD+k, block k", preset = "PRIMME_DEFAULT_MIN_MATVECS", block = 4L),
  list(label = "JDQMR, block 1", preset = "PRIMME_DEFAULT_MIN_TIME", block = 1L),
  list(label = "JDQMR, block k", preset = "PRIMME_DEFAULT_MIN_TIME", block = 4L),
  list(label = "DYNAMIC, block 1, run 1", preset = "PRIMME_DYNAMIC", block = 1L),
  list(label = "DYNAMIC, block 1, run 2", preset = "PRIMME_DYNAMIC", block = 1L))

rows <- list()
for (n in c(1e5, 1e6)) {
  ap <- well(n)
  set.seed(1)
  v0 <- rnorm(n)
  for (cf in configs) {
    invisible(gc())
    t <- system.time(f <- linop.primme::primme_eigsh(
      ap, n, 4L, which = "SA", tol = 1e-10, maxit = 20000L, block = cf$block,
      v0 = v0, preset = cf$preset))[["elapsed"]]
    s <- f$stats
    rows[[length(rows) + 1L]] <- data.frame(
      n = n, config = cf$label, seconds = round(t, 2), applies = s$numMatvecs,
      nconv = f$nconv, basis = f$params$maxBasisSize,
      matvec_share = round(s$timeMatvec / s$elapsedTime, 3),
      ortho_share = round(s$timeOrtho / s$elapsedTime, 3),
      dense_share = round(s$timeDense / s$elapsedTime, 3),
      ortho_inner_per_apply = round(s$numOrthoInnerProds / s$numMatvecs, 1),
      dense_flops_per_apply_per_n = round(s$flopsDense / s$numMatvecs / n, 1),
      value1 = sprintf("%.12f", sort(f$values)[1]))
    cat(sprintf("n = %.0e  %-24s %6.2f s  %4d applies\n", n, cf$label, t,
                as.integer(s$numMatvecs)))
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, CSV, row.names = FALSE)
print(res, row.names = FALSE)
