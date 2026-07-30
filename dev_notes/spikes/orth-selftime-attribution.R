## Whether Rprof's self time for orth_against is work, or bookkeeping.
##
## orth-allocation-share.R part 2 measured the shipped routine at 1.04x of the two
## BLAS products it cannot avoid, at the same shape where the ceiling probe
## attributes it a second of self time and calls that second compilable. Both
## cannot be true. This decides which by profiling the routine on its own, where
## the wall clock of each piece is known independently.
##
## If Rprof charges orth_against for time that direct timing says is inside
## crossprod and %*%, then linopR_frac in the ceiling probe is inflated by exactly
## that amount, and the ceiling it reports is not a ceiling.
##
## From the repository root: Rscript dev_notes/spikes/orth-selftime-attribution.R

if (!file.exists("DESCRIPTION")) stop("run this from the repository root")

devtools::load_all(".", quiet = TRUE)

OUT_DIR <- file.path("dev_notes", "spikes", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT_DIR, "orth-selftime-attribution.txt")
cat("", file = LOG)

say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOG, append = TRUE)
  flush.console()
}

bench <- function(f, seconds = 2) {
  f()
  t0 <- Sys.time(); reps <- 0L
  repeat {
    f(); reps <- reps + 1L
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) >= seconds) break
  }
  as.numeric(difftime(Sys.time(), t0, units = "secs")) / reps
}

for (cfg in list(c(1e6, 80), c(1e5, 20))) {
  n <- cfg[1]; ncv <- as.integer(cfg[2])
  say("")
  say("=== n = %.0e  ncv = %d", n, ncv)

  Q <- qr.Q(qr(matrix(stats::rnorm(n * ncv), n, ncv)))
  W <- matrix(stats::rnorm(n), n, 1L)

  ## Wall clock for each piece, timed on its own.
  t_cross <- bench(function() { crossprod(Q, W); invisible(NULL) })
  C0 <- crossprod(Q, W)
  t_prod  <- bench(function() { Q %*% C0; invisible(NULL) })
  t_norms <- bench(function() { col_norms(W); invisible(NULL) })
  t_full  <- bench(function() { orth_against(Q, W); invisible(NULL) })

  say("  wall clock, each piece timed alone")
  say("    crossprod(Q, W)   %8.2f ms  %5.1f%% of full", 1000 * t_cross, 100 * t_cross / t_full)
  say("    Q %%*%% C           %8.2f ms  %5.1f%% of full", 1000 * t_prod,  100 * t_prod  / t_full)
  say("    col_norms(W) x2    %8.2f ms  %5.1f%% of full", 2000 * t_norms, 200 * t_norms / t_full)
  say("    ----")
  say("    orth_against      %8.2f ms", 1000 * t_full)
  say("    two products      %8.2f ms  %5.1f%% of full",
      1000 * (t_cross + t_prod), 100 * (t_cross + t_prod) / t_full)

  ## Now how Rprof splits the same call.
  reps <- max(20L, ceiling(3 / t_full))
  pf <- tempfile(fileext = ".Rprof")
  Rprof(pf, interval = 0.002, gc.profiling = TRUE)
  for (i in seq_len(reps)) orth_against(Q, W)
  Rprof(NULL)
  s <- summaryRprof(pf); unlink(pf)
  self <- s$by.self; total <- s$sampling.time

  say("  Rprof over %d calls, %.3fs sampled", reps, total)
  ord <- head(self[order(-self$self.time), , drop = FALSE], 8L)
  for (i in seq_len(nrow(ord))) {
    say("    %-24s %7.3f  %5.1f%%", rownames(ord)[i], ord$self.time[i],
        100 * ord$self.time[i] / total)
  }

  oa <- self[rownames(self) == "orth_against", "self.time"]
  oa <- if (length(oa)) oa else 0
  say("  Rprof charges orth_against %.1f%% self; direct timing leaves %.1f%% outside the two products",
      100 * oa / total, 100 * (1 - (t_cross + t_prod) / t_full))

  rm(Q, W, C0); gc(FALSE)
}
