## What is inside orth_against's self time, and how much of it is reachable at all.
##
## The compilable-ceiling probe leaves one entry unexplained. On the largest cell
## it measures, orth_against's *self* time is the single largest entry in the run
## -- larger than `%*%` and `crossprod` -- while every operation in its four lines
## dispatches to a primitive that is listed separately. `-` is 6 ms there. So the
## second of self time is not the arithmetic and it is not the primitives.
##
## Two candidates, and they lead to opposite conclusions:
##
##   per-call overhead   interpreter cost of evaluating six calls. Fixed per
##                       iteration, so its share falls as n grows, and compiling
##                       is worth little at the sizes the package is for.
##   allocation and GC   the per-iteration temporaries, collected against a live
##                       basis of ncv vectors. Grows with n*ncv, so its share is
##                       flat in n, and compiling removes it.
##
## The crossover probe already settles which by its shape: orth_fraction is flat
## in n (0.67, 0.69, 0.67 for the stencil at ncv = 80 over n = 1e4, 1e5, 1e6).
## Fixed per-call cost cannot do that. This probe measures it directly instead of
## reading it off a trend.
##
## Part 1 re-runs that cell with gc.profiling on, so collection is attributed
## rather than folded into the frame that triggered it.
##
## Part 2 is the number the C question actually turns on: the BLAS floor at the
## same shape. orth_against cannot go below the two products it has to issue, so
## the gap between those and the whole function is the most a compiled
## implementation of this routine could remove, at this shape, with the algorithm
## unchanged.
##
## From the repository root: Rscript dev_notes/spikes/orth-allocation-share.R

if (!file.exists("DESCRIPTION")) stop("run this from the repository root")

devtools::load_all(".", quiet = TRUE)

OUT_DIR <- file.path("dev_notes", "spikes", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT_DIR, "orth-allocation-share.txt")
CSV <- file.path(OUT_DIR, "orth-allocation-share.csv")
cat("", file = LOG)

say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOG, append = TRUE)
  flush.console()
}

LINOP_NAMES <- ls(asNamespace("linop"), all.names = TRUE)
USER_APPLY <- c("user_operator_apply", "inner")

make_stencil <- function(n) {
  force(n)
  function(X) {
    X <- as.matrix(X)
    Y <- 2 * X
    Y[-n, ] <- Y[-n, , drop = FALSE] - X[-1L, , drop = FALSE]
    Y[-1L, ] <- Y[-1L, , drop = FALSE] - X[-n, , drop = FALSE]
    Y
  }
}

N <- 1e6
NCV <- 80L

## ------------------------------------------------- 1. gc attributed separately
say("=== 1. the same cell, with gc.profiling on")
say("    stencil n = %.0e ncv = %d", N, NCV)

inner <- make_stencil(N)
user_operator_apply <- function(X) inner(X)
A <- linop(user_operator_apply, dim = c(N, N), properties = c(hermitian = TRUE))
run <- function() eigs(A, k = 4L, ncv = NCV, maxit = NCV, tol = 0, seed = 1L)

run()
pf <- tempfile(fileext = ".Rprof")
Rprof(pf, interval = 0.002, gc.profiling = TRUE, memory.profiling = TRUE)
run()
Rprof(NULL)
s <- summaryRprof(pf, memory = "both")
unlink(pf)

self <- s$by.self
total <- s$sampling.time
say("    sampled %.3fs", total)
ord <- head(self[order(-self$self.time), , drop = FALSE], 12L)
for (i in seq_len(nrow(ord))) {
  nm <- rownames(ord)[i]
  bare <- gsub('^"|"$', "", nm)
  cls <- if (bare %in% USER_APPLY) "apply"
         else if (bare %in% LINOP_NAMES) "linopR" else "compiled"
  say("      %-10s %-26s %7.3f  %5.1f%%", cls, nm, ord$self.time[i],
      100 * ord$self.time[i] / total)
}

## summaryRprof reports GC time in its own column when gc.profiling is on.
if ("gc" %in% names(s)) {
  say("    gc total: %.3fs (%.1f%% of sampled)", s$gc$gc.time, 100 * s$gc$gc.time / total)
}
gcm <- self[grepl("gc", rownames(self), fixed = TRUE), , drop = FALSE]
if (nrow(gcm)) {
  for (i in seq_len(nrow(gcm))) {
    say("    gc entry %-20s %7.3f  %5.1f%%", rownames(gcm)[i],
        gcm$self.time[i], 100 * gcm$self.time[i] / total)
  }
}
rm(A, inner, user_operator_apply); gc(FALSE)

## ------------------------------------------------------------ 2. the BLAS floor
## orth_against has to issue crossprod(Q, W) and Q %*% C whatever language it is
## written in. Timing those two alone bounds what compiling the routine could buy
## with the algorithm left alone.
say("")
say("=== 2. the BLAS floor at the same shape")

bench <- function(f, seconds = 2) {
  f()
  t0 <- Sys.time(); reps <- 0L
  repeat {
    f(); reps <- reps + 1L
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) >= seconds) break
  }
  as.numeric(difftime(Sys.time(), t0, units = "secs")) / reps
}

rows <- list()
for (n in c(1e5, 1e6)) {
  for (ncv in c(20L, 80L)) {
    Q <- matrix(stats::rnorm(n * ncv), n, ncv)
    Q <- qr.Q(qr(Q))
    W <- matrix(stats::rnorm(n), n, 1L)

    ## The two products, and nothing else: no norms, no second pass, no
    ## subtraction. This is the floor any implementation shares.
    floor_t <- bench(function() {
      C <- crossprod(Q, W)
      Q %*% C
      invisible(NULL)
    })
    ## The whole routine as shipped.
    full_t <- bench(function() { orth_against(Q, W); invisible(NULL) })

    say("    n = %.0e ncv = %2d   floor %8.2f ms   full %8.2f ms   removable %.2fx  basis %.0f MB",
        n, ncv, 1000 * floor_t, 1000 * full_t, full_t / floor_t,
        n * ncv * 8 / 1024^2)
    rows[[length(rows) + 1L]] <- data.frame(
      n = n, ncv = ncv, floor_ms = 1000 * floor_t, full_ms = 1000 * full_t,
      removable_x = full_t / floor_t, basis_mb = n * ncv * 8 / 1024^2)
    rm(Q, W); gc(FALSE)
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, CSV, row.names = FALSE)

say("")
say("median removable factor inside orth_against: %.2fx", median(res$removable_x))
say("")
say("Reading: `removable` is the ratio of the shipped routine to the two products")
say("it cannot avoid. A compiled orth_against keeps the floor and can remove at")
say("most the rest, so this bounds the routine-level win with the algorithm fixed.")
