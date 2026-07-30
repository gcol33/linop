## The basis slice at eigs.R:455, and whether pure R can avoid it.
##
## `orth_against(V[, seq_len(i), drop = FALSE], W)` copies i columns of the live
## basis before BLAS can see them, because base R has no matrix view. The copy is
## one extra pass over the same memory each of the two products already streams,
## so it is a third of the routine rather than a rounding error, and it is charged
## to orth_against's frame because lazy evaluation forces the promise there.
##
## Three routes are timed at each i, summed over a whole round the way the solver
## runs it:
##
##   slice   what ships. copy i columns, two products over i columns.
##   zeros   orthogonalise against the whole preallocated basis. The unused
##           columns are zero, so their coefficients are zero and the result is
##           identical; no copy, but both products run over ncv columns instead
##           of i.
##   floor   the two products over i columns with no copy at all. Not reachable
##           from R; this is what a compiled routine with a leading-dimension
##           argument would issue, and it bounds the win.
##
## `zeros` is checked for bitwise agreement with `slice`, since a route that
## changes the answer is not a candidate whatever it costs.
##
## From the repository root: Rscript dev_notes/spikes/orth-basis-slice.R

if (!file.exists("DESCRIPTION")) stop("run this from the repository root")

devtools::load_all(".", quiet = TRUE)

OUT_DIR <- file.path("dev_notes", "spikes", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT_DIR, "orth-basis-slice.txt")
CSV <- file.path(OUT_DIR, "orth-basis-slice.csv")
cat("", file = LOG)

say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOG, append = TRUE)
  flush.console()
}

bench <- function(f, seconds = 1.5) {
  f()
  t0 <- Sys.time(); reps <- 0L
  repeat {
    f(); reps <- reps + 1L
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) >= seconds) break
  }
  as.numeric(difftime(Sys.time(), t0, units = "secs")) / reps
}

rows <- list()
for (cfg in list(c(1e5, 80), c(1e6, 80))) {
  n <- cfg[1]; ncv <- as.integer(cfg[2])
  say("")
  say("=== n = %.0e  ncv = %d   full basis %.0f MB", n, ncv, n * ncv * 8 / 1024^2)

  V <- matrix(0, n, ncv)
  Vfull <- qr.Q(qr(matrix(stats::rnorm(n * ncv), n, ncv)))
  W0 <- matrix(stats::rnorm(n), n, 1L)

  ## i sampled across the round rather than only at the top, since the copy grows
  ## with i and a single i would misprice the sum.
  is <- unique(as.integer(round(seq(ncv / 8, ncv, length.out = 4))))
  tot <- c(slice = 0, zeros = 0, floor = 0)

  for (i in is) {
    V[] <- 0
    V[, seq_len(i)] <- Vfull[, seq_len(i)]

    t_slice <- bench(function() {
      Q <- V[, seq_len(i), drop = FALSE]
      C <- crossprod(Q, W0); W0 - Q %*% C
      invisible(NULL)
    })
    t_zeros <- bench(function() {
      C <- crossprod(V, W0); W0 - V %*% C
      invisible(NULL)
    })
    Qi <- V[, seq_len(i), drop = FALSE]
    t_floor <- bench(function() {
      C <- crossprod(Qi, W0); W0 - Qi %*% C
      invisible(NULL)
    })

    ## Identical answer, or the route is not a candidate.
    a <- { Q <- V[, seq_len(i), drop = FALSE]; W0 - Q %*% crossprod(Q, W0) }
    b <- W0 - V %*% crossprod(V, W0)
    ok <- identical(a, b)

    say("    i = %3d   slice %8.2f ms   zeros %8.2f ms   floor %8.2f ms   zeros bitwise: %s",
        i, 1000 * t_slice, 1000 * t_zeros, 1000 * t_floor, ok)
    tot["slice"] <- tot["slice"] + t_slice
    tot["zeros"] <- tot["zeros"] + t_zeros
    tot["floor"] <- tot["floor"] + t_floor
    rows[[length(rows) + 1L]] <- data.frame(
      n = n, ncv = ncv, i = i, slice_ms = 1000 * t_slice,
      zeros_ms = 1000 * t_zeros, floor_ms = 1000 * t_floor, zeros_bitwise = ok)
    rm(Qi, a, b); gc(FALSE)
  }

  say("    summed over the sampled round:")
  say("      slice  %8.2f ms   1.00x  (ships)", 1000 * tot["slice"])
  say("      zeros  %8.2f ms   %.2fx  (pure R, no copy)",
      1000 * tot["zeros"], tot["slice"] / tot["zeros"])
  say("      floor  %8.2f ms   %.2fx  (compiled, leading dimension)",
      1000 * tot["floor"], tot["slice"] / tot["floor"])
  rm(V, Vfull, W0); gc(FALSE)
}

utils::write.csv(do.call(rbind, rows), CSV, row.names = FALSE)
