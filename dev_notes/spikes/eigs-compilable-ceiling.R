## How much of an eigs() run could compiled code remove at all?
##
## The question this answers is not "is R slow" but "what is the ceiling on
## rewriting linop in C". Self time in the profile divides three ways:
##
##   compiled   BLAS, FFT, colSums, the arithmetic primitives. Already C. A
##              rewrite issues the same calls and removes none of it.
##   linop R    self time inside linop's own closures. This is the whole of what
##              compiling linop could remove, and it is an over-estimate, since
##              part of it is primitive arithmetic that C still has to perform.
##   apply      the operator's own callback, which for a matrix-free operator is
##              the user's R function. Compiling linop cannot touch it; a C core
##              would still call back into R for every one.
##
## The ceiling is then 1 / (1 - linopR), and it is generous in linop's favour.
##
## From the repository root: Rscript dev_notes/spikes/eigs-compilable-ceiling.R

if (!file.exists("DESCRIPTION")) stop("run this from the repository root")

devtools::load_all(".", quiet = TRUE)

OUT_DIR <- file.path("dev_notes", "spikes", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
CSV <- file.path(OUT_DIR, "eigs-compilable-ceiling.csv")
LOG <- file.path(OUT_DIR, "eigs-compilable-ceiling-detail.txt")
cat("", file = LOG)

say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOG, append = TRUE)
  flush.console()
}

## Everything linop defines. Membership here is what makes a profile entry
## linop's own R code rather than base R's compiled internals.
LINOP_NAMES <- ls(asNamespace("linop"), all.names = TRUE)

## The operator's callback, named so it is distinguishable in the profile from
## both linop's code and base R's.
## `inner` is the closure the fixture builds and `user_operator_apply` the named
## wrapper around it, so the operator's own work lands under one of the two and
## never under base R's compiled internals.
USER_APPLY <- c("user_operator_apply", "inner", "f", "h$apply")

classify <- function(nm) {
  bare <- gsub('^"|"$', "", nm)
  if (bare %in% USER_APPLY) return("apply")
  if (bare %in% LINOP_NAMES) return("linopR")
  "compiled"
}

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
make_fft <- function(n) {
  cvec <- numeric(n); cvec[1L] <- 2; cvec[2L] <- -1; cvec[n] <- -1
  fc <- Re(stats::fft(cvec))
  force(n)
  function(X) {
    X <- as.matrix(X)
    Y <- X
    for (j in seq_len(ncol(X))) {
      Y[, j] <- Re(stats::fft(fc * stats::fft(X[, j]), inverse = TRUE)) / n
    }
    Y
  }
}

cell <- function(kind, n, ncv) {
  inner <- if (kind == "stencil") make_stencil(n) else make_fft(n)
  ## Named, so the profile separates the operator's own work from linop's.
  user_operator_apply <- function(X) inner(X)
  A <- linop(user_operator_apply, dim = c(n, n),
             properties = c(hermitian = TRUE))
  run <- function() eigs(A, k = 4L, ncv = ncv, maxit = ncv, tol = 0, seed = 1L)

  t0 <- Sys.time()
  run()
  one <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  reps <- max(1L, min(20L, ceiling(3 / max(one, 1e-3))))

  pf <- tempfile(fileext = ".Rprof")
  Rprof(pf, interval = 0.002, gc.profiling = FALSE)
  for (i in seq_len(reps)) run()
  Rprof(NULL)
  s <- summaryRprof(pf)
  unlink(pf)

  self <- s$by.self
  total <- s$sampling.time
  cls <- vapply(rownames(self), classify, character(1))
  agg <- tapply(self$self.time, cls, sum)
  get <- function(k) if (k %in% names(agg)) unname(agg[[k]]) else 0
  compiled <- get("compiled"); linopR <- get("linopR"); apply_t <- get("apply")
  acc <- compiled + linopR + apply_t

  say("")
  say("--- %s n = %.0e ncv = %d   %.2fs sampled over %d reps", kind, n, ncv, total, reps)
  say("    compiled %6.1f%%   linop R %6.1f%%   apply %6.1f%%   (accounted %.1f%% of sampled)",
      100 * compiled / acc, 100 * linopR / acc, 100 * apply_t / acc, 100 * acc / total)
  ord <- head(self[order(-self$self.time), , drop = FALSE], 14L)
  for (i in seq_len(nrow(ord))) {
    say("      %-10s %-26s %7.3f", classify(rownames(ord)[i]),
        rownames(ord)[i], ord$self.time[i])
  }
  ceiling_x <- acc / max(acc - linopR, 1e-9)
  say("    ceiling if every line of linop's R became C: %.2fx", ceiling_x)

  data.frame(kind = kind, n = n, ncv = ncv, reps = reps,
             sampled_seconds = total,
             compiled_frac = compiled / acc,
             linopR_frac = linopR / acc,
             apply_frac = apply_t / acc,
             compile_ceiling = ceiling_x,
             stringsAsFactors = FALSE)
}

grid <- expand.grid(ncv = c(20L, 80L), n = c(1e5, 1e6),
                    kind = c("stencil", "fft"), stringsAsFactors = FALSE)
rows <- list()
for (i in seq_len(nrow(grid))) {
  r <- cell(grid$kind[i], grid$n[i], grid$ncv[i])
  rows[[length(rows) + 1L]] <- r
  utils::write.csv(do.call(rbind, rows), CSV, row.names = FALSE)
  gc(FALSE)
}
res <- do.call(rbind, rows)

say("")
say("=============================================== summary")
say("%-8s %8s %4s %10s %9s %8s %9s", "kind", "n", "ncv",
    "compiled%", "linopR%", "apply%", "ceiling")
for (i in seq_len(nrow(res))) {
  say("%-8s %8.0e %4d %9.1f%% %8.1f%% %7.1f%% %8.2fx",
      res$kind[i], res$n[i], res$ncv[i], 100 * res$compiled_frac[i],
      100 * res$linopR_frac[i], 100 * res$apply_frac[i], res$compile_ceiling[i])
}
say("")
say("median ceiling on compiling all of linop: %.2fx", median(res$compile_ceiling))
writeLines(format(Sys.time()), file.path(OUT_DIR, "eigs-compilable-ceiling.done"))
