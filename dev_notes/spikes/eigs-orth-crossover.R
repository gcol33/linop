## What fraction of eigs() runtime is the orthogonalisation that linop.primme
## would replace?
##
## `dev_notes/eigs-svds-and-the-third-certificate.md` section 6 says the
## reference label is about storage and work rather than accuracy: a round holds
## at most `ncv` vectors and orthogonalises each new direction against all of
## them, O(ncv) storage and O(ncv^2 n) work, "which is the gap `linop.primme`
## closes in Phase 3". That is a cost claim and it has never been measured.
##
## Whatever PRIMME saves on this axis is bounded above by the fraction of
## runtime that goes into `orth_against()`. If that fraction is small at the
## sizes a matrix-free caller runs, then vendoring the C library buys a
## proportionally small amount, whatever the flop counts say -- because the
## orthogonalisation is two BLAS calls and the apply is interpreted R, so flops
## are not the currency.
##
## Not measured here, and it has to be said plainly: PRIMME's other axis is
## converging in fewer matvecs (JDQMR, GD+k, preconditioned inner solves). That
## cannot be measured without PRIMME and this script says nothing about it. What
## is bounded here is only the axis the note names.
##
## From the repository root: Rscript dev_notes/spikes/eigs-orth-crossover.R

if (!file.exists("DESCRIPTION")) {
  stop("run this from the repository root")
}

devtools::load_all(".", quiet = TRUE)

OUT_DIR <- file.path("dev_notes", "spikes", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
CSV <- file.path(OUT_DIR, "eigs-orth-crossover.csv")
LOG <- file.path(OUT_DIR, "eigs-orth-crossover-profiles.txt")

cat("", file = LOG)
say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOG, append = TRUE)
  flush.console()
}

## ------------------------------------------------------------ the two applies
##
## Both are matrix free and both are hermitian with a closed-form spectrum, so a
## profiled run can be checked against truth rather than assumed.
##
## They differ in cost per apply, which is the variable the decision turns on. A
## stencil is about the cheapest apply an operator can have, four flops an entry.
## The FFT route is a stand-in for any operator that does real work per apply,
## which is the usual reason an operator has no matrix in the first place.

## Dirichlet Laplacian, eigenvalues 4 sin^2(k pi / (2(n+1))).
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
stencil_values <- function(n, k) {
  lam <- 4 * sin(seq_len(n) * pi / (2 * (n + 1)))^2
  sort(lam, decreasing = TRUE)[seq_len(k)]
}

## Periodic Laplacian as a circulant, applied through the FFT. Eigenvalues are
## the DFT of the first column, which plan section 10 already lists as a
## closed-form fixture.
make_fft <- function(n) {
  cvec <- numeric(n)
  cvec[1L] <- 2
  cvec[2L] <- -1
  cvec[n] <- -1
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
fft_values <- function(n, k) {
  cvec <- numeric(n)
  cvec[1L] <- 2
  cvec[2L] <- -1
  cvec[n] <- -1
  sort(Re(stats::fft(cvec)), decreasing = TRUE)[seq_len(k)]
}

## An independently built matrix for each, used only by validate_kind(). Written
## from the definition rather than from the apply, so checking one against the
## other is not circular.
ref_stencil <- function(n) {
  M <- matrix(0, n, n)
  diag(M) <- 2
  i <- seq_len(n - 1L)
  M[cbind(i, i + 1L)] <- -1
  M[cbind(i + 1L, i)] <- -1
  M
}
ref_fft <- function(n) {
  cvec <- numeric(n)
  cvec[1L] <- 2
  cvec[2L] <- -1
  cvec[n] <- -1
  outer(seq_len(n), seq_len(n), function(i, j) cvec[((i - j) %% n) + 1L])
}

KINDS <- list(
  stencil = list(make = make_stencil, truth = stencil_values, ref = ref_stencil,
                 note = "tridiagonal stencil, ~4 flops an entry"),
  fft     = list(make = make_fft, truth = fft_values, ref = ref_fft,
                 note = "circulant through two FFTs, O(n log n)")
)

## ------------------------------------------------------------------ the grid
##
## maxit = ncv runs exactly one round, so the basis grows 1 .. ncv and the
## orthogonalisation cost is the sum the note describes rather than whatever a
## restart schedule happened to produce. tol = 0 stops anything converging or
## locking early, so every cell does the same amount of work and the fractions
## are comparable across the grid. Neither is a cap on a reported quantity; both
## are the control that makes the cells mean the same thing.

NS   <- c(1e4, 1e5, 1e6)
NCVS <- c(20L, 40L, 80L)
K    <- 4L
TARGET_SECONDS <- 3     # of profiling per cell, which sets the repetitions

profile_cell <- function(kind, n, ncv) {
  spec <- KINDS[[kind]]
  A <- linop(spec$make(n), dim = c(n, n),
             properties = c(hermitian = TRUE))

  run <- function() eigs(A, k = K, ncv = ncv, maxit = ncv, tol = 0, seed = 1L)

  t0 <- Sys.time()
  warm <- run()
  one <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  reps <- max(1L, min(20L, ceiling(TARGET_SECONDS / max(one, 1e-3))))

  pf <- tempfile(fileext = ".Rprof")
  Rprof(pf, interval = 0.002, gc.profiling = FALSE)
  for (i in seq_len(reps)) run()
  Rprof(NULL)
  s <- summaryRprof(pf)
  unlink(pf)

  by <- s$by.total
  total <- s$sampling.time
  grab <- function(nm) {
    row <- paste0("\"", nm, "\"")
    if (row %in% rownames(by)) by[row, "total.time"] else 0
  }
  orth <- grab("orth_against")

  ## What the profile actually saw, so the extraction above can be checked
  ## rather than trusted.
  say("")
  say("--- %s n = %.0e ncv = %d  (%d reps, %.2fs sampled)", kind, n, ncv, reps, total)
  top <- head(by[order(-by$total.time), c("total.time", "total.pct", "self.time"), drop = FALSE], 12L)
  for (i in seq_len(nrow(top))) {
    say("    %-40s %8.3f %6.1f%% %8.3f", rownames(top)[i],
        top$total.time[i], top$total.pct[i], top$self.time[i])
  }

  ## No accuracy column here on purpose. tol = 0 means these cells never
  ## converge, which is what makes them comparable; reading an error off them
  ## would be reporting how far one round happened to get. What has to be right
  ## for the profile to mean anything is that the apply is the operator its
  ## closed form describes, and that is checked once per kind by validate_kind()
  ## on a size where a converged run is cheap.
  data.frame(
    kind = kind, n = n, ncv = ncv, k = K,
    reps = reps,
    seconds_run = one,
    sampled_seconds = total,
    orth_seconds = orth,
    orth_fraction = if (total > 0) orth / total else NA_real_,
    basis_mb = n * ncv * 8 / 2^20,
    iterations = warm$iterations,
    stringsAsFactors = FALSE
  )
}

## ------------------------------------------------------------- the apply is
## the operator it claims to be
##
## Two exact checks at a small n, neither involving an iteration. The apply code
## does not depend on n, so a small n checks it.
##
## Not a converged eigs() run, which is what the first draft did: both fixtures
## are Laplacians whose extreme eigenvalues are separated by O(1/n^2), 2.5e-7 at
## n = 2000, so a run that has not resolved them is reporting clustering rather
## than a wrong apply. That check failed the fixture for a property the fixture
## does not have to have, and told us nothing about the apply either way.
validate_kind <- function(kind, n = 300L) {
  spec <- KINDS[[kind]]
  M <- spec$ref(n)
  f <- spec$make(n)

  ## 1. The callback applies the matrix the definition describes.
  X <- with_preserved_seed(7L, matrix(stats::rnorm(n * 3L), n, 3L))
  got <- f(X)
  want <- M %*% X
  apply_err <- max(abs(got - want)) / max(abs(want))

  ## 2. The closed form is that matrix's spectrum.
  lam <- sort(eigen(M, symmetric = TRUE, only.values = TRUE)$values,
              decreasing = TRUE)[seq_len(K)]
  truth <- spec$truth(n, K)
  value_err <- max(abs(lam - truth) / abs(truth))

  say("validate %-8s n = %d: apply %.2e, closed form %.2e  (%s)",
      kind, n, apply_err, value_err, spec$note)
  if (!is.finite(apply_err) || apply_err > 1e-10) {
    stop(sprintf("%s apply disagrees with its definition: %.3e", kind, apply_err))
  }
  if (!is.finite(value_err) || value_err > 1e-10) {
    stop(sprintf("%s closed form disagrees with eigen(): %.3e", kind, value_err))
  }
  invisible(c(apply = apply_err, values = value_err))
}

## ------------------------------------------------------------------- the run

for (kind in names(KINDS)) validate_kind(kind)
say("")

rows <- list()
for (kind in names(KINDS)) {
  for (n in NS) {
    for (ncv in NCVS) {
      r <- tryCatch(profile_cell(kind, n, ncv),
                    error = function(e) {
                      say("!!! %s n = %.0e ncv = %d failed: %s", kind, n, ncv,
                          conditionMessage(e))
                      NULL
                    })
      if (is.null(r)) next
      rows[[length(rows) + 1L]] <- r
      ## Written every cell, so a run that is killed still leaves what it reached.
      utils::write.csv(do.call(rbind, rows), CSV, row.names = FALSE)
      say("  -> orth %.1f%% of %.2fs, basis %.0f Mb",
          100 * r$orth_fraction, r$sampled_seconds, r$basis_mb)
      gc(FALSE)
    }
  }
}

res <- do.call(rbind, rows)
utils::write.csv(res, CSV, row.names = FALSE)

say("")
say("================================ summary")
say("%-8s %8s %5s %9s %9s", "kind", "n", "ncv", "orth %", "basis Mb")
for (i in seq_len(nrow(res))) {
  say("%-8s %8.0e %5d %8.1f%% %9.0f",
      res$kind[i], res$n[i], res$ncv[i], 100 * res$orth_fraction[i],
      res$basis_mb[i])
}
say("")
say("done: %s", CSV)
writeLines(format(Sys.time()), file.path(OUT_DIR, "eigs-orth-crossover.done"))
