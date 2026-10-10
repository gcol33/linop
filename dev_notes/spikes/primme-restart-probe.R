## Does method = "primme" remove the reference eigensolver's O(ncv^2 n)
## orthogonalisation and its O(ncv) basis, or only the basis copy that
## dev_notes/compile-ceiling-and-the-basis-that-was-copied.md priced at 1.42x?
##
## Two fixtures, both at n = 1e5 and 1e6, both through eigs() so the
## certificate's own applies are in every timing:
##
##   budget   the stencil of the compile-ceiling spike, largest algebraic, k = 4,
##            which does not converge at these sizes. Every engine spends the
##            same budget of applies, so the time and the orthogonalisation per
##            apply are compared at equal work on the operator.
##   well     the stencil minus a potential well of depth 1 and width 30, the
##            four lowest bound states, tol 1e-10. These converge, so this is
##            time to an answer.
##
## Each cell runs in its own Rscript process, one after another, so the peak
## working set ps reports belongs to that cell alone. The reference's basis is
## R memory and PRIMME's is C memory, and the peak working set sees both.
##
## From the repository root:
##   Rscript dev_notes/spikes/primme-restart-probe.R            (all cells)
##   Rscript dev_notes/spikes/primme-restart-probe.R <cell...>  (one cell, internal)

if (!file.exists("DESCRIPTION")) stop("run this from the repository root")
OUT_DIR <- file.path("dev_notes", "spikes", "results")
CSV <- file.path(OUT_DIR, "primme-restart-probe.csv")
RSCRIPT <- file.path(R.home("bin"), "Rscript")

CELLS <- expand.grid(
  engine = c("lanczos-20", "lanczos-80", "primme-default", "primme-80"),
  n = c(1e5, 1e6), fixture = c("budget", "well"), stringsAsFactors = FALSE)
BUDGET <- 400L

args <- commandArgs(trailingOnly = TRUE)

if (!length(args)) {
  if (file.exists(CSV)) file.remove(CSV)
  for (i in seq_len(nrow(CELLS))) {
    a <- c(CELLS$fixture[i], format(CELLS$n[i], scientific = FALSE), CELLS$engine[i])
    cat(sprintf("%-7s n = %.0e  %-15s ", a[1], CELLS$n[i], a[3]))
    flush.console()
    st <- system2(RSCRIPT, c(shQuote("dev_notes/spikes/primme-restart-probe.R"), a),
                  stdout = TRUE, stderr = TRUE)
    cat(tail(st, 1L), "\n")
  }
  res <- utils::read.csv(CSV)
  print(res, row.names = FALSE)
  writeLines(format(Sys.time()), file.path(OUT_DIR, "primme-restart-probe.done"))
  quit(save = "no")
}

fixture <- args[1]; n <- as.numeric(args[2]); engine <- args[3]
suppressMessages(devtools::load_all(".", quiet = TRUE))

stencil <- function(n, v = NULL) {
  force(n)
  function(X) {
    X <- as.matrix(X)
    Y <- 2 * X
    Y[-n, ] <- Y[-n, , drop = FALSE] - X[-1L, , drop = FALSE]
    Y[-1L, ] <- Y[-1L, , drop = FALSE] - X[-n, , drop = FALSE]
    if (!is.null(v)) Y <- Y - v * X
    Y
  }
}

if (fixture == "budget") {
  A <- linop(stencil(n), dim = c(n, n), properties = c(hermitian = TRUE))
  which <- "largest_algebraic"
  tol <- .Machine$double.eps
  maxit <- BUDGET
} else {
  v <- numeric(n)
  mid <- n %/% 2
  v[(mid - 14):(mid + 15)] <- 1
  A <- linop(stencil(n, v), dim = c(n, n), properties = c(hermitian = TRUE))
  which <- "smallest_algebraic"
  tol <- 1e-10
  maxit <- 5000L
}

method <- sub("-.*", "", engine)
ncv <- switch(engine, "lanczos-20" = 20L, "lanczos-80" = 80L,
              "primme-default" = NULL, "primme-80" = 80L)

invisible(gc())
t0 <- proc.time()[["elapsed"]]
fit <- eigs(A, k = 4L, which = which, tol = tol, maxit = maxit, ncv = ncv,
            method = method, seed = 1L)
secs <- proc.time()[["elapsed"]] - t0
peak <- ps::ps_memory_info()[["peak_wset"]] / 2^20

st <- fit$engine$stats
row <- data.frame(
  fixture = fixture, n = n, engine = engine, seconds = round(secs, 3),
  applies = fit$iterations, ms_per_apply = round(1000 * secs / fit$iterations, 3),
  nconv = fit$nconv, certificate = fit$certificate$overall,
  basis = if (is.null(st)) ncv else fit$engine$params$maxBasisSize,
  ortho_inner_per_apply = if (is.null(st)) NA else round(st$numOrthoInnerProds / st$numMatvecs, 2),
  ortho_share = if (is.null(st)) NA else round(st$timeOrtho / st$elapsedTime, 3),
  peak_mb = round(peak), value1 = sprintf("%.15g", fit$values[1]),
  stringsAsFactors = FALSE)
utils::write.table(row, CSV, sep = ",", row.names = FALSE,
                   col.names = !file.exists(CSV), append = file.exists(CSV))
cat(sprintf("%.2fs, %d applies, %d conv, peak %.0f MB\n", secs, fit$iterations,
            fit$nconv, peak))
