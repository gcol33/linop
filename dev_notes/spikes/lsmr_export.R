## Write the fixtures and this implementation's iterates for the SciPy
## cross-check, and read the comparison back.
##
##   Rscript dev_notes/spikes/lsmr_export.R write
##   python  dev_notes/spikes/lsmr_scipy.py
##   Rscript dev_notes/spikes/lsmr_export.R compare

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

dir <- file.path(tempdir(), "lsmr-xcheck")
dir <- Sys.getenv("LSMR_XDIR", dir)
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
cat("dir:", dir, "\n")

## Well conditioned, badly conditioned, and rank deficient, each at several step
## counts, so the comparison covers the regime where the recurrence is still
## reproducible and the regime where it is not.
cases <- list(
  list(name = "kappa1e2",  kappa = 1e2,  m = 60, n = 20, seed = 1),
  list(name = "kappa1e6",  kappa = 1e6,  m = 60, n = 20, seed = 2),
  list(name = "kappa1e10", kappa = 1e10, m = 60, n = 20, seed = 3),
  list(name = "tall",      kappa = 1e4,  m = 200, n = 30, seed = 4))
steps <- c(1L, 2L, 4L, 8L, 16L, 32L, 64L)

build <- function(cs) {
  sg <- exp(seq(log(1), log(1 / cs$kappa), length.out = cs$n))
  f <- lsq_prescribed(cs$m, cs$n, sg, seed = cs$seed)
  set.seed(5000L + cs$seed)
  list(M = f$A, b = matrix(stats::rnorm(cs$m), cs$m, 1))
}

action <- commandArgs(trailingOnly = TRUE)[1L]

if (identical(action, "write")) {
  for (cs in cases) {
    d <- build(cs)
    utils::write.table(d$M, file.path(dir, paste0(cs$name, "-A.txt")),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(d$b, file.path(dir, paste0(cs$name, "-b.txt")),
                       row.names = FALSE, col.names = FALSE)
    for (k in steps) {
      x <- linop:::lsmr_solve(linop(d$M), d$b, tol = 0, maxit = k)$x
      utils::write.table(x, file.path(dir, sprintf("%s-mine-%d.txt", cs$name, k)),
                         row.names = FALSE, col.names = FALSE)
    }
  }
  writeLines(as.character(steps), file.path(dir, "steps.txt"))
  writeLines(vapply(cases, function(cs) cs$name, character(1)),
             file.path(dir, "cases.txt"))
  cat("written\n")
}

if (identical(action, "compare")) {
  cat(sprintf("%-10s %6s  %12s  %12s\n", "case", "steps", "rel gap", "||x||"))
  for (cs in cases) {
    for (k in steps) {
      fm <- file.path(dir, sprintf("%s-mine-%d.txt", cs$name, k))
      fs <- file.path(dir, sprintf("%s-scipy-%d.txt", cs$name, k))
      if (!file.exists(fs)) next
      mine <- as.matrix(utils::read.table(fm))
      theirs <- as.matrix(utils::read.table(fs))
      nx <- max(Mod(theirs))
      cat(sprintf("%-10s %6d  %12.3e  %12.3e\n", cs$name, k,
                  if (nx > 0) max(Mod(mine - theirs)) / nx else max(Mod(mine - theirs)),
                  nx))
    }
  }
}
