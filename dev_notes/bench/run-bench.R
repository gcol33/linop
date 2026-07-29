## The benchmark harness, end to end. From the repository root:
##
##     Rscript dev_notes/bench/run-bench.R
##
## Writes one CSV per module under dev_notes/bench/results/, the machine that
## produced them to results/environment.txt, and the tables to RESULTS.md. All of
## it is committed: a benchmark whose numbers live only in a console is not a
## result anyone can check.

if (!file.exists("DESCRIPTION")) {
  stop("run this from the repository root: Rscript dev_notes/bench/run-bench.R")
}

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")
source("dev_notes/bench/bench-common.R")
source("dev_notes/bench/bench-solvers.R")
source("dev_notes/bench/bench-blocks.R")
source("dev_notes/bench/bench-matrixfree.R")
source("dev_notes/bench/bench-spectral.R")

bench_init()
started <- Sys.time()

env <- bench_environment()
writeLines(sprintf("%-20s %s", env$field, env$value),
           file.path(BENCH_RESULTS, "environment.txt"))
cat("environment\n")
cat(sprintf("  %-20s %s\n", env$field, env$value), sep = "")

cat("\nsolvers\n")
solvers <- bench_solvers()
bench_write(solvers, "solvers")

cat("\nblocks\n")
blocks <- bench_blocks()
bench_write(blocks, "blocks")

cat("\nmatrix free against a matrix\n")
matrixfree <- bench_matrixfree()
bench_write(matrixfree, "matrixfree")

cat("\nspectral\n")
spectral <- bench_spectral()
bench_write(spectral, "spectral")

elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))

## ------------------------------------------------------------------ RESULTS

section <- function(title, prose, df, cols) {
  c(sprintf("## %s", title), "", gsub(" \n", "\n", prose, fixed = TRUE), "",
    bench_md_table(df[, cols, drop = FALSE]), "")
}

md <- c(
  "# Benchmark results",
  "",
  "Produced by `dev_notes/bench/run-bench.R`, which regenerates this file and",
  "everything under `results/`. Gate 2's last line asks for a harness that runs",
  "end to end with committed results, and these are they.",
  "",
  "Two units. Operator applies is the primary one: it is what a matrix-free",
  "caller pays and it is the same number on every machine. Wall time is",
  "secondary, because it is a fact about the machine below and about no other.",
  "Applies also do not price everything, which the GMRES and eigensolver rows",
  "make visible: both orthogonalise against a stored basis, and that arithmetic",
  "is invisible to a count of applies.",
  "",
  "```",
  sprintf("%-20s %s", env$field, env$value),
  sprintf("%-20s %.1f minutes", "harness runtime", elapsed),
  "```",
  "",
  section(
    "Seven methods, eight operators",
    paste(
      "Same tolerance, same budget, same right-hand side and the same zero start",
      "for every method on an operator. The roster per fixture is what the",
      "operator supplies rather than a ranking: `convdiff_1d` is nonsymmetric and",
      "offers the forward action alone, which are two separate exclusions, CG and",
      "MINRES for the missing symmetry and LSQR and LSMR for the missing adjoint.",
      "Those two are the only methods that require an adjoint. CG is also absent",
      "wherever the operator is indefinite. The `norm2()` row",
      "under each fixture is the operator-norm estimate the certificate's",
      "arithmetic floor rests on, priced on its own and already included in every",
      "solve row above it.",
      "\n\nForward error is measured against the closed form in the `truth`",
      "column of `results/solvers.csv` and is not something the certificate",
      "claims: `forward error` reports `not_checked` in every one of these runs,",
      "for the reason in section 6.1."),
    solvers,
    c("fixture", "kappa", "method", "iterations", "applies", "seconds",
      "backward_error", "forward_error", "overall")),
  section(
    "Right-hand sides in lockstep",
    paste(
      "k columns through one solve against k columns through k solves. The block",
      "takes one apply per step where the separate route takes k, and pays for it",
      "in column-applies, because every column keeps being carried until the last",
      "one converges. The final column is how far apart the two answers are."),
    blocks,
    c("fixture", "k", "route", "applies", "column_applies", "seconds",
      "agreement_vs_lockstep")),
  section(
    "The size that stops being a problem",
    paste(
      "tridiag(-1, 3, -1) solved through its action, through a dense matrix and",
      "through a sparse one. The shift keeps the spectrum in (1, 5) at every n,",
      "so the iteration count is flat down the ladder and what moves is size.",
      "The sparse column is the fastest route at every size here, which is what a",
      "banded operator with an O(n) direct solve should look like.",
      "\n\n`object_mb` is `object.size()` of what the route holds: the closure and",
      "its counter for the operator, the matrix for the other two."),
    matrixfree,
    c("n", "route", "seconds", "iterations", "applies", "object_mb",
      "forward_error", "note")),
  section(
    "eigs and svds",
    paste(
      "Both verbs against closed-form spectra, at three subspace sizes, with a",
      "dense factorisation of the same operator underneath. The dense routes",
      "return the whole spectrum where the two verbs return k pairs, so the",
      "seconds column compares different deliverables and the value error column",
      "compares the same k values.",
      "\n\nRead `nconv` before anything else in this table. Both ends of",
      "`laplacian_1d(400)` have a relative gap of 4.6e-5, and at ncv = 20 and 40",
      "the run stalls with nothing converged. Raising `maxit` from 1000 to 5000",
      "or loosening `tol` from 1e-10 to 1e-8 changes neither the iteration count",
      "nor the value error there, so the subspace is the binding knob and the",
      "stall detector is reporting that rather than giving up early. A row with",
      "`nconv` below k is returning its best unconverged pairs, and its value",
      "error is the size of that.",
      "\n\nThe shift-invert row counts every inner MINRES step through the same",
      "counter, so its applies column is the cost of the transformation and not",
      "of the outer recurrence."),
    spectral,
    c("problem", "verb", "which", "route", "nconv", "applies", "seconds",
      "backward_error", "value_error")))

writeLines(md, file.path(BENCH_DIR, "RESULTS.md"))
cat(sprintf("\nwrote %s\n", file.path(BENCH_DIR, "RESULTS.md")))
cat(sprintf("harness runtime %.1f minutes\n", elapsed))
