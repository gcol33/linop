## Run one test file. The default reporter prints both operands of a failed
## comparison, which for the block sizes the solver suites use is tens of
## megabytes of matrix; this keeps the message and drops the operands.
##
##   Rscript dev_notes/spikes/run_one.R test-solve-lsqr.R
files <- commandArgs(trailingOnly = TRUE)
devtools::load_all(".", quiet = TRUE)

for (f in files) {
  res <- testthat::test_file(file.path("tests/testthat", f), reporter = "silent",
                             package = "linop")
  count <- function(cls) {
    sum(vapply(res, function(t) sum(vapply(t$results, inherits, logical(1), cls)),
               integer(1)))
  }
  cat(sprintf("%-28s tests %d | pass %d | fail %d | error %d\n", f,
              length(res), count("expectation_success"),
              count("expectation_failure"), count("expectation_error")))
  for (t in res) {
    bad <- Filter(function(r) inherits(r, c("expectation_failure", "expectation_error")),
                  t$results)
    if (!length(bad)) next
    cat(sprintf("\n[%s] %s\n", f, t$test))
    for (r in bad) {
      msg <- utils::head(strsplit(conditionMessage(r), "\n")[[1]], 8)
      cat(paste0("    ", substr(msg, 1, 160), collapse = "\n"), "\n")
    }
  }
}
