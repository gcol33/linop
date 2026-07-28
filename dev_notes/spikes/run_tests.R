## Run the suite and report counts plus the failing expectations only.
##
## The default reporter prints both operands of a failed comparison, which for
## the block sizes the solver suites use is tens of megabytes of matrix. This
## keeps the message and drops the operands.
##
##   Rscript dev_notes/spikes/run_tests.R                     # everything
##   Rscript dev_notes/spikes/run_tests.R test-solve-lsqr.R   # one file
files <- commandArgs(trailingOnly = TRUE)
devtools::load_all(".", quiet = TRUE)

res <- testthat::test_dir("tests/testthat", reporter = "silent",
                          stop_on_failure = FALSE, load_helpers = TRUE)

if (length(files)) {
  res <- Filter(function(t) basename(t$file) %in% files, res)
}

count <- function(cls) {
  sum(vapply(res, function(t) sum(vapply(t$results, inherits, logical(1), cls)),
             integer(1)))
}
cat(sprintf("tests %d | pass %d | fail %d | error %d | warn %d | skip %d\n",
            length(res), count("expectation_success"), count("expectation_failure"),
            count("expectation_error"), count("expectation_warning"),
            count("expectation_skip")))

for (t in res) {
  bad <- Filter(function(r) inherits(r, c("expectation_failure", "expectation_error")),
                t$results)
  if (!length(bad)) next
  cat(sprintf("\n[%s] %s\n", basename(t$file), t$test))
  for (r in bad) {
    msg <- utils::head(strsplit(conditionMessage(r), "\n")[[1]], 8)
    cat(paste0("    ", substr(msg, 1, 160), collapse = "\n"), "\n")
  }
}
