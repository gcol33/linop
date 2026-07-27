pkg <- "C:/Users/Gilles Colling/Documents/dev/linop"
suppressMessages({ library(devtools); library(testthat) })
devtools::load_all(pkg, quiet = TRUE)

cat("################ GATE 1 ################\n")
cat("R:", R.version.string, "\n")
cat("linop:", as.character(utils::packageVersion("linop")), "\n\n")

## ---------------------------------------------------------- install cost ----
cat("=== install cost (tools::package_dependencies recursive) ===\n")
d <- read.dcf(file.path(pkg, "DESCRIPTION"))
cat("Depends:  ", if ("Depends" %in% colnames(d)) d[1, "Depends"] else "-", "\n")
cat("Imports:  ", if ("Imports" %in% colnames(d)) d[1, "Imports"] else "(none)", "\n")
cat("LinkingTo:", if ("LinkingTo" %in% colnames(d)) d[1, "LinkingTo"] else "(none)", "\n")
cat("Suggests: ", if ("Suggests" %in% colnames(d)) d[1, "Suggests"] else "-", "\n")

db <- available.packages(repos = "https://cloud.r-project.org")
imports <- if ("Imports" %in% colnames(d)) trimws(strsplit(d[1, "Imports"], ",")[[1]]) else character(0)
imports <- sub("\\s*\\(.*\\)$", "", imports)
imports <- imports[nzchar(imports)]
base_pkgs <- rownames(installed.packages(priority = c("base", "recommended")))
if (length(imports)) {
  rec <- unique(unlist(tools::package_dependencies(imports, db = db, recursive = TRUE,
                        which = c("Depends", "Imports", "LinkingTo"))))
  cat("recursive hard deps:", length(unique(c(imports, rec))), "\n")
} else {
  cat("recursive hard deps: 0  (no Imports at all)\n")
}
cat("compiled code:", if (dir.exists(file.path(pkg, "src"))) "YES" else "none", "\n\n")

## ------------------------------------------------------- full fuzzer run ----
cat("=== expression fuzzer, full Gate 1 count ===\n")
Sys.setenv(LINOP_FUZZ_N = "10000")
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(file.path(pkg, "tests/testthat/test-fuzzer.R"),
                           reporter = "silent")
el <- proc.time()[["elapsed"]] - t0
df <- as.data.frame(res)
cat(sprintf("trees: 10000   failed: %d   errors: %d   time: %.1fs\n",
            sum(df$failed), sum(df$error), el))
Sys.unsetenv("LINOP_FUZZ_N")

## ---------------------------------------------------------- full suite ------
cat("\n=== full test suite ===\n")
t0 <- proc.time()[["elapsed"]]
res2 <- testthat::test_dir(file.path(pkg, "tests/testthat"), reporter = "silent",
                           stop_on_failure = FALSE)
el2 <- proc.time()[["elapsed"]] - t0
df2 <- as.data.frame(res2)
cat(sprintf("files: %d  tests: %d  passed: %d  failed: %d  errors: %d  skipped: %d  time: %.1fs\n",
            length(unique(df2$file)), nrow(df2), sum(df2$passed),
            sum(df2$failed), sum(df2$error), sum(df2$skipped), el2))

by_file <- aggregate(cbind(passed, failed, error) ~ file, data = df2, FUN = sum)
print(by_file, row.names = FALSE)

cat("\n=== node coverage ===\n")
cat("registered nodes:", paste(linop_nodes(), collapse = ", "), "\n")
