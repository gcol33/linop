## S0.5 -- function-level search. Title/Description misses implementations that
## do not advertise the algorithm name, so query r-universe's function index.
dest <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(dest)) dest <- "."

api <- function(url) {
  txt <- tryCatch(paste(readLines(url, warn = FALSE), collapse = ""),
                  error = function(e) NA_character_)
  if (is.na(txt)) return(NULL)
  tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
}

has_jsonlite <- requireNamespace("jsonlite", quietly = TRUE)
cat("jsonlite available:", has_jsonlite, "\n\n")

terms <- c("minres", "lsqr", "lsmr", "gmres", "bicgstab", "lanczos", "arnoldi",
           "golub kahan", "conjugate gradient", "matrix-free", "chebfun")

if (has_jsonlite) {
  for (t in terms) {
    u <- sprintf("https://r-universe.dev/api/search?q=%s&limit=20",
                 utils::URLencode(t, reserved = TRUE))
    j <- api(u)
    cat(sprintf("=== %-20s ", t))
    if (is.null(j) || is.null(j$results) || !length(j$results)) { cat("(no results / API unavailable)\n"); next }
    r <- j$results
    cat(sprintf("(%s total)\n", if (!is.null(j$total)) j$total else nrow(r)))
    keep <- intersect(c("Package", "package", "_user", "Title", "title"), names(r))
    print(utils::head(r[keep], 12), row.names = FALSE)
    cat("\n")
  }
}

## Direct check of the specific comparators the metadata sweep surfaced.
cat("\n=== inspect specific packages ===\n")
repos <- "https://cloud.r-project.org"
db <- available.packages(repos = repos)
for (p in c("sanic", "gmresls", "pcg", "lazymatrix", "Rlinsolve", "RSpectra", "irlba", "eigencore")) {
  if (!p %in% rownames(db)) { cat(sprintf("%-11s not on CRAN\n", p)); next }
  cat(sprintf("\n--- %s %s ---\n", p, db[p, "Version"]))
  cat("  Imports:  ", db[p, "Imports"], "\n")
  cat("  LinkingTo:", db[p, "LinkingTo"], "\n")
  ## pull the reference manual index to see exported function names
  u <- sprintf("https://cran.r-project.org/web/packages/%s/%s.pdf", p, p)
  cat("  manual:   ", u, "\n")
}
