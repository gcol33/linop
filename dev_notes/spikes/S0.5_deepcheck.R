## S0.5 -- resolve the ambiguous hits. Download candidate sources and read their
## NAMESPACE exports, so the "absent from CRAN" claim rests on exports, not blurbs.
dest <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(dest)) dest <- "."
src  <- file.path(dest, "prior"); dir.create(src, recursive = TRUE, showWarnings = FALSE)
repos <- "https://cloud.r-project.org"
db <- available.packages(repos = repos)

cand <- c("pracma", "sanic", "gmresls", "pcg", "lazymatrix", "TSSS", "diffcp",
          "withinr", "amatrix", "optR", "SparseM", "svd", "rARPACK")

cat("=== on CRAN? ===\n")
for (p in cand) cat(sprintf("  %-11s %s\n", p, if (p %in% rownames(db)) db[p, "Version"] else "NOT ON CRAN"))

get_exports <- function(p) {
  if (!p %in% rownames(db)) return(NULL)
  v <- db[p, "Version"]
  f <- file.path(src, sprintf("%s_%s.tar.gz", p, v))
  if (!file.exists(f)) {
    ok <- tryCatch({ download.file(sprintf("%s/src/contrib/%s_%s.tar.gz", repos, p, v),
                                   f, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (!ok) return(NULL)
  }
  nsfile <- file.path(p, "NAMESPACE")
  tryCatch({
    untar(f, files = nsfile, exdir = src)
    readLines(file.path(src, nsfile), warn = FALSE)
  }, error = function(e) NULL)
}

pats <- c(minres = "minres", lsqr = "lsqr", lsmr = "lsmr", gmres = "gmres",
          bicgstab = "bicgstab", cg = "\\bcg\\b|conjgrad", lanczos = "lanczos",
          arnoldi = "arnoldi")

cat("\n=== exported names matching Krylov algorithm tokens ===\n")
for (p in cand) {
  ns <- get_exports(p)
  if (is.null(ns)) { cat(sprintf("\n--- %-11s (no NAMESPACE retrieved)\n", p)); next }
  exp <- unlist(regmatches(ns, gregexpr('(?<=export\\()[^)]+', ns, perl = TRUE)))
  exp <- unlist(strsplit(gsub('["\\s]', "", exp, perl = TRUE), ","))
  exp <- exp[nzchar(exp)]
  hits <- character()
  for (nm in names(pats)) {
    h <- grep(pats[[nm]], exp, ignore.case = TRUE, perl = TRUE, value = TRUE)
    if (length(h)) hits <- c(hits, sprintf("%s: %s", nm, paste(h, collapse = " ")))
  }
  cat(sprintf("\n--- %-11s (%d exports) ---\n", p, length(exp)))
  if (length(hits)) cat("   ", paste(hits, collapse = "\n    "), "\n") else cat("    no Krylov-token exports\n")
  ## surface any S4 class exports, for the operator-abstraction question
  cls <- unlist(regmatches(ns, gregexpr('(?<=exportClasses\\()[^)]+', ns, perl = TRUE)))
  if (length(cls)) cat("    exportClasses:", paste(cls, collapse = ", "), "\n")
}
