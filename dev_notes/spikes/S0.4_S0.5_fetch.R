## S0.4 + S0.5 -- fetch eigencore source, and sweep CRAN for prior art.
dest <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(dest)) dest <- "."
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
repos <- "https://cloud.r-project.org"

## ---------------------------------------------------- S0.4: eigencore source
src <- file.path(dest, "src")
dir.create(src, showWarnings = FALSE)
db <- available.packages(repos = repos)
for (p in c("eigencore", "Rlinsolve")) {
  v <- db[p, "Version"]
  url <- sprintf("%s/src/contrib/%s_%s.tar.gz", repos, p, v)
  f <- file.path(src, basename(url))
  cat("downloading", basename(url), "... ")
  r <- tryCatch({ download.file(url, f, mode = "wb", quiet = TRUE); "ok" },
                error = function(e) paste("FAIL", conditionMessage(e)))
  cat(r, "\n")
  if (r == "ok") { untar(f, exdir = src); cat("  extracted to", file.path(src, p), "\n") }
}

## ---------------------------------------------------- S0.5: CRAN prior art
cat("\n=== CRAN metadata sweep ===\n")
fields <- c("Package", "Version", "Title", "Description")
meta <- tools::CRAN_package_db()
cat("CRAN packages in db:", nrow(meta), "\n")
saveRDS(meta, file.path(dest, "cran_db.rds"))

hit <- function(pat, where = c("Title", "Description")) {
  txt <- do.call(paste, c(lapply(where, function(w) meta[[w]]), sep = " || "))
  i <- grep(pat, txt, ignore.case = TRUE, perl = TRUE)
  if (!length(i)) return(data.frame(Package = character(), Version = character(), Title = character()))
  data.frame(Package = meta$Package[i], Version = meta$Version[i],
             Title = gsub("\\s+", " ", meta$Title[i]), stringsAsFactors = FALSE)
}

terms <- list(
  MINRES        = "\\bMINRES\\b|minimum residual method",
  LSQR          = "\\bLSQR\\b",
  LSMR          = "\\bLSMR\\b",
  GMRES         = "\\bGMRES\\b",
  BiCGSTAB      = "\\bBiCGSTAB\\b|bi-?conjugate gradient",
  Krylov        = "\\bKrylov\\b",
  matrix_free   = "matrix[- ]free",
  Lanczos       = "\\bLanczos\\b",
  Arnoldi       = "\\bArnoldi\\b",
  preconditioner= "precondition",
  linear_operat = "linear operator",
  Chebfun_like  = "\\bChebfun\\b|\\bApproxFun\\b|spectral method|Chebyshev.*(approximation|spectral)",
  Hilbert_space = "Hilbert space",
  fun_space     = "function space|functional analysis",
  operator_disc = "discreti[sz]ation of.*operator|operator.*discreti[sz]ation"
)
res <- list()
for (nm in names(terms)) {
  h <- hit(terms[[nm]])
  res[[nm]] <- h
  cat(sprintf("\n--- %-14s (%d hits) ---\n", nm, nrow(h)))
  if (nrow(h)) print(utils::head(h[order(h$Package), ], 25), row.names = FALSE)
}
saveRDS(res, file.path(dest, "S0.5_hits.rds"))

## Numerical Mathematics task view
cat("\n=== Numerical Mathematics task view ===\n")
tv <- tryCatch({
  u <- "https://cran.r-project.org/web/views/NumericalMathematics.ctv"
  readLines(u, warn = FALSE)
}, error = function(e) { cat("fetch failed:", conditionMessage(e), "\n"); character() })
if (length(tv)) {
  writeLines(tv, file.path(dest, "NumericalMathematics.ctv"))
  cat("saved ctv,", length(tv), "lines\n")
  pk <- unique(unlist(regmatches(tv, gregexpr('(?<=<pkg>)[^<]+', tv, perl = TRUE))))
  cat("packages listed in view:", length(pk), "\n")
  writeLines(sort(pk), file.path(dest, "ctv_packages.txt"))
}
cat("\ndone\n")
