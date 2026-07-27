## S0.2 -- R callback overhead for a `fun` leaf.
## Question (plan section 9): is the R round trip expensive relative to the matvec?
## If yes, block widths become mandatory and some solver designs change.
##
## Timing uses microbenchmark's high-resolution counter. proc.time() was tried
## first and is unusable here: ~1 ms resolution on Windows against operations of
## a few hundred nanoseconds.

set.seed(20260727)
suppressMessages(library(microbenchmark))
cat("R:", R.version.string, "\n")
cat("LAPACK:", tryCatch(La_version(), error = function(e) "unknown"), "\n")
cat("platform:", R.version$platform, "\n")
cat("microbenchmark:", as.character(packageVersion("microbenchmark")), "\n\n")

## median seconds for f(X); microbenchmark reports nanoseconds
tmed <- function(f, X, reps) {
  f(X)
  median(microbenchmark::microbenchmark(f(X), times = reps, unit = "ns")$time) / 1e9
}

## ---------------------------------------------------------------------------
## A. pure round trip: closures that do no arithmetic
## ---------------------------------------------------------------------------
cat("=== A. bare round trip, no arithmetic ===\n")
X1 <- matrix(0, 4, 1); v1 <- numeric(4)
noop  <- function(X) X
wrap2 <- function(X) noop(X)
wrap3 <- function(X) wrap2(X)
dimchk  <- function(X) if (is.null(dim(X))) as.matrix(X) else X
asmat   <- function(X) as.matrix(X)
dimset  <- function(X) { dim(X) <- c(length(X), 1L); X }

for (nm in c("noop", "wrap2", "wrap3", "dimchk")) {
  cat(sprintf("  %-32s %7.1f ns\n",
      switch(nm, noop = "1 closure call (identity)", wrap2 = "2 nested (leaf -> user fn)",
             wrap3 = "3 nested calls", dimchk = "dim() check, already a matrix"),
      tmed(get(nm), X1, 3e4) * 1e9))
}
cat(sprintf("  %-32s %7.1f ns\n", "as.matrix() on length-4 vector", tmed(asmat,  v1, 3e4) * 1e9))
cat(sprintf("  %-32s %7.1f ns\n", "dim<-() on length-4 vector",     tmed(dimset, v1, 3e4) * 1e9))
CALL <- tmed(wrap2, X1, 3e4)

## ---------------------------------------------------------------------------
## B. realistic operators: direct vs wrapped in a `fun` leaf
## ---------------------------------------------------------------------------
mk <- list(
  diag = function(n) { d <- runif(n) + 1; force(d); function(X) d * X },
  tri  = function(n) function(X) {
    m <- nrow(X)
    Y <- 2 * X
    Y[-1, ] <- Y[-1, ] - X[-m, ]
    Y[-m, ] <- Y[-m, ] - X[-1, ]
    Y
  },
  dense = function(n) { M <- matrix(rnorm(n * n), n, n); force(M); function(X) M %*% X }
)

grid <- list(
  list(n = 1e3, bs = c(1, 4, 16, 64), kinds = c("diag", "tri", "dense"), reps = 500),
  list(n = 1e5, bs = c(1, 4, 16, 64), kinds = c("diag", "tri"),          reps = 50),
  list(n = 1e7, bs = c(1, 4, 16),     kinds = c("diag", "tri"),          reps = 7)
)

rows <- list()
for (g in grid) for (kind in g$kinds) {
  f_user <- mk[[kind]](g$n)
  leaf   <- function(X) f_user(X)
  for (b in g$bs) {
    X  <- matrix(rnorm(g$n * b), g$n, b)
    td <- tmed(f_user, X, g$reps)
    tw <- tmed(leaf,   X, g$reps)
    rows[[length(rows) + 1L]] <- data.frame(
      n = g$n, kind = kind, b = b,
      direct_us  = td * 1e6,
      wrapped_us = tw * 1e6,
      per_col_us = tw * 1e6 / b,
      call_pct   = 100 * CALL / tw)
    rm(X); invisible(gc(verbose = FALSE))
  }
  rm(f_user, leaf); invisible(gc(verbose = FALSE))
}
res <- do.call(rbind, rows)
cat("\n=== B. apply cost; call_pct = round trip as % of one wrapped apply ===\n")
print(res, row.names = FALSE, digits = 4)

## ---------------------------------------------------------------------------
## C. does block width amortise anything?
## ---------------------------------------------------------------------------
cat("\n=== C. per-column cost vs block width ===\n")
for (k in c("diag", "tri", "dense")) for (nn in sort(unique(res$n[res$kind == k]))) {
  s <- res[res$kind == k & res$n == nn, ]
  if (nrow(s) < 2) next
  cat(sprintf("  %-5s n=%-8g %s\n    -> b=%d vs b=1: %.2fx cheaper per column\n",
      k, nn, paste(sprintf("b%d=%.3fus", s$b, s$per_col_us), collapse = "  "),
      max(s$b), s$per_col_us[s$b == 1] / s$per_col_us[which.max(s$b)]))
}

## ---------------------------------------------------------------------------
## D. Krylov-scale extrapolation
## ---------------------------------------------------------------------------
cat("\n=== D. cumulative round-trip cost across a Krylov run ===\n")
for (it in c(1e2, 1e3, 1e4, 1e5))
  cat(sprintf("  %6g applies x %.0f ns = %8.4f s of pure call overhead\n", it, CALL * 1e9, it * CALL))

cat("\n=== E. smallest operator where the round trip is under 1% of an apply ===\n")
sub1 <- res[res$call_pct < 1, ]
if (nrow(sub1)) {
  s <- sub1[order(sub1$wrapped_us), ][1, ]
  cat(sprintf("  %s n=%g b=%d: apply %.3f us, round trip %.4f%%\n",
              s$kind, s$n, s$b, s$wrapped_us, s$call_pct))
}
cat(sprintf("  worst case in grid: %.4f%%\n", max(res$call_pct)))

out <- commandArgs(trailingOnly = TRUE)
write.csv(res, if (length(out)) out[1] else "S0.2_results.csv", row.names = FALSE)
cat("\nwrote results csv\n")
