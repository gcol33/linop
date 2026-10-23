## Applies, time and basis memory of the reference eigs() against the subspace
## size, for a default ncv that grows with n.
##
## Run from the repository root:
##   Rscript dev_notes/spikes/eigs-ncv-rule.R out.csv
##
## Fixtures: the 1-D Laplacian (a crowded end, gap ratio ~ 1/n^2) and the 2-D
## five-point Laplacian on an N x N grid (gap ratio ~ 1/n), largest end, k = 1
## and 5. Candidates: ncv 20, 40, 80, and min(n, max(2k + 1, 20, c log2(n) k))
## for c = 2 and 4. Memory is the basis alone, n * ncv doubles. Each run is
## capped at the default budget or CAP applies, whichever is smaller.
suppressMessages(devtools::load_all(".", quiet = TRUE, helpers = FALSE))
out_file <- commandArgs(TRUE)[1]
if (is.na(out_file)) out_file <- "eigs-ncv-rule.csv"
CAP <- 40000L

lap1 <- function(n) {
  f <- function(X) {
    Y <- 2 * X
    Y[-n, ] <- Y[-n, , drop = FALSE] - X[-1, , drop = FALSE]
    Y[-1, ] <- Y[-1, , drop = FALSE] - X[-n, , drop = FALSE]
    Y
  }
  linop(f, adjoint = f, dim = c(n, n), properties = c(hermitian = TRUE))
}
lap2 <- function(N) {
  n <- N * N
  f <- function(X) {
    Y <- 4 * X
    for (j in seq_len(ncol(X))) {
      G <- matrix(X[, j], N, N)
      H <- 4 * G
      H[-N, ] <- H[-N, ] - G[-1, ]
      H[-1, ] <- H[-1, ] - G[-N, ]
      H[, -N] <- H[, -N] - G[, -1]
      H[, -1] <- H[, -1] - G[, -N]
      Y[, j] <- as.vector(H)
    }
    Y
  }
  linop(f, adjoint = f, dim = c(n, n), properties = c(hermitian = TRUE))
}
rule <- function(n, k, c) as.integer(min(n, max(2 * k + 1, 20, ceiling(c * log2(n) * k))))

cases <- list()
for (n in c(500, 2000, 10000)) cases[[length(cases) + 1L]] <- list(name = "lap1d", n = n, A = lap1(n))
for (N in c(23, 45, 100, 317)) cases[[length(cases) + 1L]] <- list(name = "lap2d", n = N * N, A = lap2(N))

rows <- list()
for (cs in cases) for (k in c(1L, 5L)) {
  cands <- list(`20` = 20L, `40` = 40L, `80` = 80L,
                `rule c=2` = rule(cs$n, k, 2), `rule c=4` = rule(cs$n, k, 4))
  for (lab in names(cands)) {
    ncv <- max(cands[[lab]], k + 1L)
    t0 <- proc.time()[["elapsed"]]
    fit <- eigs(cs$A, k = k, which = "largest_algebraic", method = "lanczos",
                ncv = ncv, maxit = min(20L * cs$n, CAP))
    el <- proc.time()[["elapsed"]] - t0
    r <- data.frame(problem = cs$name, n = cs$n, k = k, candidate = lab, ncv = ncv,
                    applies = fit$iterations, nconv = fit$nconv, secs = round(el, 2),
                    basis_MB = round(cs$n * ncv * 8 / 2^20, 2))
    print(r, row.names = FALSE)
    rows[[length(rows) + 1L]] <- r
    write.csv(do.call(rbind, rows), out_file, row.names = FALSE)
  }
}
