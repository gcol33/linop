## What PRIMME's block size buys and costs, called on the engine directly.
##
## A Krylov-type method grown from one vector reaches one vector of each
## eigenspace in exact arithmetic, and only rounding supplies the others. The
## reference recovers the missing copies by restarting from a random vector once
## the locked set has deflated its start; PRIMME's remedy is a block of starting
## directions. This measures, on 2-, 3-, 5- and 8-fold eigenvalues at the top of
## prescribed spectra of order 30, 60 and 200 over ten seeds, how many copies a
## block of one and a block of k return, against the reference, and on the
## n = 1e5 well of primme-restart-probe.R what each
## block size costs in applies and seconds. Last, the smallest basis (ncv)
## PRIMME accepts with a block of k, for each target, since below it PRIMME's
## own parameter check fails.
##
## From the repository root: Rscript dev_notes/spikes/primme-block-probe.R

if (!file.exists("DESCRIPTION")) stop("run this from the repository root")
OUT_DIR <- file.path("dev_notes", "spikes", "results")
LOG <- file.path(OUT_DIR, "primme-block-probe.txt")
sink(LOG, split = TRUE)

eng <- linop.primme::primme_eigsh
blocks <- c(1L, 2L, 4L)
suppressMessages(devtools::load_all(".", quiet = TRUE, helpers = FALSE))

## Copies of a k-fold eigenvalue at the top of a prescribed spectrum, asked for
## with k = its multiplicity: PRIMME's GD+k (the engine's default preset) with
## a block of one and a block of k, and the reference through eigs().
cat("copies of a k-fold top eigenvalue returned, over ten seeds
")
for (m in c(30L, 60L, 200L)) for (mult in c(2L, 3L, 5L, 8L)) {
  got <- matrix(0L, 3L, 10L, dimnames = list(c("primme block 1", "primme block k",
                                                 "lanczos"), NULL))
  for (seed in 1:10) {
    set.seed(seed)
    Q <- qr.Q(qr(matrix(rnorm(m * m), m)))
    lam <- c(rep(10, mult), seq(2, 0.5, length.out = m - mult))
    M <- Q %*% diag(lam) %*% t(Q)
    v0 <- rnorm(m)
    ap <- function(X) M %*% X
    got[1L, seed] <- sum(abs(eng(ap, m, mult, which = "LA", tol = 1e-10, block = 1L,
                                 v0 = v0)$values - 10) < 1e-8)
    got[2L, seed] <- sum(abs(eng(ap, m, mult, which = "LA", tol = 1e-10, block = mult,
                                 v0 = v0)$values - 10) < 1e-8)
    fit <- eigs(linop(M, properties = c(hermitian = TRUE)), k = mult,
                which = "largest_algebraic", method = "lanczos", v0 = v0)
    got[3L, seed] <- sum(abs(fit$values - 10) < 1e-8)
  }
  cat(sprintf("  n = %3d, multiplicity %d
", m, mult))
  for (i in 1:3) cat(sprintf("    %-15s %s
", rownames(got)[i], paste(got[i, ], collapse = " ")))
}

n <- 1e5
v <- numeric(n); mid <- n %/% 2; v[(mid - 14):(mid + 15)] <- 1
well <- function(X) {
  X <- as.matrix(X)
  Y <- 2 * X
  Y[-n, ] <- Y[-n, , drop = FALSE] - X[-1L, , drop = FALSE]
  Y[-1L, ] <- Y[-1L, , drop = FALSE] - X[-n, , drop = FALSE]
  Y - v * X
}
cat("\nwell, n = 1e5, k = 4 smallest algebraic, tol 1e-10\n")
for (b in blocks) {
  t <- system.time(f <- eng(well, n, 4L, which = "SA", tol = 1e-10, block = b,
                            maxit = 20000L))[["elapsed"]]
  cat(sprintf("  block %d: %4d applies, %.2f s, %d converged, basis %d\n", b,
              as.integer(f$stats$numMatvecs), t, f$nconv, f$params$maxBasisSize))
}

cat("
smallest ncv PRIMME runs with, block = k (and every larger one up to 4k + 4 runs)
")
nl <- 300
lap <- function(X) 2 * X - rbind(X[-1, , drop = FALSE], 0) - rbind(0, X[-nl, , drop = FALSE])
set.seed(1)
Am <- matrix(rnorm(300 * 200), 300)
rect <- function(X, trans) if (trans) crossprod(Am, X) else Am %*% X
runs <- function(f) tryCatch({ f(); TRUE }, error = function(e) FALSE)
smallest <- function(k, call) {
  grid <- (k + 1L):(4L * k + 4L)
  ok <- vapply(grid, function(ncv) runs(function() call(ncv)), logical(1))
  first <- which(ok)[1L]
  c(ncv = grid[first], monotone = all(ok[first:length(ok)]))
}
for (w in c("LA", "SA", "LM", "SM")) {
  r <- vapply(1:8, function(k) smallest(k, function(ncv) {
    f <- eng(lap, nl, k, which = w, ncv = ncv, block = k, maxit = 50L)
    if (!f$info %in% c(0L, -3L)) stop(f$info)
  }), numeric(2))
  cat(sprintf("  eigs %s  k = 1..8: %s%s
", w, paste(r["ncv", ], collapse = " "),
              if (all(r["monotone", ] == 1)) "" else "  (not monotone)"))
}
for (w in c("LM", "SM")) {
  r <- vapply(1:6, function(k) smallest(k, function(ncv) {
    f <- linop.primme::primme_svds(rect, 300, 200, k, which = w, ncv = ncv,
                                   block = k, maxit = 50L)
    if (!f$info %in% c(0L, -103L, -203L)) stop(f$info)
  }), numeric(2))
  cat(sprintf("  svds %s  k = 1..6: %s%s
", w, paste(r["ncv", ], collapse = " "),
              if (all(r["monotone", ] == 1)) "" else "  (not monotone)"))
}
sink()
