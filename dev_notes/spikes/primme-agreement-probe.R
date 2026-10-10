## How many steps can a reference test assert between the Lanczos reference and
## PRIMME, and on what?
##
## The roster's ceilings are properties of the method: four steps for LSQR and
## LSMR, eight for MINRES, every step tried for BiCGSTAB. Here the two
## implementations are different methods by default (PRIMME's default preset
## switches between GD+k and JDQMR, with inner iterations), so only the answers
## can agree, never the iterates. With PRIMME's "Arnoldi" preset, no restart and
## a block of one, its basis is the Krylov space of the starting vector, which is
## the space the reference's first round builds, so the Ritz values after m
## applies are the same numbers in exact arithmetic. This measures how long they
## stay the same in floating point.
##
## Both start from the same unit vector. The reference runs one round of m
## applies (ncv > m, tol = 0, so nothing converges and nothing restarts); PRIMME
## runs maxMatvecs = m with a basis wider than m and its arithmetic floor at
## zero, so that it cannot stop early either. The quantity compared is the
## Rayleigh quotient of the top Ritz vector, measured on A in both cases, so the
## comparison is of the spaces and not of the two codes' bookkeeping.
##
## From the repository root: Rscript dev_notes/spikes/primme-agreement-probe.R

if (!file.exists("DESCRIPTION")) stop("run this from the repository root")
suppressMessages(devtools::load_all(".", quiet = TRUE, helpers = TRUE))
OUT_DIR <- file.path("dev_notes", "spikes", "results")
CSV <- file.path(OUT_DIR, "primme-agreement-probe.csv")

rq_top <- function(A, x) {
  x <- as_block(x)
  Ax <- linop_apply(A, x, "N")
  Re(sum(Conj(x) * Ax)) / Re(sum(Conj(x) * x))
}

ref_after <- function(A, m, v) {
  fit <- eigs(A, k = 1L, which = "largest_algebraic", method = "lanczos",
              maxit = m, ncv = m + 1L, tol = 0, v0 = v)
  list(value = fit$values[1], vec = fit$vectors[, 1], applies = fit$iterations)
}

primme_after <- function(A, m, v) {
  out <- linop.primme::primme_eigsh(
    function(X) linop_apply(A, X, "N"), A$dim[1], 1L, which = "LA",
    tol = .Machine$double.eps, maxit = m, ncv = 80L, min_restart = 1L, block = 1L,
    v0 = v, complex = A$dtype == "complex", floor = 0, preset = "PRIMME_Arnoldi")
  list(value = rq_top(A, out$vectors[, 1]), vec = out$vectors[, 1],
       applies = out$stats$numMatvecs)
}

fixtures <- list(
  list(name = "spd kappa 1e2", make = function(seed)
    as_spd_linop(spd_prescribed(60, exp(seq(log(100), 0, length.out = 60)), seed = seed))),
  list(name = "spd kappa 1e6", make = function(seed)
    as_spd_linop(spd_prescribed(60, exp(seq(log(1e6), 0, length.out = 60)), seed = seed))),
  list(name = "hpd kappa 1e4", make = function(seed) {
    M <- hpd_prescribed(60, exp(seq(log(1e4), 0, length.out = 60)), seed = seed)
    linop(M, properties = c(hermitian = TRUE))
  }),
  list(name = "laplacian 400", make = function(seed) laplacian_1d(400))
)

steps <- c(1L, 2L, 4L, 6L, 8L, 12L, 16L, 24L, 32L, 40L)
rows <- list()
for (fx in fixtures) for (seed in 1:5) {
  A <- fx$make(seed)
  v <- eigen_start_vector(A$dim[1], A$dtype, NULL, 100L + seed)
  for (m in steps) {
    r <- ref_after(A, m, v)
    p <- primme_after(A, m, v)
    cosang <- Mod(sum(Conj(r$vec) * p$vec)) /
      sqrt(Re(sum(Conj(r$vec) * r$vec)) * Re(sum(Conj(p$vec) * p$vec)))
    rows[[length(rows) + 1L]] <- data.frame(
      fixture = fx$name, seed = seed, m = m, ref_applies = r$applies,
      primme_applies = p$applies,
      value_rel = abs(r$value - p$value) / abs(r$value),
      vector_sin = sqrt(max(0, 1 - cosang^2)))
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, CSV, row.names = FALSE)

agg <- aggregate(cbind(value_rel, vector_sin) ~ fixture + m, res, max)
agg <- agg[order(agg$fixture, agg$m), ]
print(format(agg, digits = 3), row.names = FALSE)
print(unique(res[, c("m", "ref_applies", "primme_applies")]), row.names = FALSE)
