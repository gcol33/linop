## BiCGSTAB: does it solve, on all three sides, and where do the three inner
## products that become denominators actually go on a healthy solve.
##
##   Rscript dev_notes/spikes/bicgstab_probe.R

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

bicgstab <- linop:::bicgstab_solve
gmres <- linop:::gmres_solve

cat("=== 1. recovery against closed-form truth ===\n")

nn <- 40
cases <- list(
  list(name = "convdiff mu=.3", A = convdiff_1d(nn, 0.3),
       truth = function(b) convdiff_1d_solve(nn, 0.3, b), n = nn),
  list(name = "laplacian", A = laplacian_1d(nn),
       truth = function(b) {
         V <- laplacian_1d_eigenvectors(nn)
         V %*% (crossprod(V, b) / laplacian_1d_eigenvalues(nn))
       }, n = nn),
  list(name = "kms rho=.5", A = as_spd_linop(kms_matrix(nn, 0.5)),
       truth = function(b) kms_inverse(nn, 0.5) %*% b, n = nn),
  list(name = "shifted lap", A = shifted_laplacian_1d(nn, 0.9),
       truth = function(b) shifted_laplacian_solve(nn, 0.9, b), n = nn))

for (cs in cases) {
  set.seed(1)
  b <- matrix(stats::rnorm(cs$n), cs$n, 1)
  fit <- bicgstab(cs$A, b, tol = 1e-11, maxit = 2000L)
  tr <- cs$truth(b)
  reason <- fit$certificate$checks$detail[fit$certificate$checks$check == "convergence"]
  cat(sprintf("  %-16s %4d it  conv %-5s  rel err %.2e  residual %-6s  breakdown reported %s\n",
              cs$name, fit$iterations, fit$converged,
              max(Mod(fit$x - tr)) / max(Mod(tr)),
              linop:::cert_status(fit$certificate, "residual"),
              grepl("broke down", reason)))
}

cat("\n=== 2. complex ===\n")
set.seed(2)
Mz <- zmat(30, 30) + diag(8, 30)
Az <- linop(Mz)
bz <- matrix(complex(real = stats::rnorm(30), imaginary = stats::rnorm(30)), 30, 1)
fz <- bicgstab(Az, bz, tol = 1e-11, maxit = 2000L)
cat(sprintf("  %4d it  conv %s  rel err %.2e\n", fz$iterations, fz$converged,
            max(Mod(fz$x - solve(Mz, bz))) / max(Mod(solve(Mz, bz)))))

cat("\n=== 3. the three sides, and whether they differ ===\n")
## The diagonal of convdiff_1d is constant, so a Jacobi preconditioner on it is a
## scalar and all three sides coincide trivially. A preconditioner has to have a
## non-constant diagonal before the question can even be asked.
n <- 40
A <- convdiff_1d(n, 0.3)
d <- exp(seq(-1.5, 1.5, length.out = n))
set.seed(3)
b <- matrix(stats::rnorm(n), n, 1)
truth <- convdiff_1d_solve(n, 0.3, b)
mk <- function(sd) preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                                  positive_definite = TRUE, side = sd)
xs <- list()
for (sd in c("left", "right", "split")) {
  fit <- bicgstab(A, b, tol = 1e-11, preconditioner = mk(sd), maxit = 2000L)
  cat(sprintf("  %-6s %4d it  conv %-5s  rel err %.2e\n", sd, fit$iterations,
              fit$converged, max(Mod(fit$x - truth)) / max(Mod(truth))))
}
## at a fixed small budget, do the three produce different iterates
for (sd in c("left", "right", "split")) {
  xs[[sd]] <- bicgstab(A, b, tol = 0, preconditioner = mk(sd), maxit = 3L)$x
}
cat(sprintf("  at 3 steps: left vs right %.3e, left vs split %.3e, right vs split %.3e\n",
            max(Mod(xs$left - xs$right)), max(Mod(xs$left - xs$split)),
            max(Mod(xs$right - xs$split))))
cat(sprintf("  and a scalar preconditioner instead: left vs right %.3e\n", {
  ds <- rep(2, n)
  mk2 <- function(sd) preconditioner(function(R) R / ds, fixed = TRUE, hermitian = TRUE,
                                     positive_definite = TRUE, side = sd)
  l <- bicgstab(A, b, tol = 0, preconditioner = mk2("left"), maxit = 3L)$x
  r <- bicgstab(A, b, tol = 0, preconditioner = mk2("right"), maxit = 3L)$x
  max(Mod(l - r))
}))

cat("\n=== 4. where the three denominators go on a healthy solve ===\n")
## instrumented single-column BiCGSTAB, unpreconditioned, recording the relative
## size of each inner product that becomes a denominator
probe <- function(M, b, steps) {
  n <- ncol(M)
  x <- rep(0, n)
  if (is.complex(M) || is.complex(b)) storage.mode(x) <- "complex"
  r <- as.vector(b); rhat <- r; p <- r
  nrm <- function(v) sqrt(sum(Mod(v)^2))
  rho <- sum(Conj(rhat) * r)
  worst <- c(rho = Inf, denom = Inf, ts = Inf)
  for (j in seq_len(steps)) {
    v <- as.vector(M %*% p)
    den <- sum(Conj(rhat) * v)
    sc <- nrm(rhat) * nrm(v)
    if (sc > 0) worst["denom"] <- min(worst["denom"], Mod(den) / sc)
    if (Mod(den) == 0) break
    alpha <- rho / den
    s <- r - alpha * v
    if (nrm(s) <= 1e-14 * nrm(b)) { x <- x + alpha * p; break }
    tt <- as.vector(M %*% s)
    ts <- sum(Conj(tt) * s)
    sc <- nrm(tt) * nrm(s)
    if (sc > 0) worst["ts"] <- min(worst["ts"], Mod(ts) / sc)
    if (Mod(ts) == 0) break
    omega <- ts / sum(Conj(tt) * tt)
    x <- x + alpha * p + omega * s
    r <- s - omega * tt
    rho_new <- sum(Conj(rhat) * r)
    sc <- nrm(rhat) * nrm(r)
    if (sc > 0) worst["rho"] <- min(worst["rho"], Mod(rho_new) / sc)
    if (nrm(r) <= 1e-13 * nrm(b)) break
    beta <- (rho_new / rho) * (alpha / omega)
    p <- r + beta * (p - omega * v)
    rho <- rho_new
  }
  worst
}

fixtures <- list(
  "convdiff mu=.3" = as.matrix(convdiff_1d(40, 0.3)),
  "convdiff mu=.7" = as.matrix(convdiff_1d(40, 0.7)),
  "laplacian" = as.matrix(laplacian_1d(40)),
  "kms rho=.7" = kms_matrix(40, 0.7),
  "shifted lap" = as.matrix(shifted_laplacian_1d(40, 0.9)))
set.seed(9)
Mzz <- zmat(30, 30) + diag(8, 30)
fixtures[["complex dense"]] <- Mzz

cat("  smallest relative value each denominator reached, over 12 seeds\n")
cat(sprintf("  %-16s %10s %10s %10s\n", "fixture", "rho", "<rhat,v>", "<t,s>"))
for (nm in names(fixtures)) {
  M <- fixtures[[nm]]
  w <- c(rho = Inf, denom = Inf, ts = Inf)
  for (seed in 1:12) {
    set.seed(100L + seed)
    bb <- if (is.complex(M))
      matrix(complex(real = stats::rnorm(nrow(M)), imaginary = stats::rnorm(nrow(M))), nrow(M), 1)
    else matrix(stats::rnorm(nrow(M)), nrow(M), 1)
    w <- pmin(w, probe(M, bb, 300L))
  }
  cat(sprintf("  %-16s %10.2e %10.2e %10.2e\n", nm, w["rho"], w["denom"], w["ts"]))
}

cat("\n=== 5. against gmres, on the same fixtures ===\n")
for (cs in cases) {
  set.seed(1)
  b <- matrix(stats::rnorm(cs$n), cs$n, 1)
  bb <- bicgstab(cs$A, b, tol = 1e-11, maxit = 2000L)
  gg <- gmres(cs$A, b, tol = 1e-11, maxit = 2000L, restart = 30L)
  tr <- cs$truth(b)
  cat(sprintf("  %-16s bicgstab %4d it (%4d applies) err %.2e | gmres %4d it err %.2e\n",
              cs$name, bb$iterations, 2 * bb$iterations,
              max(Mod(bb$x - tr)) / max(Mod(tr)),
              gg$iterations, max(Mod(gg$x - tr)) / max(Mod(tr))))
}

cat("\n=== 6. non-monotone residual ===\n")
n <- 60
A <- convdiff_1d(n, 0.6)
set.seed(5)
b <- matrix(stats::rnorm(n), n, 1)
M <- as.matrix(A)
tr <- function(k) {
  x <- bicgstab(A, b, tol = 0, maxit = k)$x
  sqrt(sum((b - M %*% x)^2))
}
v <- vapply(seq_len(25), tr, numeric(1))
cat(sprintf("  rises: %d of 24 steps, worst +%.1f%%\n", sum(diff(v) > 0),
            100 * max(c(0, diff(v) / utils::head(v, -1)))))
cat(sprintf("  %s\n", paste(sprintf("%.2e", v[1:12]), collapse = " ")))

cat("\n=== 7. lockstep and vector in, vector out ===\n")
n <- 40
A <- convdiff_1d(n, 0.3)
set.seed(6)
BB <- matrix(stats::rnorm(n * 4), n, 4)
blk <- bicgstab(A, BB, tol = 1e-11, maxit = 2000L)
ok <- TRUE
for (j in 1:4) {
  one <- bicgstab(A, BB[, j, drop = FALSE], tol = 1e-11, maxit = 2000L)
  ok <- ok && identical(blk$x[, j], one$x[, 1L])
}
cat(sprintf("  lockstep bitwise: %s\n", ok))
vv <- bicgstab(A, as.numeric(BB[, 1]), tol = 1e-11, maxit = 2000L)
cat(sprintf("  vector in vector out: %s\n", is.null(dim(vv$x))))

cat("\n=== 9. a guaranteed breakdown ===\n")
## For a real skew-symmetric A, <A z, z> = 0 for every real z, so <rhat0, v> is
## exactly zero at the first step of every round however the shadow vector is
## re-seeded. BiCGSTAB cannot solve this system, and the question is what comes
## back rather than whether it converges.
n <- 20
Sk <- diag(0, n)
for (i in seq_len(n - 1L)) { Sk[i, i + 1L] <- 1; Sk[i + 1L, i] <- -1 }
Ask <- linop(Sk)
set.seed(8)
bsk <- matrix(stats::rnorm(n), n, 1)
fsk <- bicgstab(Ask, bsk, tol = 1e-10, maxit = 200L)
cat(sprintf("  converged %s, iterations %d, finite x %s, any NaN %s\n",
            fsk$converged, fsk$iterations, all(is.finite(fsk$x)), any(is.na(fsk$x))))
cat(sprintf("  convergence row: %s\n",
            fsk$certificate$checks$detail[fsk$certificate$checks$check == "convergence"]))
cat(sprintf("  overall %s, residual %s\n", fsk$certificate$overall,
            linop:::cert_status(fsk$certificate, "residual")))
cat(sprintf("  and gmres on the same system: err %.2e\n", {
  g <- gmres(Ask, bsk, tol = 1e-10, maxit = 200L, restart = n)
  max(Mod(g$x - solve(Sk, bsk))) / max(Mod(solve(Sk, bsk)))
}))

cat("\n=== 8. no adjoint needed ===\n")
S <- diag(0, 25); for (i in seq_len(24)) S[i, i + 1L] <- 1
Mno <- diag(4, 25) + S
Ano <- linop(function(X) Mno %*% X, dim = c(25L, 25L))
set.seed(7)
bno <- matrix(stats::rnorm(25), 25, 1)
fno <- bicgstab(Ano, bno, tol = 1e-11, maxit = 2000L)
cat(sprintf("  conv %s, rel err %.2e\n", fno$converged,
            max(Mod(fno$x - solve(Mno, bno))) / max(Mod(solve(Mno, bno)))))
