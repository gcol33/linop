## Gate 2 asks that MINRES agree with a trusted reference on ill-conditioned
## problems, including breakdown and near-breakdown. What has to be measured
## before any of that can be asserted:
##
##   1. where the dense reference stops being ground truth. It orthogonalises a
##      basis built by repeated application of A, so it drifts for the same
##      reason the recurrence does, and a reference whose own error is unmeasured
##      is not a reference.
##   2. how far this implementation's iterate sits from that reference as the
##      condition number grows, which is what sets a tolerance rather than
##      choosing one.
##   3. what happens at an exact breakdown, beta_{j+1} = 0, and 4. through the
##      near-breakdown window where beta_{j+1} is small and nonzero.
##   5. the same against the published recurrence rather than the definition, and
##      the two references against each other, which is where the ceiling turns
##      out to belong to the method rather than to any one program.
##   6. the complex case, which has no external reference at all.
##
## The external validation of the transcription against SciPy 1.17.1 is in
## minres_export.R with minres_scipy.py.
##
##   Rscript dev_notes/spikes/minres_reference.R

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

minres <- linop:::minres_solve
cert_status <- linop:::cert_status

cdot <- function(x, y) sum(Conj(x) * y)
nrm <- function(x) sqrt(Re(cdot(x, x)))

## The exact minimiser of ||b - A x|| over K_m(A, b), computed from the
## definition. Two independent pieces: an orthonormal basis for the space, and a
## least-squares solve over it. The second is a truncated SVD rather than a QR so
## that the reference degrades gracefully when the first has lost rank, and so
## that the same code serves real and complex.
krylov_basis <- function(M, b, m, passes = 2L, drop_tol = 1e-12) {
  n <- length(b)
  cx <- is.complex(M) || is.complex(b)
  K <- matrix(if (cx) complex(1) else 0, n, m)
  v <- b / nrm(b)
  used <- 0L
  for (j in seq_len(m)) {
    K[, j] <- v
    used <- j
    v <- as.vector(M %*% v)
    n0 <- nrm(v)
    for (p in seq_len(passes)) {
      for (i in seq_len(j)) v <- v - cdot(K[, i], v) * K[, i]
    }
    nv <- nrm(v)
    ## Relative, so scaling the operator does not move where the space is
    ## declared exhausted.
    if (n0 <= 0 || nv <= drop_tol * n0) break
    v <- v / nv
  }
  K[, seq_len(used), drop = FALSE]
}

pinv_solve <- function(G, b) {
  s <- svd(G)
  keep <- s$d > max(dim(G)) * .Machine$double.eps * max(s$d)
  y <- rep(if (is.complex(G) || is.complex(b)) complex(1) else 0, ncol(G))
  if (any(keep)) {
    y[keep] <- as.vector(s$v[, keep, drop = FALSE] %*%
                         (as.vector(crossprod(Conj(s$u[, keep, drop = FALSE]), b)) /
                          s$d[keep]))
  }
  y
}

krylov_argmin <- function(M, b, m, passes = 2L) {
  K <- krylov_basis(M, b, m, passes)
  as.vector(K %*% pinv_solve(M %*% K, b))
}

## The indefinite prescribed spectrum used throughout: half below zero and half
## above, magnitudes from 1 down to 1/kappa, so kappa_2 is exactly kappa and the
## operator is one no definite method may be given.
indef_spectrum <- function(n, kappa) {
  mag <- exp(seq(log(1), log(1 / kappa), length.out = n))
  sign <- rep(c(-1, 1), length.out = n)
  sign * mag
}

indef_operator <- function(n, kappa, seed) {
  spd_prescribed(n, indef_spectrum(n, kappa), seed = seed)
}

rhs <- function(n, seed) {
  set.seed(seed)
  matrix(stats::rnorm(n), n, 1L)
}

## ------------------------------------------------------------------ part 1 --
## Is the reference ground truth? Three independent readings, none of which
## depends on this package's recurrence.

cat("=== 1. the dense reference measured against itself ===\n")
cat("    orth   ||K^H K - I||_max, whether the basis is a basis\n")
cat("    normeq ||(AK)^H r|| / (||AK|| ||r||), whether y minimises\n")
cat("    1v2    the same reference with one orthogonalisation pass instead of\n")
cat("           two, which is the drift in the space rather than in the solve\n\n")
cat(sprintf("%-10s %5s  %10s  %10s  %10s\n", "kappa", "steps", "orth", "normeq", "1v2"))

n <- 40
for (kappa in c(1e2, 1e4, 1e6, 1e8, 1e10)) {
  M <- indef_operator(n, kappa, seed = 1L)
  b <- rhs(n, 700L)
  for (m in c(2L, 4L, 8L, 16L)) {
    K <- krylov_basis(M, b, m)
    G <- M %*% K
    y <- pinv_solve(G, b)
    r <- b - G %*% y
    orth <- max(Mod(crossprod(Conj(K), K) - diag(1, ncol(K))))
    gn <- max(svd(G, nu = 0L, nv = 0L)$d)
    normeq <- nrm(as.vector(crossprod(Conj(G), r))) / (gn * max(nrm(r), .Machine$double.xmin))
    x2 <- krylov_argmin(M, b, m, passes = 2L)
    x1 <- krylov_argmin(M, b, m, passes = 1L)
    cat(sprintf("%-10.0e %5d  %10.2e  %10.2e  %10.2e\n", kappa, m, orth, normeq,
                max(Mod(x1 - x2)) / max(Mod(x2))))
  }
}

## ------------------------------------------------------------------ part 2 --
## This implementation against that reference, as the conditioning grows. tol = 0
## removes every stopping test so both run the same fixed number of steps and the
## gap is the iteration rather than the budget.

cat("\n=== 2. this implementation against the reference, worst over 5 seeds ===\n")
cat("    relative gap in the iterate at step m, tol = 0 so the step count is\n")
cat("    the same on both sides\n\n")
cat(sprintf("%-10s %5s  %10s  %10s\n", "kappa", "steps", "worst gap", "at seed"))

for (kappa in c(1e2, 1e4, 1e6, 1e8, 1e10)) {
  for (m in c(1L, 2L, 4L, 8L, 16L)) {
    worst <- 0
    at <- NA_integer_
    for (seed in 1:5) {
      M <- indef_operator(n, kappa, seed = seed)
      A <- linop(M, properties = c(hermitian = TRUE))
      b <- rhs(n, 700L + seed)
      mine <- minres(A, b, tol = 0, maxit = m)
      ref <- krylov_argmin(M, b, m)
      g <- max(Mod(as.vector(mine$x) - ref)) / max(Mod(ref))
      if (g > worst) { worst <- g; at <- seed }
    }
    cat(sprintf("%-10.0e %5d  %10.2e  %10d\n", kappa, m, worst, at))
  }
}

cat("\n=== 2b. run to convergence instead of to a step count ===\n")
cat("    the other half of the same fact: what the two agree on is the\n")
cat("    solution, whatever they did on the way to it\n\n")
cat(sprintf("%-10s %5s  %8s  %12s  %12s\n", "kappa", "seed", "iters", "vs truth", "kappa*eps"))

for (kappa in c(1e2, 1e4, 1e6, 1e8)) {
  for (seed in 1:3) {
    M <- indef_operator(n, kappa, seed = seed)
    A <- linop(M, properties = c(hermitian = TRUE))
    x_true <- rhs(n, 800L + seed)
    b <- M %*% x_true
    fit <- minres(A, b, tol = 1e-13, maxit = 40L * n)
    cat(sprintf("%-10.0e %5d  %8d  %12.2e  %12.2e\n", kappa, seed, fit$iterations,
                max(Mod(as.vector(fit$x) - as.vector(x_true))) / max(Mod(x_true)),
                kappa * .Machine$double.eps))
  }
}

## ------------------------------------------------------------------ part 3 --
## Exact breakdown. b inside a d-dimensional invariant subspace makes
## beta_{d+1} = 0 to rounding: the Krylov space is exhausted and the iterate
## there is the exact solution. The shifted Laplacian has closed-form
## eigenvectors, so the subspace is constructed rather than found.

cat("\n=== 3. exact breakdown, b inside a d-dimensional invariant subspace ===\n")
cat("    what the recurrence should do is terminate at d with the answer, and\n")
cat("    what the reference should do is find the space exhausted at the same\n")
cat("    step\n\n")
cat(sprintf("%-4s %8s  %8s  %12s  %12s  %10s\n",
            "d", "iters", "ref dim", "vs truth", "ref vs truth", "cert"))

nl <- 60
sigma <- 1.5
V <- laplacian_1d_eigenvectors(nl)
Ash <- shifted_laplacian_1d(nl, sigma)
Msh <- as.matrix(Ash)

for (d in c(1L, 2L, 3L, 5L, 8L)) {
  idx <- seq_len(d) * 5L
  set.seed(900L + d)
  b <- matrix(V[, idx, drop = FALSE] %*% stats::rnorm(d), nl, 1L)
  truth <- shifted_laplacian_solve(nl, sigma, b)
  fit <- minres(Ash, b, tol = 1e-12, maxit = 10L * nl)
  K <- krylov_basis(Msh, b, nl)
  ref <- krylov_argmin(Msh, b, d)
  cat(sprintf("%-4d %8d  %8d  %12.2e  %12.2e  %10s\n", d, fit$iterations, ncol(K),
              max(Mod(as.vector(fit$x) - as.vector(truth))) / max(Mod(truth)),
              max(Mod(ref - as.vector(truth))) / max(Mod(truth)),
              cert_status(fit$certificate, "residual")))
}

## The off-diagonal sequence itself, by a dense Lanczos with no reorthogonalisation
## at all, which is the recurrence MINRES runs. In exact arithmetic beta_{d+1} = 0
## and every later one is undefined; what happens instead is the question.
cat("\n=== 3b. the beta sequence on the same fixtures ===\n")
cat("    beta_{j+1} / ||A v_j||, so the row is scale free. Exact arithmetic\n")
cat("    puts a zero at j = d and nothing after it\n\n")

beta_sequence <- function(M, b, steps) {
  out <- rep(NA_real_, steps)
  v <- b / nrm(b)
  v_old <- rep(0, length(b))
  beta <- 0
  for (j in seq_len(steps)) {
    p <- as.vector(M %*% v)
    al <- Re(cdot(v, p))
    w <- p - al * v - beta * v_old
    bn <- nrm(w)
    out[j] <- bn / nrm(p)
    if (bn <= 0) break
    v_old <- v
    v <- w / bn
    beta <- bn
  }
  out
}

for (d in c(1L, 2L, 3L, 5L, 8L)) {
  idx <- seq_len(d) * 5L
  set.seed(900L + d)
  b <- matrix(V[, idx, drop = FALSE] %*% stats::rnorm(d), nl, 1L)
  bs <- beta_sequence(Msh, b, 12L)
  cat(sprintf("d = %d:  %s\n", d,
              paste(sprintf("%.1e", bs[seq_len(min(12L, d + 4L))]), collapse = " ")))
}

## If the sequence does not stop at d, the question is whether continuing costs
## anything. tol = 0 forces the extra steps to be taken rather than skipped.
cat("\n=== 3c. running past the breakdown, tol = 0 ===\n")
cat("    forward error at maxit = d, d+1, ... d+8\n\n")

for (d in c(1L, 3L, 5L, 8L)) {
  idx <- seq_len(d) * 5L
  set.seed(900L + d)
  b <- matrix(V[, idx, drop = FALSE] %*% stats::rnorm(d), nl, 1L)
  truth <- shifted_laplacian_solve(nl, sigma, b)
  errs <- vapply(d + 0:8, function(m) {
    x <- minres(Ash, b, tol = 0, maxit = m)$x
    max(Mod(as.vector(x) - as.vector(truth))) / max(Mod(truth))
  }, numeric(1))
  cat(sprintf("d = %d:  %s\n", d, paste(sprintf("%.1e", errs), collapse = " ")))
}

## ------------------------------------------------------------------ part 4 --
## Near-breakdown. b = v_i + delta v_j is inside a one-dimensional invariant
## subspace to relative accuracy delta, so beta_2 is O(delta): the window between
## a breakdown the guard treats as exhaustion and an ordinary small number. The
## guard is ||w|| <= sqrt(eps) ||A v||, so 1e-8 is where the behaviour changes
## and the sweep has to cross it.

cat("\n=== 4. near-breakdown, b = v_i + delta v_j ===\n")
cat("    beta2 is the first off-diagonal the recurrence forms, measured by a\n")
cat("    dense Lanczos step; sqrt(eps) = 1.49e-08 is where the exhaustion guard\n")
cat("    sits\n\n")
cat(sprintf("%-10s %10s  %8s  %12s  %12s  %10s\n",
            "delta", "beta2/||Av||", "iters", "vs truth", "residual", "cert"))

beta2_of <- function(M, b) {
  v <- b / nrm(b)
  p <- as.vector(M %*% v)
  al <- Re(cdot(v, p))
  w <- p - al * v
  nrm(w) / nrm(p)
}

for (delta in c(1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16, 0)) {
  b <- matrix(V[, 7L] + delta * V[, 23L], nl, 1L)
  truth <- shifted_laplacian_solve(nl, sigma, b)
  fit <- minres(Ash, b, tol = 1e-12, maxit = 10L * nl)
  cat(sprintf("%-10.0e %10.2e  %8d  %12.2e  %12.2e  %10s\n", delta,
              beta2_of(Msh, b), fit$iterations,
              max(Mod(as.vector(fit$x) - as.vector(truth))) / max(Mod(truth)),
              max(fit$certificate$values$residual),
              cert_status(fit$certificate, "residual")))
}

cat("\n=== 4b. the same window against the reference, step by step ===\n")
cat("    at a fixed step count rather than run to convergence, so what is\n")
cat("    compared is the iterate and not the answer\n\n")
cat(sprintf("%-10s %5s  %12s\n", "delta", "steps", "gap"))

for (delta in c(1e-4, 1e-8, 1e-12, 1e-16)) {
  b <- matrix(V[, 7L] + delta * V[, 23L], nl, 1L)
  for (m in c(1L, 2L, 3L)) {
    mine <- minres(Ash, b, tol = 0, maxit = m)
    ref <- krylov_argmin(Msh, b, m)
    cat(sprintf("%-10.0e %5d  %12.2e\n", delta, m,
                max(Mod(as.vector(mine$x) - ref)) / max(Mod(ref))))
  }
}

## ------------------------------------------------------------------ part 5 --
## The published recurrence, transcribed from the Paige and Saunders form and
## checked against SciPy by minres_export.R. What is measured here is the gap
## between this implementation and that transcription, over the same ladder, so
## the test suite's tolerances come from a reference it can carry itself rather
## than from one that needs Python.

reference_minres <- function(M, b, steps) {
  n <- length(b)
  eps <- .Machine$double.eps
  cx <- is.complex(M) || is.complex(b)
  zero <- if (cx) complex(n) else numeric(n)
  x <- zero
  r1 <- as.vector(b)
  y <- r1
  beta1 <- sqrt(Re(cdot(r1, y)))
  if (beta1 <= 0) return(x)
  oldb <- 0; beta <- beta1; dbar <- 0; epsln <- 0
  phibar <- beta1
  cs <- -1; sn <- 0
  w <- zero; w2 <- zero
  r2 <- r1
  for (itn in seq_len(steps)) {
    if (beta <= 0) break
    v <- y / beta
    y <- as.vector(M %*% v)
    if (itn >= 2) y <- y - (beta / oldb) * r1
    ## alpha is real for a hermitian operator; taking the real part is the
    ## statement of that rather than a repair of it.
    alfa <- Re(cdot(v, y))
    y <- y - (alfa / beta) * r2
    r1 <- r2; r2 <- y
    oldb <- beta
    beta <- sqrt(max(Re(cdot(r2, y)), 0))
    oldeps <- epsln
    delta <- cs * dbar + sn * alfa
    gbar <- sn * dbar - cs * alfa
    epsln <- sn * beta
    dbar <- -cs * beta
    gamma <- max(sqrt(gbar^2 + beta^2), eps)
    cs <- gbar / gamma
    sn <- beta / gamma
    phi <- cs * phibar
    phibar <- sn * phibar
    w1 <- w2; w2 <- w
    w <- (v - oldeps * w1 - delta * w2) / gamma
    x <- x + phi * w
  }
  x
}

cat("\n=== 5. against the published recurrence, worst over 5 seeds ===\n")
cat("    def  the dense minimiser, the definition\n")
cat("    pub  the transcribed Paige and Saunders recurrence\n\n")
cat(sprintf("%-10s %5s  %11s  %11s\n", "kappa", "steps", "def", "pub"))

for (kappa in c(1e2, 1e4, 1e6, 1e8, 1e10)) {
  for (m in c(1L, 2L, 4L, 6L, 8L, 12L)) {
    wd <- 0; wp <- 0
    for (seed in 1:5) {
      M <- indef_operator(n, kappa, seed = seed)
      A <- linop(M, properties = c(hermitian = TRUE))
      b <- rhs(n, 700L + seed)
      mine <- as.vector(minres(A, b, tol = 0, maxit = m)$x)
      rd <- krylov_argmin(M, b, m)
      rp <- reference_minres(M, as.vector(b), m)
      wd <- max(wd, max(Mod(mine - rd)) / max(Mod(rd)))
      wp <- max(wp, max(Mod(mine - rp)) / max(Mod(rp)))
    }
    cat(sprintf("%-10.0e %5d  %11.2e  %11.2e\n", kappa, m, wd, wp))
  }
}

cat("\n=== 5b. the two references against each other ===\n")
cat("    a variational characterisation and a short recurrence, sharing no\n")
cat("    code, so the step where they part is the step past which neither is\n")
cat("    ground truth for the other\n\n")
cat(sprintf("%-10s %5s  %11s\n", "kappa", "steps", "def vs pub"))
for (kappa in c(1e2, 1e6, 1e10)) {
  for (m in c(2L, 4L, 8L, 12L, 16L)) {
    w <- 0
    for (seed in 1:5) {
      M <- indef_operator(n, kappa, seed = seed)
      b <- rhs(n, 700L + seed)
      rd <- krylov_argmin(M, b, m)
      rp <- reference_minres(M, as.vector(b), m)
      w <- max(w, max(Mod(rd - rp)) / max(Mod(rd)))
    }
    cat(sprintf("%-10.0e %5d  %11.2e\n", kappa, m, w))
  }
}

## ------------------------------------------------------------------ part 6 --
## Complex hermitian, where no external reference is available: SciPy's minres
## casts a complex operator to real. The dense definition is the only reference
## this case has, which is a reason to have measured it in part 1.

cat("\n=== 6. complex hermitian indefinite against both references ===\n")
cat(sprintf("%-10s %5s  %11s  %11s\n", "kappa", "steps", "def", "pub"))

for (kappa in c(1e2, 1e6)) {
  for (m in c(1L, 2L, 4L, 6L, 8L)) {
    wd <- 0; wp <- 0
    for (seed in 1:3) {
      M <- hpd_prescribed(30, indef_spectrum(30, kappa), seed = seed)
      A <- linop(M, properties = c(hermitian = TRUE))
      set.seed(950L + seed)
      b <- matrix(complex(real = stats::rnorm(30), imaginary = stats::rnorm(30)), 30, 1L)
      mine <- as.vector(minres(A, b, tol = 0, maxit = m)$x)
      rd <- krylov_argmin(M, b, m)
      rp <- reference_minres(M, as.vector(b), m)
      wd <- max(wd, max(Mod(mine - rd)) / max(Mod(rd)))
      wp <- max(wp, max(Mod(mine - rp)) / max(Mod(rp)))
    }
    cat(sprintf("%-10.0e %5d  %11.2e  %11.2e\n", kappa, m, wd, wp))
  }
}
