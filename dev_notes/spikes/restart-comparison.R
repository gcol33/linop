## Does the thick restart earn its place, or would a restart from one vector do?
##
## The only difference between the two runs below is the block that chooses what
## the next round starts from. Everything else -- locking, full
## reorthogonalisation, the stall counter, the budget, the seed -- is the shipped
## code, reached by putting a modified lanczos_core() into the namespace and
## taking it out again.
##
## simple: the next round starts from the sum of the unconverged Ritz vectors and
##         nothing else, which is what the first draft of eigs() did.
## thick:  the shipped version, which retains those Ritz vectors themselves and
##         couples them to the direction the run ended on.

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

simple_restart_core <- function(Op, A, k, sel, tol, maxit, ncv, v0, seed) {
  n <- Op$dim[1L]
  eps <- .Machine$double.eps
  Qv <- NULL
  lock_theta <- numeric(0); lock_value <- numeric(0); lock_res <- numeric(0)
  v <- linop:::eigen_start_vector(n, Op$dtype, v0, seed)
  steps <- 0L; rounds <- 0L; anorm_lb <- 0
  stop_reason <- NULL; prev_best <- Inf; stalled <- 0L
  best_theta <- numeric(0); best_value <- numeric(0)
  best_vec <- NULL; best_res <- numeric(0)

  repeat {
    nl <- if (is.null(Qv)) 0L else ncol(Qv)
    want <- k - nl
    if (want <= 0L) break
    if (steps >= maxit) { stop_reason <- "budget exhausted"; break }
    if (nl >= n) { stop_reason <- "the locked set spans the operator"; break }
    rounds <- rounds + 1L

    m <- min(ncv - nl, n - nl)
    v <- linop:::orth_against(Qv, v)
    nv <- linop:::col_norms(v)
    if (nv <= eps * sqrt(n)) {
      v <- linop:::orth_against(Qv, linop:::eigen_random_vector(n, Op$dtype, seed, rounds))
      nv <- linop:::col_norms(v)
      if (nv <= eps * sqrt(n)) { stop_reason <- "no direction left"; break }
    }
    v <- v / nv

    ## keep = NULL every round: this is the whole difference.
    run <- linop:::lanczos_run(Op, v, Qv, m, maxit - steps, NULL)
    steps <- steps + run$applies
    if (run$applies == 0L) { stop_reason <- "no step"; break }

    ee <- eigen(run$Tm, symmetric = TRUE)
    ord <- linop:::ritz_order(ee$values, sel)
    take <- ord[seq_len(min(want, run$size))]
    theta <- ee$values[take]
    Xr <- run$V %*% ee$vectors[, take, drop = FALSE]

    AX <- linop:::linop_apply(A, Xr, "N")
    xn <- linop:::col_norms(Xr)
    scale_x <- ifelse(xn > 0, xn, 1)
    mu <- linop:::col_dot(Xr, AX) / scale_x^2
    rr <- linop:::col_norms(AX - linop:::scale_cols(Xr, mu)) / scale_x
    anorm_lb <- max(anorm_lb, max(linop:::col_norms(AX) / scale_x), max(abs(mu)))

    conv <- rr <= tol * anorm_lb
    if (any(conv)) {
      Qv <- cbind(Qv, Xr[, conv, drop = FALSE])
      lock_theta <- c(lock_theta, theta[conv])
      lock_value <- c(lock_value, mu[conv])
      lock_res <- c(lock_res, rr[conv])
    }
    best_theta <- theta; best_value <- mu; best_vec <- Xr; best_res <- rr
    if (sum(conv) + nl >= k) break

    best <- min(rr)
    if (!any(conv) && best >= prev_best) {
      stalled <- stalled + 1L
      if (stalled >= 2L) { stop_reason <- "the subspace stopped improving"; break }
    } else {
      stalled <- 0L
    }
    prev_best <- min(prev_best, best)

    left <- !conv
    v <- if (any(left)) Xr[, left, drop = FALSE] %*% matrix(1, sum(left), 1L)
         else linop:::eigen_random_vector(n, Op$dtype, seed, rounds)
  }

  values <- lock_value; vectors <- Qv
  residuals <- lock_res; thetas <- lock_theta
  short <- k - length(values)
  if (short > 0L && length(best_value)) {
    kp <- seq_len(min(short, length(best_value)))
    values <- c(values, best_value[kp])
    vectors <- cbind(vectors, best_vec[, kp, drop = FALSE])
    residuals <- c(residuals, best_res[kp])
    thetas <- c(thetas, best_theta[kp])
  }
  ordf <- linop:::ritz_order(thetas, sel)
  list(values = values[ordf], vectors = vectors[, ordf, drop = FALSE],
       residuals = residuals[ordf], steps = steps, rounds = rounds,
       anorm_lb = anorm_lb, stop_reason = stop_reason)
}

shipped <- linop:::lanczos_core

report <- function(label, fit, anorm) {
  cat(sprintf("%-8s nconv %d/%d  iterations %3d  restarts %d  worst backward %.3e  %s\n",
              label, fit$nconv, fit$k, fit$iterations, fit$restarts,
              max(fit$certificate$values$backward_error),
              fit$certificate$checks$detail[
                fit$certificate$checks$check == "convergence"]))
}

cases <- list(
  list(label = "laplacian_1d(60), 4 smallest algebraic, ncv 24",
       A = laplacian_1d(60), k = 4L, which = "smallest_algebraic", ncv = 24L),
  list(label = "laplacian_1d(60), 4 largest algebraic, ncv 24",
       A = laplacian_1d(60), k = 4L, which = "largest_algebraic", ncv = 24L),
  list(label = "spd_prescribed(40), 3 largest, ncv 12",
       A = as_spd_linop(spd_prescribed(40, c(100, 99.5, 40, 10,
                                             seq(1, 0.1, length.out = 36)), seed = 7L)),
       k = 3L, which = "largest_algebraic", ncv = 12L)
)

for (cs in cases) {
  cat("\n##", cs$label, "\n")
  assignInNamespace("lanczos_core", simple_restart_core, ns = "linop")
  a <- eigs(cs$A, k = cs$k, which = cs$which, ncv = cs$ncv, maxit = 300L)
  assignInNamespace("lanczos_core", shipped, ns = "linop")
  b <- eigs(cs$A, k = cs$k, which = cs$which, ncv = cs$ncv, maxit = 300L)
  report("simple", a)
  report("thick", b)
}
