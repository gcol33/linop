## Section 7.2, the reference eigensolver. Lanczos with full reorthogonalisation,
## locking, and explicit restarts.
##
## What it requires. A hermitian operator and a forward action, and nothing else:
## every apply below is mode "N". That makes it the second thing in the package,
## after GMRES and BiCGSTAB, that runs on an operator supplying no adjoint at all,
## and for the same structural reason -- a hermitian operator IS its own adjoint,
## so there is no second callback to ask for.
##
## What it does not do. There is no evidence minimum on the hermitian
## requirement, and that is a decision rather than an omission. Every other
## dispatcher in the package that applies one has somewhere to fall back to: auto
## refuses CG and takes MINRES, refuses MINRES and takes GMRES. Here there is one
## method and no non-hermitian fallback until Arnoldi arrives with the backend, so
## a minimum would not select between methods, it would only refuse to run. The
## declaration is carried into the certificate instead: the forward-error row
## records it in depends_on, so a bound resting on a bare user_declaration fails
## any requirement the declaration would have failed directly. The value gates the
## run and the evidence is reported, which is section 5.3 applied where it can
## still say something.
##
## The Rayleigh quotient is what gets reported, not the Ritz value. For a plain
## Lanczos run the two agree to rounding, since a Ritz value of a hermitian
## operator is the Rayleigh quotient of its own Ritz vector. Under shift-invert
## they do not: sigma + 1/theta inherits whatever the inner solve got wrong, and
## x^H A x / (x^H x) is the minimiser of ||A x - mu x|| over mu, so it is the value
## for which the certificate's exhibited perturbation is smallest. Measuring the
## eigenvalue on A rather than reading it off the transformed problem is the same
## discipline as the outer loop measuring b - A x, and it is what makes an inexact
## shift-invert reportable rather than merely fast.
##
## Arithmetic. alpha is a Rayleigh quotient and beta is a norm, so the projected
## tridiagonal is real and symmetric whether the operator is real or complex, and
## its eigenvectors are real. Only the basis carries complex storage.

## The wanted eigenvalues under shift-invert are the ones nearest sigma, which is
## the largest end of the transformed spectrum whatever the caller asked about the
## original one. Selecting by anything else would need a second transformation and
## would not mean what it said.
SIGMA_SELECT <- "LM"

#' Eigenvalues and eigenvectors of a hermitian operator
#'
#' Lanczos with full reorthogonalisation, run against the operator's action. The
#' result records what converged, what residual each pair reached, and what
#' argument establishes the bound on the eigenvalues.
#'
#' The eigenvalue reported for each pair is the Rayleigh quotient
#' `x^H A x / (x^H x)`, measured on `A` itself. It is the value that minimises
#' `||A x - mu x||`, so no other choice of `mu` certifies better.
#'
#' `target identity` is `not_checked`, always. A small residual proves the pair is
#' an approximate eigenpair; that it is the largest, smallest or nearest one
#' requires an inertia count, an enclosure or a separation bound, and matrix-free
#' there is generally none. See `vignette("solvers")`.
#'
#' @param A A square hermitian `linop`. No adjoint is required.
#' @param k How many eigenpairs.
#' @param which Which end of the spectrum: `"largest"` or `"smallest"` by
#'   magnitude, `"largest_algebraic"` or `"smallest_algebraic"` by value.
#'   RSpectra's `"LM"`, `"SM"`, `"LA"` and `"SA"` are accepted for the same four.
#' @param sigma A real shift. Given one, the eigenvalues nearest `sigma` are
#'   sought, by running the same recurrence on `(A - sigma I)^-1` through the
#'   package's own solvers. `which` does not apply and passing both is an error.
#' @param B Reserved for the generalized problem `A x = lambda B x`, which is not
#'   in this version.
#' @param tol Convergence tolerance on the backward error
#'   `||A x - mu x|| / (||A|| ||x||)`.
#' @param maxit Total Lanczos steps across all restarts.
#' @param ncv Subspace size per restart. Storage is `ncv` vectors.
#' @param v0 Starting vector, or `NULL` for a seeded random one.
#' @param seed Seed for the starting vector. The caller's random stream is
#'   restored afterwards.
#' @param method `"auto"` or `"lanczos"`. There is one method in this version and
#'   `"auto"` resolves to it, recording the reason in the certificate.
#' @param inner Arguments for the inner solve when `sigma` is given, such as
#'   `method` or `tol`.
#' @param floor_const `c` in the arithmetic floor of the certificate.
#' @param norm_control Arguments for the `||A||` estimate.
#' @return An object of class `linop_eigen`, with `values`, `vectors` and the
#'   `certificate`. A run that exhausts its budget returns the best pairs it
#'   reached and a `fail` certificate rather than an error, and it may hold fewer
#'   than `k` of them: nothing is padded.
#' @examples
#' n <- 40
#' L <- linop(function(X) rbind(X[-1, , drop = FALSE], 0) +
#'                        rbind(0, X[-n, , drop = FALSE]) - 2 * X,
#'            dim = c(n, n), properties = c(hermitian = TRUE))
#' fit <- eigs(L, k = 3, which = "largest_algebraic")
#' fit$values
#' @export
eigs <- function(A, k, which = "largest", sigma = NULL, B = NULL,
                 tol = 1e-10, maxit = NULL, ncv = NULL, v0 = NULL, seed = 1L,
                 method = "auto", inner = list(),
                 floor_const = SOLVE_FLOOR_CONST, norm_control = list()) {
  if (!is_linop(A)) stopf("eigs() expects a linop")
  n <- A$dim[1L]
  if (A$dim[1L] != A$dim[2L]) {
    stopf(paste0("eigs() needs a square operator; this one is %d x %d.\n",
                 "  The singular triplets of a rectangular operator are svds()."),
          A$dim[1L], A$dim[2L])
  }
  if (!is.null(B)) {
    stopf(paste0("eigs() does not take B in this version.\n",
                 "  A x = lambda B x is Lanczos in the B inner product, which is a different\n",
                 "  recurrence rather than a flag on this one: every inner product changes and a\n",
                 "  solve against B enters each step. It arrives with the backend of section 3."))
  }
  require_capability(A, "hermitian", "eigs")
  if (!is_scalar_string(method) || !method %in% c("auto", "lanczos")) {
    stopf("method is 'auto' or 'lanczos'; Arnoldi for a non-hermitian operator arrives with the backend")
  }
  which_code <- normalise_which(which, EIGEN_WHICH, "eigs")
  s <- eigen_setup(A, k, ncv, maxit, "eigs", n)
  k <- s$k; ncv <- s$ncv; maxit <- s$maxit

  ## The operator the recurrence runs on, which is A itself unless a shift makes
  ## it the solve. Selection happens on that operator's spectrum; measurement
  ## always happens on A.
  if (is.null(sigma)) {
    Op <- A
    sel <- which_code
    reason <- "hermitian operator, Lanczos with full reorthogonalisation"
  } else {
    if (!identical(which, "largest")) {
      stopf(paste0("sigma and which cannot both be given.\n",
                   "  A shift asks for the eigenvalues nearest sigma, which fixes the end of the\n",
                   "  transformed spectrum the recurrence works at. Drop which, or drop sigma and\n",
                   "  ask for an end of the spectrum of A directly."))
    }
    Op <- shift_invert_operator(A, sigma, inner)
    sel <- SIGMA_SELECT
    reason <- sprintf("hermitian operator, Lanczos on (A - %g I)^-1", sigma)
  }

  run <- lanczos_core(Op, A, k, sel, tol, maxit, ncv, v0, seed)

  norm_est <- do.call(norm2, c(list(A = A), norm_control))
  cert <- eigen_certificate(A, run$values, run$vectors, tol = tol,
                            norm_estimate = norm_est,
                            iterations = run$steps, maxit = maxit,
                            floor_const = floor_const,
                            norm_lower_bound = run$anorm_lb,
                            hermitian_evidence = cape(A, "hermitian"),
                            requested = k,
                            stop_reason = run$stop_reason, subject = "eigen")

  new_eigen_result(values = run$values, vectors = run$vectors,
                   certificate = cert, k = k, which = which,
                   method = if (is.null(sigma)) "lanczos" else "lanczos shift-invert",
                   iterations = run$steps, maxit = maxit, rounds = run$rounds,
                   converged = cert$values$converged, residuals = run$residuals,
                   dispatch = list(requested = method,
                                   chosen = if (is.null(sigma)) "lanczos" else "lanczos shift-invert",
                                   reason = reason))
}

## (A - sigma I)^-1 as an operator, built out of the package's own solvers.
##
## Section 4.2 has this the right way round: shift-invert holds a solve object,
## not an operator, and what makes it usable as one here is that the inner method
## is a fixed iteration to a fixed tolerance. Whatever it gets wrong shows up in
## the residual the outer loop measures on A, which is why nothing downstream has
## to trust it.
##
## Hermitian by construction given that A is: the shifted operator inherits the
## declaration through caps_sum() and caps_scale(), and its inverse is hermitian
## with it. That is also why the adjoint callback is the apply -- a hermitian
## operator is its own adjoint, so norm2() gets its power iteration for nothing.
shift_invert_operator <- function(A, sigma, inner) {
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma)) {
    stopf(paste0("sigma must be a single finite real number.\n",
                 "  A complex shift would leave A - sigma I non-hermitian, which this recurrence\n",
                 "  is not defined for."))
  }
  n <- A$dim[1L]
  shifted <- A - sigma * linop_eye(n, A$dtype)
  ## Tighter than the outer tolerance by default: what the inner solve leaves
  ## behind is amplified by the transformation, and the outer loop can only
  ## report it, not remove it.
  args <- utils::modifyList(list(method = "minres", tol = 1e-12), inner)
  apply_inverse <- function(X) {
    z <- do.call(solve, c(list(shifted, X), args))
    attr(z, "certificate") <- NULL
    z
  }
  new_linop("fun", c(n, n), A$dtype,
            do.call(new_caps, list(
              hermitian = capability(TRUE, ev_construction(
                depends_on = Filter(Negate(is.null),
                                    list(cape(shifted, "hermitian"))))))),
            list(apply = apply_inverse, adjoint = apply_inverse))
}

## The outer loop. A round builds a subspace, extracts, measures on A, locks what
## converged and restarts thickly from what did not.
##
## Locking is what full reorthogonalisation buys beyond stability: a converged
## eigenvector joins the set every later basis is orthogonalised against, so the
## next round searches the orthogonal complement and cannot return the same pair
## twice. A degenerate eigenvalue is the case that shows this is deflation rather
## than exclusion -- a second, orthogonal eigenvector for the same eigenvalue is
## still reachable, and is still correct to return.
##
## Restarting thickly rather than from one vector is not a refinement, it is what
## makes the method converge at all. A restart from the sum of the unconverged
## Ritz vectors throws away the subspace that produced them: on laplacian_1d(60)
## at ncv = 24, asking for four pairs, that spends 240 of a 300-iteration budget
## over 9 rounds, converges none of them and stalls at a backward error of 5.7e-4,
## where this one reaches 5.1e-13 in 96. Each round was rediscovering from one
## direction what the last one spent 24 applies learning.
##
## Wu and Simon's thick restart keeps the wanted Ritz vectors themselves as the
## start of the next basis. They span an approximately invariant subspace already,
## so the projected problem keeps them diagonal at the values they were extracted
## with, and the arrow coupling below is all that ties them to the new directions.
## Where a request fits in one round there is no restart to differ over and the
## two agree to the last digit, which is the control that confines the difference
## to this block; dev_notes/spikes/restart-comparison.R runs both.
lanczos_core <- function(Op, A, k, sel, tol, maxit, ncv, v0, seed) {
  n <- Op$dim[1L]
  eps <- .Machine$double.eps

  Qv <- NULL
  lock_theta <- numeric(0)
  lock_value <- numeric(0)
  lock_res <- numeric(0)

  v <- eigen_start_vector(n, Op$dtype, v0, seed)
  keep <- NULL
  steps <- 0L
  rounds <- 0L
  anorm_lb <- 0
  stop_reason <- NULL
  prev_best <- Inf
  stalled <- 0L

  ## What comes back when the budget runs out before k pairs converge: the best
  ## pairs the last round produced, with their measured residuals, so a failed
  ## run reports a result and a certificate rather than an error.
  best_theta <- numeric(0); best_value <- numeric(0)
  best_vec <- NULL; best_res <- numeric(0)

  repeat {
    nl <- if (is.null(Qv)) 0L else ncol(Qv)
    want <- k - nl
    if (want <= 0L) break
    if (steps >= maxit) { stop_reason <- "budget exhausted"; break }
    if (nl >= n) { stop_reason <- "the locked set spans the operator"; break }
    rounds <- rounds + 1L

    ## The subspace this round may occupy, and how much of the budget is left to
    ## fill it with. Retained vectors take room in the first and cost nothing
    ## against the second, which is the whole economy of a thick restart.
    m <- min(ncv - nl, n - nl)
    if (!is.null(keep)) {
      p <- min(ncol(keep$Y), m - 2L)
      keep <- if (p >= 1L) list(Y = keep$Y[, seq_len(p), drop = FALSE],
                                theta = keep$theta[seq_len(p)],
                                b = keep$b[seq_len(p)]) else NULL
    }

    ## The start of a round lives in the complement of what is locked. A start
    ## that lies entirely inside the locked space carries no new direction, and a
    ## fresh random one does.
    v <- orth_against(Qv, v)
    if (!is.null(keep)) v <- orth_against(keep$Y, v)
    nv <- col_norms(v)
    if (nv <= eps * sqrt(n)) {
      v <- eigen_random_vector(n, Op$dtype, seed, rounds)
      v <- orth_against(Qv, v)
      if (!is.null(keep)) v <- orth_against(keep$Y, v)
      nv <- col_norms(v)
      if (nv <= eps * sqrt(n)) {
        stop_reason <- "no direction left outside the locked subspace"
        break
      }
    }
    v <- v / nv

    run <- lanczos_run(Op, v, Qv, m, maxit - steps, keep)
    steps <- steps + run$applies
    if (run$applies == 0L) { stop_reason <- "the recurrence made no step"; break }

    ## The projected problem, real and symmetric whatever the operator is.
    ee <- eigen(run$Tm, symmetric = TRUE)
    ord <- ritz_order(ee$values, sel)
    take <- ord[seq_len(min(want, run$size))]
    theta <- ee$values[take]
    Xr <- run$V %*% ee$vectors[, take, drop = FALSE]

    ## Measured on A, in every case. Under a shift the recurrence has been
    ## working on a different operator entirely, and this is where that stops
    ## mattering.
    AX <- linop_apply(A, Xr, "N")
    xn <- col_norms(Xr)
    scale_x <- ifelse(xn > 0, xn, 1)
    mu <- col_dot(Xr, AX) / scale_x^2
    rr <- col_norms(AX - scale_cols(Xr, mu)) / scale_x

    ## Two lower bounds on ||A||, both read off applies this round was making
    ## anyway: ||A x|| / ||x|| for every x, and |x^H A x| / ||x||^2, which is a
    ## Rayleigh quotient and so lies inside the spectrum.
    anorm_lb <- max(anorm_lb, max(col_norms(AX) / scale_x), max(abs(mu)))

    conv <- rr <= tol * anorm_lb
    if (any(conv)) {
      Qv <- cbind(Qv, Xr[, conv, drop = FALSE])
      lock_theta <- c(lock_theta, theta[conv])
      lock_value <- c(lock_value, mu[conv])
      lock_res <- c(lock_res, rr[conv])
    }
    best_theta <- theta; best_value <- mu; best_vec <- Xr; best_res <- rr

    if (sum(conv) + nl >= k) break

    ## A round that locked nothing and did not improve on the previous one twice
    ## running has nothing left to recover. Two rather than one, because a thick
    ## restart can spend a round rebuilding the directions it dropped.
    best <- min(rr)
    if (!any(conv) && best >= prev_best) {
      stalled <- stalled + 1L
      if (stalled >= 2L) {
        stop_reason <- "the subspace stopped improving"
        break
      }
    } else {
      stalled <- 0L
    }
    prev_best <- min(prev_best, best)

    ## The thick restart. Everything the round produced that was not locked is
    ## available to retain, in the order the request puts it; the wanted ones plus
    ## about half the subspace, which is the usual setting. The coupling
    ##
    ##     y_i^H A v_next = beta_last * s_i[last]
    ##
    ## is exactly what ties each retained Ritz vector to the direction the run
    ## ended on, and it is the only off-diagonal the retained block needs.
    left_want <- want - sum(conv)
    avail <- setdiff(ord, take[conv])
    p <- min(length(avail), max(left_want + 3L, m %/% 2L), m - 2L)
    if (is.null(run$vnext) || p < 1L) {
      ## An invariant subspace, or no room to retain anything: there is no
      ## direction to continue from and the next round starts fresh.
      keep <- NULL
      v <- eigen_random_vector(n, Op$dtype, seed, rounds)
    } else {
      idx <- avail[seq_len(p)]
      Sk <- ee$vectors[, idx, drop = FALSE]
      keep <- list(Y = run$V %*% Sk, theta = ee$values[idx],
                   b = run$beta_last * Sk[run$size, ])
      v <- run$vnext
    }
  }

  ## What is returned: everything locked, then whatever the last round had left,
  ## up to k. An unconverged pair comes back with its residual and is reported as
  ## unconverged by the certificate rather than withheld.
  values <- lock_value
  vectors <- Qv
  residuals <- lock_res
  thetas <- lock_theta
  short <- k - length(values)
  if (short > 0L && length(best_value)) {
    keep <- seq_len(min(short, length(best_value)))
    values <- c(values, best_value[keep])
    vectors <- cbind(vectors, best_vec[, keep, drop = FALSE])
    residuals <- c(residuals, best_res[keep])
    thetas <- c(thetas, best_theta[keep])
  }
  if (!length(values)) {
    stopf("eigs() produced no pair at all; the operator applied to the starting vector is zero")
  }

  ord <- ritz_order(thetas, sel)
  list(values = values[ord], vectors = vectors[, ord, drop = FALSE],
       residuals = residuals[ord], steps = steps, rounds = rounds,
       anorm_lb = anorm_lb, stop_reason = stop_reason)
}

## One run of the three-term recurrence, from a unit v orthogonal to the locked
## set Q and to whatever is retained, filling a basis of at most m vectors and
## spending at most max_applies of them on the operator.
##
##     beta_j v_{j+1} = A v_j - alpha_j v_j - beta_{j-1} v_{j-1},
##     alpha_j = <v_j, A v_j>,   beta_j = ||...||,
##
## with every new direction orthogonalised against Q and against the whole of V
## afterwards. In exact arithmetic that second step subtracts nothing; in floating
## point it is the difference between a Lanczos basis and a set of vectors that
## have quietly become dependent, and the projected problem cannot tell.
##
## The projected matrix is built here rather than reassembled afterwards, because
## after a thick restart it is not tridiagonal. The retained block is diagonal,
## carrying the Ritz values it was extracted with, and it couples to the rest only
## through one row and column:
##
##     [ theta_1                b_1 ]
##     [        ...             ... ]
##     [            theta_p     b_p ]
##     [ b_1    ...     b_p   alpha ]  ->  tridiagonal from here on
##
## which is why the first new step subtracts Y b rather than beta v_{j-1}. From
## the second step the recurrence is the ordinary one again: the retained
## directions have already been removed and full reorthogonalisation keeps them
## removed.
##
## alpha is the real part of the inner product, which is the whole of it: a
## Rayleigh quotient of a hermitian operator is real, so col_dot() rather than
## col_cdot() is not an approximation here. The correction inside orth_against()
## is a different matter and uses the complex product, because the coefficient
## <v_i, w> for a general i has no reason to be real.
lanczos_run <- function(Op, v, Q, m, max_applies, keep = NULL) {
  n <- nrow(v)
  eps <- .Machine$double.eps
  p <- if (is.null(keep)) 0L else ncol(keep$Y)
  cplx <- is.complex(v) || (!is.null(keep) && is.complex(keep$Y)) ||
          Op$dtype == "complex"
  V <- matrix(if (cplx) complex(1) else 0, n, m)
  Tm <- matrix(0, m, m)
  if (p > 0L) {
    V[, seq_len(p)] <- keep$Y
    Tm[cbind(seq_len(p), seq_len(p))] <- keep$theta
    Tm[seq_len(p), p + 1L] <- keep$b
    Tm[p + 1L, seq_len(p)] <- keep$b
  }
  V[, p + 1L] <- v

  i <- p
  applies <- 0L
  beta_last <- 0
  vnext <- NULL

  repeat {
    i <- i + 1L
    applies <- applies + 1L
    vi <- V[, i, drop = FALSE]
    W <- linop_apply(Op, vi, "N")
    ## The size the direction had before anything was taken out of it, which is
    ## the scale the breakdown test is relative to. The same test GMRES makes,
    ## with the same constant, and for the same reason: what it guards is a
    ## division, and normalising by a norm at rounding level would put an
    ## overflowed vector into the basis.
    w0 <- col_norms(W)
    a <- col_dot(vi, W)
    Tm[i, i] <- a
    W <- W - scale_cols(vi, a)
    if (i == p + 1L) {
      if (p > 0L) W <- W - keep$Y %*% matrix(keep$b, p, 1L)
    } else {
      W <- W - scale_cols(V[, i - 1L, drop = FALSE], Tm[i - 1L, i])
    }
    W <- orth_against(Q, W)
    W <- orth_against(V[, seq_len(i), drop = FALSE], W)

    b <- col_norms(W)
    beta_last <- b
    ## An invariant subspace, which is convergence rather than failure: every
    ## Ritz pair of it is an exact eigenpair of the operator restricted to it.
    ## There is no direction left to restart from either, and vnext says so.
    if (b <= eps * w0) break
    if (i >= m || applies >= max_applies) { vnext <- W / b; break }
    Tm[i, i + 1L] <- b
    Tm[i + 1L, i] <- b
    V[, i + 1L] <- W / b
  }

  list(V = V[, seq_len(i), drop = FALSE],
       Tm = Tm[seq_len(i), seq_len(i), drop = FALSE],
       size = i, applies = applies, vnext = vnext, beta_last = beta_last)
}
