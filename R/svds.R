## Section 7.2, the reference singular value solver. Golub-Kahan-Lanczos with
## full reorthogonalisation, locking, and explicit restarts.
##
## It is not eigs() on A^H A, and the difference is the whole point of
## bidiagonalising. The eigenvalues of the Gram operator are the squares of the
## singular values, so a singular value at 1e-9 relative to the largest sits at
## 1e-18 there and is gone. The bidiagonal keeps them at their own scale.
##
## It is also not the same code as solvers-bidiag.R, and that is worth saying
## because the recurrence is. LSQR and LSMR build this bidiagonalisation and store
## none of it: neither method keeps a basis, which is why both lose orthogonality
## after a handful of steps and why dev_notes/lsqr-and-the-least-squares-
## certificate.md can only assert four. A singular triplet is read off the basis
## rather than off the projection alone, so here the basis has to be kept and kept
## orthogonal, and the reorthogonalisation that a least-squares solve has no use
## for is the load-bearing part. The step is shared in shape and not in code:
## bidiag_step() takes no basis to orthogonalise against and adding one would put
## storage the least-squares methods do not want into their inner loop.
##
## The certificate is the eigen certificate, literally. The augmented operator
##
##     H = [ 0    A  ]
##         [ A^H  0  ]
##
## is hermitian for every A, with no assumption to record, and its eigenpairs are
## (sigma_i, [u_i; v_i]/sqrt(2)) and (-sigma_i, [u_i; -v_i]/sqrt(2)). So a singular
## triplet is an eigenpair of a hermitian operator, its residual is that
## eigenpair's residual, and everything eigen_certificate() proves about a
## hermitian eigenpair holds here unchanged. One difference is worth noticing: the
## hermitian-ness of H is construction with an empty depends_on, so the forward
## bound on a singular value rests on nothing the caller declared, where the same
## bound from eigs() rests on whatever established the operator's own symmetry.

#' Singular values and vectors of an operator
#'
#' Golub-Kahan-Lanczos with full reorthogonalisation, run against the operator's
#' action and its adjoint. The result records what converged, what residual each
#' triplet reached, and what argument establishes the bound on the singular
#' values.
#'
#' The certificate is the certificate of the corresponding eigenpair of the
#' hermitian augmented operator, so `residual`, `backward error` and
#' `forward error` mean there exactly what they mean for [eigs()].
#'
#' `which = "smallest"` is accepted and is the hard end of the problem: a Krylov
#' space reaches the largest singular values first, so the smallest need a
#' subspace approaching the rank of the operator. The certificate reports what was
#' reached rather than the request being refused.
#'
#' @param A A `linop`. Both the forward action and the adjoint are needed.
#' @param k How many singular triplets.
#' @param which `"largest"` or `"smallest"`. RSpectra's `"LM"` and `"SM"` are
#'   accepted for the same two. The algebraic forms are not separate names here
#'   because a singular value is non-negative.
#' @param tol Convergence tolerance on the backward error and the relative
#'   residual of the augmented eigenpair, each plus the arithmetic floor: the
#'   certificate's `backward error` and `residual` rows.
#' @param maxit Total applies of the operator, in either mode, across all
#'   restarts. `NULL` takes `20 min(m, n)` for an `m x n` operator, the rule
#'   [eigs()] uses on the dimension the singular vectors are sought in.
#' @param ncv Subspace size per restart. Storage is `ncv` vectors of each
#'   length. `NULL` takes the engine's default.
#' @param v0 Starting vector in the domain, or `NULL` for a seeded random one.
#'   `"golub-kahan"` only.
#' @param seed Seed for the starting vector. The caller's random stream is
#'   restored afterwards.
#' @param method `"auto"`, `"golub-kahan"` or `"primme"`. `"golub-kahan"` is the
#'   reference implementation in this package; `"primme"` is the PRIMME library,
#'   available when the `linop.primme` package is installed. `"auto"` takes
#'   PRIMME when it is available and Golub-Kahan otherwise, and the
#'   certificate's `dispatch` records every engine considered and why. The
#'   reported singular values are Rayleigh quotients `Re(u^H A v)` of the
#'   returned unit vectors, measured on `A`, whichever engine ran.
#' @param floor_const `c` in the arithmetic floor of the certificate.
#' @param norm_control Arguments for the `||A||` estimate.
#' @return An object of class `linop_svd`, with `d`, `u`, `v` and the
#'   `certificate`. A run that exhausts its budget returns the best triplets it
#'   reached and a `fail` certificate rather than an error, and it may hold fewer
#'   than `k` of them: nothing is padded.
#' @examples
#' A <- linop(rbind(diag(c(3, 2, 1)), 0))
#' fit <- svds(A, k = 2)
#' fit$d
#' @export
svds <- function(A, k, which = "largest", tol = 1e-10, maxit = NULL,
                 ncv = NULL, v0 = NULL, seed = 1L, method = "auto",
                 floor_const = SOLVE_FLOOR_CONST, norm_control = list()) {
  if (!is_linop(A)) stopf("svds() expects a linop")
  require_finite_dim(A, "svds")
  m <- A$dim[1L]; n <- A$dim[2L]
  if (!expr_has_adjoint(A)) {
    stopf(paste0("svds() needs the adjoint as well as the forward action.\n",
                 "  A bidiagonalisation spends one apply in each mode per step, so an operator\n",
                 "  supplying only a forward action has no singular triplets reachable here.\n",
                 "  Supply adjoint = to linop()."))
  }
  which_code <- normalise_which(which, SVD_WHICH, "svds")
  s <- eigen_setup(A, k, maxit, "svds", min(m, n))
  k <- s$k; maxit <- s$maxit

  pick <- select_backend("svds", A, list(k = k, which = which_code, v0 = !is.null(v0),
                                    tol = tol, ncv = ncv),
                         method, "svds")
  ncv <- eigen_ncv(ncv %||% pick$backend$default_ncv(A, k), k, min(m, n))
  run <- pick$backend$run(A, list(k = k, sel = which_code, tol = tol, maxit = maxit,
                                 ncv = ncv, v0 = v0, seed = seed,
                                 floor_const = floor_const), pick$engine)

  cert <- svd_certificate(A, run$d, run$u, run$v, tol = tol, iterations = run$steps,
                          maxit = maxit, floor_const = floor_const,
                          norm_control = norm_control, norm_lower_bound = run$anorm_lb,
                          requested = k, stop_reason = run$stop_reason)

  new_svd_result(d = run$d, u = run$u, v = run$v, certificate = cert,
                 k = k, which = which, method = pick$backend$name,
                 iterations = run$steps, maxit = maxit, rounds = run$rounds,
                 converged = cert$values$converged, residuals = run$residuals,
                 dispatch = backend_dispatch(method, pick, pick$backend$label),
                 engine = run$engine)
}

## The hermitian operator whose eigenpairs are the singular triplets. Hermitian
## unconditionally: swapping the two blocks and conjugating gives back the same
## map for any A whatever, so the evidence is construction with nothing under it.
## Its own adjoint, for the same reason, so it needs no second callback.
augmented_operator <- function(A) {
  m <- A$dim[1L]; n <- A$dim[2L]
  top <- seq_len(m)
  bot <- m + seq_len(n)
  ap <- function(Z) {
    rbind(linop_apply(A, Z[bot, , drop = FALSE], "N"),
          linop_apply(A, Z[top, , drop = FALSE], "C"))
  }
  new_linop("fun", c(m + n, m + n), A$dtype,
            do.call(new_caps, list(hermitian = capability(TRUE, ev_construction()))),
            list(apply = ap, adjoint = ap))
}

## The certificate of a set of triplets: the eigen certificate of the augmented
## eigenproblem they are read from.
svd_certificate <- function(A, d, u, v, tol, iterations, maxit, floor_const,
                            norm_control, norm_lower_bound, requested, stop_reason) {
  ## The augmented eigenproblem the certificate is read from. z is unit because u
  ## and v are, so no rescaling enters the residual.
  H <- augmented_operator(A)
  Z <- rbind(u, v) / sqrt(2)

  ## ||H||_2 = ||A||_2 exactly: the spectrum of H is the singular values of A with
  ## both signs, plus zeros. A structural rule over an estimated child, so the
  ## estimate underneath stays visible at the top (section 5.3).
  inner <- do.call(norm2, c(list(A = A), norm_control))
  norm_est <- new_norm_estimate(
    inner$value, ev_construction(depends_on = list(inner$evidence)), inner$method,
    sprintf("the augmented operator has the singular values of A with both signs; [%s]",
            inner$detail))

  eigen_certificate(H, d, Z, tol = tol, norm_estimate = norm_est,
                    iterations = iterations, maxit = maxit,
                    floor_const = floor_const,
                    norm_lower_bound = norm_lower_bound,
                    hermitian_evidence = cape(H, "hermitian"),
                    requested = requested,
                    stop_reason = stop_reason, subject = "svd")
}

## The bidiagonalisation runs on whichever of A and A^H maps the smaller space
## into the larger. Started in a domain wider than the codomain, the right
## Krylov space picks up the start vector's component in null(A), which has
## dimension n - m and is not a singular direction: the m singular values of a
## wide A belong to triplets (sigma, u, v) with A^H u = sigma v, and no u is paired
## with a v in null(A). The projected matrix still has a zero singular value for
## it, the Ritz v satisfies A v = 0 to rounding and the Ritz u cannot satisfy
## A^H u = 0, so the smallest request returns that non-triplet and the
## certificate fails it: diff_1d(100), 99 x 100, one smallest, stopped at a
## backward error of 0.21 with d = 2.7e-17, where the smallest singular value is
## 2 sin(pi / 200) = 3.1e-2. Run on A^H, the start lives in R^m and the right
## space of A is reached only through A^H, which has no null(A) component.
##
## A caller's v0 is in the domain of A; on A^H the start is A v0, the first left
## vector the bidiagonalisation of A would have produced from it.
gkl_svds <- function(A, ctl) {
  if (A$dim[1L] >= A$dim[2L]) {
    return(gkl_core(A, ctl$k, ctl$sel, ctl$tol, ctl$maxit, ctl$ncv, ctl$v0,
                    ctl$seed, ctl$floor_const))
  }
  u0 <- if (is.null(ctl$v0)) NULL else linop_apply(A, as_block(ctl$v0), "N")
  run <- gkl_core(adjoint(A), ctl$k, ctl$sel, ctl$tol, ctl$maxit, ctl$ncv, u0,
                  ctl$seed, ctl$floor_const)
  u <- run$u
  run$u <- run$v
  run$v <- u
  run
}

## The outer loop, the same shape as lanczos_core(): build, extract, measure,
## lock what converged, restart thickly from what did not. Two bases instead of
## one and two orthogonalisation sets, because a left singular vector and a right
## one live in different spaces and a triplet converges in both at once or not at
## all.
gkl_core <- function(A, k, sel, tol, maxit, ncv, v0, seed, floor_const) {
  m <- A$dim[1L]; n <- A$dim[2L]
  eps <- .Machine$double.eps

  Ul <- NULL; Vl <- NULL
  lock_d <- numeric(0); lock_res <- numeric(0)

  v <- eigen_start_vector(n, A$dtype, v0, seed)
  keep <- NULL
  steps <- 0L
  rounds <- 0L
  anorm_lb <- 0
  stop_reason <- NULL
  progress <- NULL

  best_d <- numeric(0); best_u <- NULL; best_v <- NULL; best_res <- numeric(0)

  repeat {
    nl <- if (is.null(Vl)) 0L else ncol(Vl)
    want <- k - nl
    if (want <= 0L) break
    if (steps >= maxit) { stop_reason <- "budget exhausted"; break }
    if (nl >= min(m, n)) { stop_reason <- "the locked set spans the operator"; break }
    rounds <- rounds + 1L

    j <- min(ncv - nl, min(m, n) - nl)
    if (!is.null(keep)) {
      p <- min(ncol(keep$Yv), j - 2L)
      keep <- if (p >= 1L) list(Yu = keep$Yu[, seq_len(p), drop = FALSE],
                                Yv = keep$Yv[, seq_len(p), drop = FALSE],
                                H = keep$H[seq_len(p), seq_len(p), drop = FALSE])
              else NULL
    }

    v <- orth_against(Vl, v)
    if (!is.null(keep)) v <- orth_against(keep$Yv, v)
    nv <- col_norms(v)
    if (nv <= eps * sqrt(n)) {
      v <- eigen_random_vector(n, A$dtype, seed, rounds)
      v <- orth_against(Vl, v)
      if (!is.null(keep)) v <- orth_against(keep$Yv, v)
      nv <- col_norms(v)
      if (nv <= eps * sqrt(n)) {
        stop_reason <- "no direction left outside the locked subspace"
        break
      }
    }
    v <- v / nv

    run <- gkl_run(A, v, Ul, Vl, j, maxit - steps, keep)
    steps <- steps + run$applies
    if (run$size == 0L) { stop_reason <- "the recurrence made no step"; break }

    ## The projected problem, U^H A V as measured, real for a real operator.
    sv <- svd(run$Bm)
    ord <- ritz_order(sv$d, sel)
    take <- ord[seq_len(min(want, run$size))]
    sigma <- sv$d[take]
    Ur <- run$U %*% sv$u[, take, drop = FALSE]
    Vr <- run$V %*% sv$v[, take, drop = FALSE]

    ## Measured on A, in both directions. A triplet is converged when the
    ## augmented eigenpair is, which is the pair of residuals combined in the
    ## norm the augmented vector [u; v]/sqrt(2) is unit in.
    AV <- linop_apply(A, Vr, "N")
    AtU <- linop_apply(A, Ur, "C")
    ru <- col_norms(AV - scale_cols(Ur, sigma))
    rv <- col_norms(AtU - scale_cols(Vr, sigma))
    rr <- sqrt(ru^2 + rv^2) / sqrt(2)

    vn <- col_norms(Vr)
    anorm_lb <- max(anorm_lb, max(col_norms(AV) / ifelse(vn > 0, vn, 1)), max(sigma))

    ## The augmented pair's two certificate rows, as lanczos_core() locks on.
    conv <- eigen_pair_bounds(rr, sigma, anorm_lb, 1, tol, floor_const)$converged
    if (any(conv)) {
      Ul <- cbind(Ul, Ur[, conv, drop = FALSE])
      Vl <- cbind(Vl, Vr[, conv, drop = FALSE])
      lock_d <- c(lock_d, sigma[conv])
      lock_res <- c(lock_res, rr[conv])
    }
    ## The triplets this round did not lock, which are what a run ending here
    ## adds after the locked set; a locked triplet is already in it.
    open <- !conv
    best_d <- sigma[open]; best_u <- Ur[, open, drop = FALSE]
    best_v <- Vr[, open, drop = FALSE]; best_res <- rr[open]

    if (sum(conv) + nl >= k) break

    progress <- restart_progress(progress, rr, sigma, anorm_lb, floor_const, conv)
    if (progress$stalled >= 2L) {
      stop_reason <- "the subspace stopped improving"
      break
    }

    ## The thick restart, in the augmented form. In exact arithmetic the retained
    ## triplets stay diagonal in the projected matrix and couple to the direction
    ## the run ended on through one column,
    ##
    ##     u_i^H A v_next = beta_last * s_u,i[last],
    ##
    ## which follows from A^H U = V B^H + beta_last v_next e^H and v_next being
    ## orthogonal to V. The retained block is carried as Z_u^H B Z_v and the
    ## coupling is measured by the next run, for the reason lanczos_run() gives.
    left_want <- want - sum(conv)
    avail <- setdiff(ord, take[conv])
    p <- min(length(avail), max(left_want + 3L, j %/% 2L), j - 2L)
    if (is.null(run$vnext) || p < 1L) {
      keep <- NULL
      v <- eigen_random_vector(n, A$dtype, seed, rounds)
    } else {
      idx <- avail[seq_len(p)]
      Zu <- sv$u[, idx, drop = FALSE]
      Zv <- sv$v[, idx, drop = FALSE]
      keep <- list(Yu = run$U %*% Zu, Yv = run$V %*% Zv,
                   H = cross_adjoint(Zu, run$Bm %*% Zv))
      v <- run$vnext
    }
  }

  d <- lock_d; U <- Ul; V <- Vl; residuals <- lock_res
  short <- k - length(d)
  if (short > 0L && length(best_d)) {
    keep <- seq_len(min(short, length(best_d)))
    d <- c(d, best_d[keep])
    U <- cbind(U, best_u[, keep, drop = FALSE])
    V <- cbind(V, best_v[, keep, drop = FALSE])
    residuals <- c(residuals, best_res[keep])
  }
  if (!length(d)) stopf("svds() produced no triplet at all; the operator is zero on the starting vector")

  ord <- ritz_order(d, sel)
  list(d = d[ord], u = U[, ord, drop = FALSE], v = V[, ord, drop = FALSE],
       residuals = residuals[ord], steps = steps, rounds = rounds,
       anorm_lb = anorm_lb, stop_reason = stop_reason)
}

## One run of the bidiagonalisation, from a unit v orthogonal to the locked right
## vectors and to whatever is retained, filling bases of at most j vectors and
## spending at most max_applies steps.
##
##     alpha_1 u_1     = A v_1,
##     beta_i v_{i+1}  = A^H u_i - alpha_i v_i,
##     alpha_{i+1} u_{i+1} = A v_{i+1} - beta_i u_i,
##
## so A V = U B with B upper bidiagonal, alpha on the diagonal and beta above it.
## Each new left vector is orthogonalised against the locked left vectors and the
## whole of U, each new right vector against the locked right vectors and the
## whole of V. Two applies per step, one in each mode, which is the same cost the
## least-squares methods pay for the same recurrence without the storage.
##
## In exact arithmetic, after a thick restart, the projection is a broken arrow
## rather than a bidiagonal: the retained triplets are diagonal, the column that
## follows them carries the couplings rho, and from the step after that it is
## bidiagonal again.
##
##     [ sigma_1               rho_1              ]
##     [        ...            ...                ]
##     [            sigma_p    rho_p              ]
##     [                       alpha   beta       ]
##     [                               alpha  ... ]
##
## What is stored is B = U^H A V as the orthogonalisations measure it: the
## column above each new left vector from its projection against U, the row left
## of each new right vector from its projection against V. That is the matrix
## above to rounding, and the retained block as it is rather than as exact
## arithmetic would have it.
gkl_run <- function(A, v, Ul, Vl, j, max_applies, keep = NULL) {
  m <- A$dim[1L]; n <- A$dim[2L]
  eps <- .Machine$double.eps
  p <- if (is.null(keep)) 0L else ncol(keep$Yv)
  cplx <- is.complex(v) || A$dtype == "complex" ||
          (!is.null(keep) && is.complex(keep$Yv))
  zed <- if (cplx) complex(1) else 0
  U <- matrix(zed, m, j)
  V <- matrix(zed, n, j)
  Bm <- matrix(zed, j, j)
  if (p > 0L) {
    U[, seq_len(p)] <- keep$Yu
    V[, seq_len(p)] <- keep$Yv
    Bm[seq_len(p), seq_len(p)] <- keep$H
  }
  V[, p + 1L] <- v

  ## The left vector matching the new right one. Under a restart its projection
  ## on the retained left vectors is the arrow column above.
  W <- linop_apply(A, V[, p + 1L, drop = FALSE], "N")
  w0 <- col_norms(W)
  W <- orth_against(Ul, W)
  if (p > 0L) {
    pr <- orth_project(U[, seq_len(p), drop = FALSE], W)
    W <- pr$W
    Bm[seq_len(p), p + 1L] <- pr$h[, 1L]
  }
  a <- col_norms(W)
  if (a <= eps * w0) {
    return(list(U = U, V = V, Bm = Bm, size = 0L, applies = 1L,
                vnext = NULL, beta_last = 0))
  }
  Bm[p + 1L, p + 1L] <- a
  U[, p + 1L] <- W / a

  i <- p + 1L
  applies <- 1L
  beta_last <- 0
  vnext <- NULL

  repeat {
    ## The direction the next right vector would come from, computed whether or
    ## not there is room for it: it is what a thick restart continues from.
    Z <- linop_apply(A, U[, i, drop = FALSE], "C")
    applies <- applies + 1L
    z0 <- col_norms(Z)
    Z <- Z - scale_cols(V[, i, drop = FALSE], Conj(Bm[i, i]))
    Z <- orth_against(Vl, Z)
    pr <- orth_project(V[, seq_len(i), drop = FALSE], Z)
    Z <- pr$W
    g <- Conj(pr$h[, 1L])
    g[i] <- g[i] + Bm[i, i]
    Bm[i, seq_len(i)] <- g
    b <- col_norms(Z)
    beta_last <- b
    if (b <= eps * z0) break
    vnext <- Z / b
    if (i >= j || applies >= max_applies) break

    V[, i + 1L] <- vnext
    Bm[i, i + 1L] <- b
    W <- linop_apply(A, V[, i + 1L, drop = FALSE], "N")
    applies <- applies + 1L
    w0 <- col_norms(W)
    W <- W - scale_cols(U[, i, drop = FALSE], b)
    W <- orth_against(Ul, W)
    pr <- orth_project(U[, seq_len(i), drop = FALSE], W)
    W <- pr$W
    h <- pr$h[, 1L]
    h[i] <- h[i] + b
    Bm[seq_len(i), i + 1L] <- h
    a <- col_norms(W)
    if (a <= eps * w0) break
    Bm[i + 1L, i + 1L] <- a
    U[, i + 1L] <- W / a
    i <- i + 1L
  }

  Bm <- Bm[seq_len(i), seq_len(i), drop = FALSE]
  if (!cplx) Bm <- Re(Bm)
  list(U = U[, seq_len(i), drop = FALSE], V = V[, seq_len(i), drop = FALSE],
       Bm = Bm,
       size = i, applies = applies, vnext = vnext, beta_last = beta_last)
}
