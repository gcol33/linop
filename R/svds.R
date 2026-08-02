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
#' @param tol Convergence tolerance on the backward error of the augmented
#'   eigenpair.
#' @param maxit Total bidiagonalisation steps across all restarts.
#' @param ncv Subspace size per restart. Storage is `ncv` vectors of each length.
#' @param v0 Starting vector in the domain, or `NULL` for a seeded random one.
#' @param seed Seed for the starting vector. The caller's random stream is
#'   restored afterwards.
#' @param method `"auto"` or `"golub-kahan"`. There is one method in this version
#'   and `"auto"` resolves to it, recording the reason in the certificate.
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
  if (!is_scalar_string(method) || !method %in% c("auto", "golub-kahan")) {
    stopf("method is 'auto' or 'golub-kahan'")
  }
  which_code <- normalise_which(which, SVD_WHICH, "svds")
  s <- eigen_setup(A, k, ncv, maxit, "svds", min(m, n))
  k <- s$k; ncv <- s$ncv; maxit <- s$maxit

  run <- gkl_core(A, k, which_code, tol, maxit, ncv, v0, seed)

  ## The augmented eigenproblem the certificate is read from. z is unit because u
  ## and v are, so no rescaling enters the residual.
  H <- augmented_operator(A)
  Z <- rbind(run$u, run$v) / sqrt(2)

  ## ||H||_2 = ||A||_2 exactly: the spectrum of H is the singular values of A with
  ## both signs, plus zeros. A structural rule over an estimated child, so the
  ## estimate underneath stays visible at the top (section 5.3).
  inner <- do.call(norm2, c(list(A = A), norm_control))
  norm_est <- new_norm_estimate(
    inner$value, ev_construction(depends_on = list(inner$evidence)), inner$method,
    sprintf("the augmented operator has the singular values of A with both signs; [%s]",
            inner$detail))

  cert <- eigen_certificate(H, run$d, Z, tol = tol, norm_estimate = norm_est,
                            iterations = run$steps, maxit = maxit,
                            floor_const = floor_const,
                            norm_lower_bound = run$anorm_lb,
                            hermitian_evidence = cape(H, "hermitian"),
                            requested = k,
                            stop_reason = run$stop_reason, subject = "svd")

  new_svd_result(d = run$d, u = run$u, v = run$v, certificate = cert,
                 k = k, which = which, method = "golub-kahan",
                 iterations = run$steps, maxit = maxit, rounds = run$rounds,
                 converged = cert$values$converged, residuals = run$residuals,
                 dispatch = list(requested = method, chosen = "golub-kahan",
                                 reason = "Golub-Kahan-Lanczos with full reorthogonalisation"))
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

## The outer loop, the same shape as lanczos_core(): build, extract, measure,
## lock what converged, restart thickly from what did not. Two bases instead of
## one and two orthogonalisation sets, because a left singular vector and a right
## one live in different spaces and a triplet converges in both at once or not at
## all.
gkl_core <- function(A, k, sel, tol, maxit, ncv, v0, seed) {
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
  prev_best <- Inf
  stalled <- 0L

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
                                sigma = keep$sigma[seq_len(p)],
                                rho = keep$rho[seq_len(p)]) else NULL
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

    ## The projected problem is real whatever A is, since alpha, beta and the
    ## restart couplings are all norms. Only the two bases carry complex storage.
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

    conv <- rr <= tol * anorm_lb
    if (any(conv)) {
      Ul <- cbind(Ul, Ur[, conv, drop = FALSE])
      Vl <- cbind(Vl, Vr[, conv, drop = FALSE])
      lock_d <- c(lock_d, sigma[conv])
      lock_res <- c(lock_res, rr[conv])
    }
    best_d <- sigma; best_u <- Ur; best_v <- Vr; best_res <- rr

    if (sum(conv) + nl >= k) break

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

    ## The thick restart, in the augmented form: the retained triplets stay
    ## diagonal in the projected matrix and couple to the direction the run ended
    ## on through one column,
    ##
    ##     u_i^H A v_next = beta_last * s_u,i[last],
    ##
    ## which follows from A^H U = V B^H + beta_last v_next e^H and v_next being
    ## orthogonal to V. That is the broken arrow the next round starts from.
    left_want <- want - sum(conv)
    avail <- setdiff(ord, take[conv])
    p <- min(length(avail), max(left_want + 3L, j %/% 2L), j - 2L)
    if (is.null(run$vnext) || p < 1L) {
      keep <- NULL
      v <- eigen_random_vector(n, A$dtype, seed, rounds)
    } else {
      idx <- avail[seq_len(p)]
      keep <- list(Yu = run$U %*% sv$u[, idx, drop = FALSE],
                   Yv = run$V %*% sv$v[, idx, drop = FALSE],
                   sigma = sv$d[idx],
                   rho = run$beta_last * sv$u[run$size, idx])
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
## After a thick restart the projection is a broken arrow rather than a
## bidiagonal: the retained triplets are diagonal, and the column that follows
## them carries the couplings rho. From the step after that it is bidiagonal
## again.
##
##     [ sigma_1               rho_1              ]
##     [        ...            ...                ]
##     [            sigma_p    rho_p              ]
##     [                       alpha   beta       ]
##     [                               alpha  ... ]
gkl_run <- function(A, v, Ul, Vl, j, max_applies, keep = NULL) {
  m <- A$dim[1L]; n <- A$dim[2L]
  eps <- .Machine$double.eps
  p <- if (is.null(keep)) 0L else ncol(keep$Yv)
  cplx <- is.complex(v) || A$dtype == "complex" ||
          (!is.null(keep) && is.complex(keep$Yv))
  zed <- if (cplx) complex(1) else 0
  U <- matrix(zed, m, j)
  V <- matrix(zed, n, j)
  Bm <- matrix(0, j, j)
  if (p > 0L) {
    U[, seq_len(p)] <- keep$Yu
    V[, seq_len(p)] <- keep$Yv
    Bm[cbind(seq_len(p), seq_len(p))] <- keep$sigma
    Bm[seq_len(p), p + 1L] <- keep$rho
  }
  V[, p + 1L] <- v

  ## The left vector matching the new right one. Under a restart it has the
  ## retained couplings removed from it, which is the arrow column above.
  W <- linop_apply(A, V[, p + 1L, drop = FALSE], "N")
  w0 <- col_norms(W)
  if (p > 0L) W <- W - keep$Yu %*% matrix(keep$rho, p, 1L)
  W <- orth_against(Ul, W)
  if (p > 0L) W <- orth_against(U[, seq_len(p), drop = FALSE], W)
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
    Z <- Z - scale_cols(V[, i, drop = FALSE], Bm[i, i])
    Z <- orth_against(Vl, Z)
    Z <- orth_against(V[, seq_len(i), drop = FALSE], Z)
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
    W <- orth_against(U[, seq_len(i), drop = FALSE], W)
    a <- col_norms(W)
    if (a <= eps * w0) break
    Bm[i + 1L, i + 1L] <- a
    U[, i + 1L] <- W / a
    i <- i + 1L
  }

  list(U = U[, seq_len(i), drop = FALSE], V = V[, seq_len(i), drop = FALSE],
       Bm = Bm[seq_len(i), seq_len(i), drop = FALSE],
       size = i, applies = applies, vnext = vnext, beta_last = beta_last)
}
