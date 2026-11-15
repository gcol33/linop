## Section 7.2's eigensolver for an operator that is not established hermitian:
## Arnoldi with Krylov-Schur restarting and locking.
##
## What it requires. A forward action and nothing else, so it is to eigs() what
## GMRES is to solve(): the method that has no declaration to contradict, and
## which method = "auto" therefore falls back to for any square operator.
##
## The restart is Stewart's (A Krylov-Schur algorithm for large eigenproblems,
## SIAM J. Matrix Anal. Appl. 23, 2001), the default of SLEPc (EPSKRYLOVSCHUR)
## and of KrylovKit.jl's eigsolve. ARPACK's implicit restart filters the same
## subspace with exact shifts and reaches the same retained space in exact
## arithmetic; Stewart's form keeps it by an ordered Schur decomposition of the
## projected matrix instead, which is easier to make robust and is what makes
## locking a matter of not touching a block. The method keeps a Krylov
## decomposition
##
##     A V = V M + v b^H,
##
## V orthonormal and M general. Expansion appends Arnoldi steps; the restart
## takes the Schur form of M with the wanted Ritz values leading, and keeps that
## leading invariant subspace, so the retained block of the next M is the
## leading block of a triangular matrix and b is the last row of the Schur basis.
## For the hermitian case that is exactly the thick restart lanczos_core() runs.
##
## Locking. A converged pair's Schur vectors stay at the front of V and its
## entries of b are set to zero, which perturbs the decomposition by no more than
## that block's residual, so M stays block upper triangular with the locked block
## frozen. Every later direction is still orthogonalised against the locked
## vectors, and the eigenvector of a later pair is formed from the whole of M,
## because for a non-normal operator it has a component in the locked space. The
## locked block is never changed again, so the vectors measured when a pair
## locked are still eigenvectors of the final M.
##
## Real arithmetic for a real operator. The Schur form of a real M is complex,
## but its leading invariant subspace has a real basis whenever the chosen set is
## closed under conjugation, and that real basis is what the next round keeps:
## the real and imaginary parts of the complex Schur vectors span it, and their
## leading left singular vectors are an orthonormal basis of it. The cuts between
## locked, kept and discarded never separate a conjugate pair; which Schur
## eigenvalues are real is decided by the count LAPACK's real driver returns, so
## no threshold on an imaginary part enters.

arnoldi_core <- function(Op, A, k, sel, tol, maxit, ncv, v0, seed, floor_const,
                         hermitian = FALSE) {
  n <- Op$dim[1L]
  eps <- .Machine$double.eps

  v <- eigen_start_vector(n, Op$dtype, v0, seed)
  cplx <- Op$dtype == "complex" || is.complex(v)
  m <- min(ncv, n)

  ## The Krylov decomposition: V, M and the row b, with l locked columns leading.
  V <- NULL
  M <- NULL
  brow <- NULL
  l <- 0L

  lock_theta <- complex(0); lock_value <- complex(0); lock_res <- numeric(0)
  lock_vec <- NULL

  best_theta <- complex(0); best_value <- complex(0)
  best_vec <- NULL; best_res <- numeric(0)

  steps <- 0L
  rounds <- 0L
  anorm_lb <- 0
  ## Every expansion step computes ||A v_j|| for a unit v_j, a lower bound on
  ## ||A|| that costs nothing, and the stopping rule's floor scales with that
  ## bound. Measured pairs alone give a poor one at an end of the spectrum near
  ## zero: on convection-diffusion at n = 1000 the eigenvectors at the largest
  ## real part have ||A x|| about 1e-4 ||A||, the floor shrank with it, and a
  ## pair at a residual of 1.3e-14 never locked. Under a shift the iteration
  ## runs on the solve and its applies say nothing about A.
  op_is_a <- identical(Op, A)
  stop_reason <- NULL
  progress <- NULL

  repeat {
    want <- k - length(lock_value)
    if (want <= 0L) break
    if (steps >= maxit) { stop_reason <- "budget exhausted"; break }
    if (l >= n) { stop_reason <- "the locked set spans the operator"; break }
    rounds <- rounds + 1L

    ex <- arnoldi_expand(Op, V, M, brow, v, m, maxit - steps, cplx)
    steps <- steps + ex$applies
    if (ex$applies == 0L) { stop_reason <- "the recurrence made no step"; break }
    if (op_is_a) anorm_lb <- max(anorm_lb, ex$apply_norm)
    V <- ex$V
    M <- ex$M
    j <- ncol(V)
    act <- (l + 1L):j
    ma <- length(act)

    ## The Schur form of the active block, with every Ritz value the request
    ## ranks inside the wanted count leading, and their conjugates beside them.
    Ma <- M[act, act, drop = FALSE]
    sch <- schur_complex(Ma)
    lam <- diag(sch$T)
    real_ev <- schur_real_flags(Ma, lam, cplx)
    partner <- conjugate_partners(lam, real_ev)
    ord <- ritz_order(lam, sel)
    q <- closed_cut(ord, min(want, ma), partner, grow = TRUE)
    ord <- q$ord
    q <- q$cut
    re <- schur_reorder(sch$T, sch$U, ord[seq_len(q)])
    T <- re$T
    U <- re$U
    lam <- diag(T)
    real_ev <- real_ev[re$pos]
    partner <- match(partner[re$pos], re$pos)

    ## The Ritz pairs at the leading q positions, eigenvectors of the whole of M.
    smin <- eps * max(Mod(T), if (l) max(Mod(M[seq_len(l), seq_len(l)])) else 0,
                      .Machine$double.xmin / eps)
    Y <- vapply(seq_len(q), function(i) schur_eigvec(T, i, smin), complex(ma))
    Y <- matrix(Y, ma, q)
    Za <- U %*% Y
    if (l) {
      Tl <- M[seq_len(l), seq_len(l), drop = FALSE]
      Xc <- M[seq_len(l), act, drop = FALSE] %*% Za
      Zl <- vapply(seq_len(q), function(i)
                     as.vector(floored_shift_solve(Tl, lam[i], -Xc[, i, drop = FALSE], smin)),
                   complex(l))
      Z <- rbind(matrix(Zl, l, q), Za)
    } else {
      Z <- Za
    }
    Xr <- V %*% Z
    theta <- lam[seq_len(q)]
    if (!cplx) {
      ## A real eigenvalue of a real matrix has a real eigenvector, which the
      ## complex Schur form returns times a unit phase. Removing the phase of
      ## the largest entry and the rounding left in the imaginary part gives it.
      for (i in which(real_ev[seq_len(q)])) {
        z <- Xr[, i]
        p <- z[which.max(Mod(z))]
        Xr[, i] <- Re(z * Conj(p) / Mod(p))
        theta[i] <- Re(theta[i])
      }
    }
    Xr <- scale_cols(Xr, 1 / col_norms(Xr))

    meas <- eigen_measure(A, Xr, hermitian = hermitian)
    mu <- meas$mu
    rr <- meas$rr
    anorm_lb <- max(anorm_lb, meas$anorm_lb)
    conv <- eigen_pair_bounds(rr, mu, anorm_lb, 1, tol, floor_const)$converged
    ## A conjugate pair locks together or not at all.
    pq <- partner[seq_len(q)]
    has <- !is.na(pq) & pq <= q
    conv[has] <- conv[has] & conv[pq[has]]
    conv[!is.na(pq) & pq > q] <- FALSE
    ## Only the leading run of converged pairs locks, SLEPc's rule. A pair the
    ## request ranks lower may converge first, and locking it would make it part
    ## of the answer before a better-ranked value has had the chance to appear: on
    ## random 120 x 120 matrices, six eigenvalues at the largest magnitude or real
    ## part, locking every converged pair returned a set that was not the leading
    ## six in 5 of 60 runs.
    conv <- cumprod(conv) == 1
    if (any(conv)) {
      cl <- closed_cut(seq_len(q), sum(conv), partner, grow = FALSE)$cut
      conv <- seq_len(q) <= cl
    }
    c_new <- sum(conv)
    if (c_new) {
      lock_vec <- cbind(lock_vec, Xr[, conv, drop = FALSE])
      lock_theta <- c(lock_theta, theta[conv])
      lock_value <- c(lock_value, mu[conv])
      lock_res <- c(lock_res, rr[conv])
    }
    open <- !conv
    best_theta <- theta[open]; best_value <- mu[open]
    best_vec <- Xr[, open, drop = FALSE]; best_res <- rr[open]

    if (length(lock_value) >= k) break

    progress <- restart_progress(progress, rr, mu, anorm_lb, floor_const, conv)
    if (progress$stalled >= 2L) {
      stop_reason <- "the subspace stopped improving"
      break
    }

    ## The restart: the newly locked Schur vectors first, then the unconverged
    ## wanted ones and as many more as the rule lanczos_core() measured for its
    ## thick restart keeps, then nothing.
    left_want <- want - c_new
    rest <- setdiff(seq_len(ma), which(conv))
    rest <- rest[ritz_order(lam[rest], sel)]
    room <- m - l - c_new - 2L
    p_u <- max(0L, min(length(rest), max(left_want + 3L, ma %/% 2L), room))
    front <- c(which(conv), rest)
    kept <- closed_cut(front, c_new + p_u, partner, grow = FALSE)$cut
    re2 <- schur_reorder(T, U, front[seq_len(kept)])
    beta <- ex$beta
    if (kept == 0L) {
      V <- if (l) V[, seq_len(l), drop = FALSE] else NULL
      M <- if (l) M[seq_len(l), seq_len(l), drop = FALSE] else NULL
      brow <- if (l) numeric(l) else NULL
    } else {
      Q <- if (cplx) re2$U[, seq_len(kept), drop = FALSE]
           else real_nested_basis(re2$U, c_new, kept)
      S <- if (cplx) re2$T[seq_len(kept), seq_len(kept), drop = FALSE]
           else crossprod(Q, Ma %*% Q)
      if (!cplx && c_new && kept > c_new) {
        S[(c_new + 1L):kept, seq_len(c_new)] <- 0
      }
      Va <- V[, act, drop = FALSE] %*% Q
      if (l) {
        Mtop <- cbind(M[seq_len(l), seq_len(l), drop = FALSE],
                      M[seq_len(l), act, drop = FALSE] %*% Q)
        Mbot <- cbind(matrix(0, kept, l), S)
        M <- rbind(Mtop, Mbot)
        V <- cbind(V[, seq_len(l), drop = FALSE], Va)
      } else {
        M <- S
        V <- Va
      }
      brow <- c(numeric(l), beta * Q[ma, ])
      if (c_new) brow[l + seq_len(c_new)] <- 0
    }
    l <- l + c_new

    if (is.null(ex$vnext)) {
      ## An invariant subspace: every Ritz pair of it is exact and there is no
      ## direction to continue from, so the next round starts from a random one
      ## outside it with nothing coupling it to the retained block.
      v <- eigen_random_vector(n, Op$dtype, seed, rounds)
      if (!is.null(V)) v <- orth_against(V, v)
      nv <- col_norms(v)
      if (nv <= eps * sqrt(n)) {
        stop_reason <- "no direction left outside the retained subspace"
        break
      }
      v <- v / nv
      if (!is.null(brow)) brow[] <- 0
    } else {
      v <- ex$vnext
    }
  }

  values <- lock_value
  vectors <- lock_vec
  residuals <- lock_res
  thetas <- lock_theta
  short <- k - length(values)
  if (short > 0L && length(best_value)) {
    take <- seq_len(min(short, length(best_value)))
    values <- c(values, best_value[take])
    vectors <- cbind(vectors, best_vec[, take, drop = FALSE])
    residuals <- c(residuals, best_res[take])
    thetas <- c(thetas, best_theta[take])
  }
  if (!length(values)) {
    stopf("eigs() produced no pair at all; the operator applied to the starting vector is zero")
  }
  if (length(values) > k) {
    keep <- ritz_order(thetas, sel)[seq_len(k)]
    values <- values[keep]; vectors <- vectors[, keep, drop = FALSE]
    residuals <- residuals[keep]; thetas <- thetas[keep]
  }
  ## A spectrum with no non-real member comes back in real storage, as base
  ## eigen() returns it.
  if (all(Im(values) == 0) && all(Im(vectors) == 0)) {
    values <- Re(values)
    thetas <- Re(thetas)
    if (!cplx) vectors <- Re(vectors)
  }

  ord <- ritz_order(thetas, sel)
  list(values = values[ord], vectors = vectors[, ord, drop = FALSE],
       residuals = residuals[ord], steps = steps, rounds = rounds,
       anorm_lb = anorm_lb, stop_reason = stop_reason)
}

## One expansion of the Krylov decomposition from size p to at most m, spending
## at most max_applies applies. The unit vector v is the direction b couples the
## retained block to; its row of M is b, and each new column of M holds the
## coefficients the orthogonalisation removed, measured, plus the norm of what
## was left on the subdiagonal.
arnoldi_expand <- function(Op, V, M, brow, v, m, max_applies, cplx) {
  n <- nrow(v)
  eps <- .Machine$double.eps
  p <- if (is.null(V)) 0L else ncol(V)
  zero <- if (cplx) 0 + 0i else 0
  Vn <- matrix(zero, n, m)
  Mn <- matrix(zero, m, m)
  if (p) {
    Vn[, seq_len(p)] <- V
    Mn[seq_len(p), seq_len(p)] <- M
    Mn[p + 1L, seq_len(p)] <- brow
  }
  Vn[, p + 1L] <- v

  j <- p
  applies <- 0L
  beta <- 0
  vnext <- NULL
  apply_norm <- 0
  repeat {
    j <- j + 1L
    applies <- applies + 1L
    W <- linop_apply(Op, Vn[, j, drop = FALSE], "N")
    ## The breakdown test is relative to the size of the direction before the
    ## projection, as in lanczos_run() and GMRES.
    w0 <- col_norms(W)
    apply_norm <- max(apply_norm, w0)
    h <- orth_project(Vn[, seq_len(j), drop = FALSE], W)
    W <- h$W
    Mn[seq_len(j), j] <- h$h[, 1L]
    beta <- col_norms(W)
    if (beta <= eps * w0) { beta <- 0; break }
    if (j >= m || applies >= max_applies) { vnext <- W / beta; break }
    Vn[, j + 1L] <- W / beta
    Mn[j + 1L, j] <- beta
  }
  list(V = Vn[, seq_len(j), drop = FALSE], M = Mn[seq_len(j), seq_len(j), drop = FALSE],
       beta = beta, vnext = vnext, applies = applies, apply_norm = apply_norm)
}

## Which diagonal entries of the complex Schur form of a real matrix are real
## eigenvalues. LAPACK's real driver returns real eigenvalues exactly real and
## non-real ones in exact conjugate pairs, so its count of real ones is exact;
## the Schur entries with the smallest imaginary parts are those. Everything is
## non-real-capable for a complex operator, where conjugation pairs nothing.
schur_real_flags <- function(Ma, lam, cplx) {
  if (cplx) return(rep(FALSE, length(lam)))
  ev <- eigen(Ma, only.values = TRUE)$values
  nreal <- if (is.complex(ev)) sum(Im(ev) == 0) else length(ev)
  rank(Mod(Im(lam)), ties.method = "first") <= nreal
}

## The conjugate partner of each non-real eigenvalue of a real matrix, NA for a
## real one or for a complex operator. Paired greedily by nearest conjugate.
conjugate_partners <- function(lam, real_ev) {
  partner <- rep(NA_integer_, length(lam))
  open <- which(!real_ev)
  open <- open[order(-Im(lam[open]))]
  while (length(open) >= 2L) {
    i <- open[1L]
    others <- open[-1L]
    jj <- others[which.min(Mod(lam[others] - Conj(lam[i])))]
    partner[i] <- jj
    partner[jj] <- i
    open <- setdiff(open, c(i, jj))
  }
  partner
}

## A cut of the ordering ord after `cut` entries that separates no conjugate
## pair. grow = TRUE brings a missing partner forward to just after the cut and
## includes it; grow = FALSE drops the member whose partner falls outside.
closed_cut <- function(ord, cut, partner, grow) {
  repeat {
    if (cut <= 0L) return(list(ord = ord, cut = 0L))
    inside <- ord[seq_len(cut)]
    pi <- partner[inside]
    miss <- which(!is.na(pi) & !pi %in% inside)
    if (!length(miss)) return(list(ord = ord, cut = cut))
    if (grow && cut < length(ord)) {
      mate <- pi[miss[1L]]
      ord <- c(inside, mate, setdiff(ord[-seq_len(cut)], mate))
      cut <- cut + 1L
    } else {
      drop <- inside[miss[length(miss)]]
      ord <- c(setdiff(inside, drop), drop, ord[-seq_len(cut)])
      cut <- cut - 1L
    }
  }
}

## A real orthonormal basis of the span of the leading `kept` columns of the
## complex Schur basis U, nested so that its first `c` columns span the leading
## c. Both spans are closed under conjugation, so each has a real basis of the
## same dimension, read off the left singular vectors of the real and imaginary
## parts side by side.
real_nested_basis <- function(U, c, kept) {
  parts <- function(cols) cbind(Re(U[, cols, drop = FALSE]), Im(U[, cols, drop = FALSE]))
  B1 <- if (c) svd(parts(seq_len(c)), nu = c, nv = 0L)$u else NULL
  if (kept == c) return(B1)
  W <- parts(seq_len(kept))
  if (c) {
    W <- W - B1 %*% crossprod(B1, W)
    W <- W - B1 %*% crossprod(B1, W)
  }
  B2 <- svd(W, nu = kept - c, nv = 0L)$u
  cbind(B1, B2)
}
