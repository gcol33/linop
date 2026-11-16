## With C = [H - tau I; beta e_m^H], harmonic extraction solves
## C^H C z = (theta - tau) (H - tau I)^H z (SLEPc STR-9, eq. 26).
## C^+ [I; 0] has eigenvalues 1 / (theta - tau). Its SVD construction
## avoids normal equations and admits a singular H - tau I. At an invariant
## subspace (beta = 0), every eigenpair of H is already an eigenpair of A.
harmonic_projection <- function(H, beta, target) {
  if (beta == 0) return(list(M = H, reciprocal = FALSE))
  m <- nrow(H)
  C <- H
  diag(C) <- diag(C) - target
  C <- rbind(C, c(rep(0, m - 1L), beta))
  s <- svd(C)
  M <- s$v %*% (cross_adjoint(s$u[seq_len(m), , drop = FALSE], diag(m)) / s$d)
  list(M = M, reciprocal = TRUE)
}

harmonic_values <- function(values, target, reciprocal) {
  if (!reciprocal) return(values)
  theta <- rep(if (is.complex(values)) as.complex(Inf) else Inf, length(values))
  finite <- Mod(values) > 0
  theta[finite] <- target + 1 / values[finite]
  theta
}

## The leading harmonic Schur space is invariant for the reciprocal projected
## matrix. Its residual in [V, vnext] has rank one. Factoring that small
## residual recovers the translated direction of STR-9, eqs. 46 and 49,
## including when H - tau I is singular, without a solve with that matrix.
harmonic_restart <- function(V, H, beta, vnext, Q) {
  p <- ncol(Q)
  S <- cross_adjoint(Q, H %*% Q)
  D <- rbind(H %*% Q - Q %*% S, beta * Q[nrow(Q), ])
  d <- svd(D, nu = 1L, nv = 1L)
  Y <- V %*% Q
  w <- V %*% d$u[seq_len(ncol(V)), , drop = FALSE] +
    vnext * d$u[ncol(V) + 1L, 1L]
  b <- d$d[1L] * Conj(d$v[, 1L])
  h <- orth_project(Y, w)
  S <- S + h$h %*% matrix(b, 1L, p)
  nw <- col_norms(h$W)
  list(Y = Y, H = S, v = h$W / nw, brow = nw * b)
}

## Schur vectors span the retained harmonic space even when the individual
## eigenvectors are poorly conditioned. Real and imaginary parts of a
## conjugation-closed selection give a real orthonormal basis.
harmonic_keep_basis <- function(M, target, p, cplx, reciprocal, locked = 0L) {
  sch <- schur_complex(M)
  lam <- diag(sch$T)
  partner <- if (cplx) rep(NA_integer_, length(lam))
             else conjugate_partners(lam, schur_real_flags(M, lam, cplx))
  theta <- harmonic_values(lam, target, reciprocal)
  cut <- closed_cut(ritz_order(theta, "LM", target), locked + p, partner, grow = FALSE)
  re <- schur_reorder(sch$T, sch$U, cut$ord[seq_len(cut$cut)])
  Q <- if (cplx) re$U[, seq_len(cut$cut), drop = FALSE]
       else real_nested_basis(re$U, locked, cut$cut)
  Q[, locked + seq_len(cut$cut - locked), drop = FALSE]
}
