## The complex Schur form of a small dense matrix, M = U T U^H with U unitary and
## T upper triangular, the reordering that brings chosen eigenvalues to its
## leading block, and the eigenvectors of a triangular T.
##
## The Krylov-Schur restart needs an orthonormal basis of the invariant subspace
## a chosen set of Ritz values belongs to, and the leading columns of a reordered
## Schur basis are exactly that. Base R reaches LAPACK's eigenvalue drivers
## through eigen() but not its Schur routines (zgees, ztrexc), so the three are
## written here, in the form LAPACK gives them: Householder reduction to
## Hessenberg form, the shifted QR iteration with Wilkinson shifts and an
## exceptional shift every tenth sweep (Golub and Van Loan, Matrix Computations,
## 4th ed., algorithms 7.4.2 and 7.5.1; zlahqr), and a swap of adjacent diagonal
## entries by one rotation (ztrexc). The matrices are the projected ones of an
## eigensolver, of order ncv, so the cost is that of the order rather than of n.
##
## Eigenvectors of a matrix whose eigenvectors are nearly parallel are themselves
## ill conditioned, while the invariant subspace of a separated group of
## eigenvalues is not. That is why the restart uses the Schur basis rather than
## orthonormalised eigenvectors, and why eigenvectors are formed only for the pairs
## that are measured.

## A rotation G = [[c, s], [-conj(s), c]], c real, with G (a, b)^T = (r, 0)^T.
givens_complex <- function(a, b) {
  if (b == 0) return(c(1, 0))
  if (a == 0) return(c(0, Conj(b) / Mod(b)))
  r <- sqrt(Mod(a)^2 + Mod(b)^2)
  c(Mod(a) / r, (a / Mod(a)) * Conj(b) / r)
}

## The eigenvalue of the trailing 2 x 2 block [[a, b], [c, d]] closer to d.
wilkinson_shift <- function(a, b, c, d) {
  half <- (a - d) / 2
  disc <- sqrt(half * half + b * c)
  l1 <- (a + d) / 2 + disc
  l2 <- (a + d) / 2 - disc
  if (Mod(l1 - d) <= Mod(l2 - d)) l1 else l2
}

## M = U T U^H. M may be real or complex; T and U are complex.
schur_complex <- function(M) {
  m <- nrow(M)
  H <- M
  storage.mode(H) <- "complex"
  U <- diag(1 + 0i, m)
  if (m == 1L) return(list(T = H, U = U))
  eps <- .Machine$double.eps
  scale <- max(Mod(H))
  if (scale == 0) return(list(T = H, U = U))

  ## Householder reduction: P = I - 2 v v^H maps the column below the
  ## subdiagonal onto its first entry, and H <- P H P, U <- U P.
  for (j in seq_len(m - 2L)) {
    rows <- (j + 1L):m
    x <- H[rows, j]
    if (all(x[-1L] == 0)) next
    xn <- sqrt(sum(Mod(x)^2))
    ph <- if (x[1L] != 0) x[1L] / Mod(x[1L]) else 1 + 0i
    v <- x
    v[1L] <- v[1L] + ph * xn
    v <- v / sqrt(sum(Mod(v)^2))
    vh <- Conj(v)
    H[rows, ] <- H[rows, , drop = FALSE] - 2 * outer(v, as.vector(vh %*% H[rows, , drop = FALSE]))
    H[, rows] <- H[, rows, drop = FALSE] - 2 * outer(as.vector(H[, rows, drop = FALSE] %*% v), vh)
    U[, rows] <- U[, rows, drop = FALSE] - 2 * outer(as.vector(U[, rows, drop = FALSE] %*% v), vh)
    H[(j + 2L):m, j] <- 0
  }

  ## The shifted QR iteration on the active window lo..hi. A subdiagonal entry
  ## below eps times its two diagonal neighbours is set to zero, which splits the
  ## problem; the window then closes from the bottom as eigenvalues converge.
  hi <- m
  iter <- 0L
  total <- 0L
  limit <- 30L * m
  while (hi > 1L) {
    lo <- hi
    while (lo > 1L) {
      s <- Mod(H[lo - 1L, lo - 1L]) + Mod(H[lo, lo])
      if (s == 0) s <- scale
      if (Mod(H[lo, lo - 1L]) <= eps * s) {
        H[lo, lo - 1L] <- 0
        break
      }
      lo <- lo - 1L
    }
    if (lo == hi) {
      hi <- hi - 1L
      iter <- 0L
      next
    }
    iter <- iter + 1L
    total <- total + 1L
    if (total > limit) {
      stopf("the Schur iteration did not converge in %d sweeps on a %d x %d matrix", limit, m, m)
    }
    mu <- if (iter %% 10L == 0L) H[hi, hi] + 0.75 * Mod(H[hi, hi - 1L])
          else wilkinson_shift(H[hi - 1L, hi - 1L], H[hi - 1L, hi],
                               H[hi, hi - 1L], H[hi, hi])

    ## One explicit shifted QR step on the window: H_w - mu I = Q R, then
    ## R Q + mu I. Rows lo..hi take Q^H across every column to their right and
    ## columns lo..hi take Q down every row above, so the whole of H stays
    ## similar to M rather than only the window.
    win <- lo:hi
    dg <- cbind(win, win)
    H[dg] <- H[dg] - mu
    nrot <- hi - lo
    rot <- vector("list", nrot)
    for (r in seq_len(nrot)) {
      i <- lo + r - 1L
      g <- givens_complex(H[i, i], H[i + 1L, i])
      rot[[r]] <- g
      cols <- i:m
      top <- H[i, cols]
      bot <- H[i + 1L, cols]
      H[i, cols] <- g[1L] * top + g[2L] * bot
      H[i + 1L, cols] <- -Conj(g[2L]) * top + g[1L] * bot
    }
    for (r in seq_len(nrot)) {
      i <- lo + r - 1L
      g <- rot[[r]]
      rows <- seq_len(i + 1L)
      a <- H[rows, i]
      b <- H[rows, i + 1L]
      H[rows, i] <- g[1L] * a + Conj(g[2L]) * b
      H[rows, i + 1L] <- -g[2L] * a + g[1L] * b
      a <- U[, i]
      b <- U[, i + 1L]
      U[, i] <- g[1L] * a + Conj(g[2L]) * b
      U[, i + 1L] <- -g[2L] * a + g[1L] * b
    }
    H[dg] <- H[dg] + mu
  }
  H[lower.tri(H)] <- 0
  list(T = H, U = U)
}

## Swap the diagonal entries k and k + 1 of an upper triangular T by one
## rotation, carrying U with it. The first column of the rotation is the
## eigenvector of the 2 x 2 block for its second eigenvalue, (T[k, k+1],
## T[k+1, k+1] - T[k, k]), so after it that eigenvalue leads.
schur_swap <- function(T, U, k) {
  m <- nrow(T)
  k1 <- k + 1L
  a <- T[k, k]
  d <- T[k1, k1]
  x1 <- T[k, k1]
  x2 <- d - a
  nx <- sqrt(Mod(x1)^2 + Mod(x2)^2)
  if (nx == 0) return(list(T = T, U = U))
  x1 <- x1 / nx
  x2 <- x2 / nx
  cols <- k:m
  top <- T[k, cols]
  bot <- T[k1, cols]
  T[k, cols] <- Conj(x1) * top + Conj(x2) * bot
  T[k1, cols] <- -x2 * top + x1 * bot
  rows <- seq_len(k1)
  a1 <- T[rows, k]
  b1 <- T[rows, k1]
  T[rows, k] <- x1 * a1 + x2 * b1
  T[rows, k1] <- -Conj(x2) * a1 + Conj(x1) * b1
  a1 <- U[, k]
  b1 <- U[, k1]
  U[, k] <- x1 * a1 + x2 * b1
  U[, k1] <- -Conj(x2) * a1 + Conj(x1) * b1
  T[k1, k] <- 0
  T[k, k] <- d
  T[k1, k1] <- a
  list(T = T, U = U)
}

## Bring the eigenvalues at positions front[1], front[2], ... of the diagonal to
## positions 1, 2, ..., in that order, by adjacent swaps. pos records which
## original position sits where afterwards, so a caller can carry per-eigenvalue
## information through the permutation.
schur_reorder <- function(T, U, front) {
  pos <- seq_len(nrow(T))
  for (t in seq_along(front)) {
    j <- match(front[t], pos)
    while (j > t) {
      sw <- schur_swap(T, U, j - 1L)
      T <- sw$T
      U <- sw$U
      pos[c(j - 1L, j)] <- pos[c(j, j - 1L)]
      j <- j - 1L
    }
  }
  list(T = T, U = U, pos = pos)
}

## x solving (A - theta I) x = b, where a singular value of A - theta I below
## smin is raised to smin. That is ztrevc's treatment of a small pivot, written
## for a general matrix: theta repeated in A makes the system singular, and for
## a semisimple repetition any solution gives an eigenvector, which the floor
## picks one of.
floored_shift_solve <- function(A, theta, b, smin) {
  S <- A
  diag(S) <- diag(S) - theta
  sv <- svd(S)
  d <- pmax(sv$d, smin)
  sv$v %*% ((cross_adjoint(sv$u, b)) / d)
}

## The unit eigenvector of upper triangular T for its i-th diagonal entry.
schur_eigvec <- function(T, i, smin) {
  y <- complex(nrow(T))
  y[i] <- 1
  if (i > 1L) {
    idx <- seq_len(i - 1L)
    y[idx] <- floored_shift_solve(T[idx, idx, drop = FALSE], T[i, i],
                                  -T[idx, i, drop = FALSE], smin)
  }
  y / sqrt(sum(Mod(y)^2))
}
