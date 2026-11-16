## Section 8.10's first rigorous class, and the first operator in the package
## whose index set is infinite.
##
##     (H u)_j = b_{j-1} u_{j-1} + a_j u_j + b_j u_{j+1}          on l^2(Z)
##
## with both coefficient sequences real and eventually constant: a_j = a_inf and
## b_j = b_inf outside a finite window. That restriction is what every closed form
## below rests on, and it is not a convenience.
##
## A finitely supported perturbation of the free operator is finite rank, hence
## compact, so by Weyl's theorem the essential spectrum is the free one. The free
## operator is Laurent with symbol a_inf + 2 b_inf cos(theta), so
##
##     sigma_ess(H) = [a_inf - 2|b_inf|, a_inf + 2|b_inf|]
##
## exactly, from the coefficients and with nothing estimated. Everything the
## finite section certifies about H -- where the eigenvector decays, how far a
## truncated eigenvalue can be from the true spectrum, whether a value has
## separated from the continuum at all -- comes back to that interval being known
## rather than approximated.
##
## Outside the window the eigenvalue equation is the free recurrence, so an
## eigenfunction is exactly geometric there rather than asymptotically geometric.
## That is what makes the truncation error a closed form instead of an existence
## statement, and it is why this class is the first one and not a general banded
## Toeplitz operator.

## Where the coefficients sit. b_j couples j and j+1, so the couplings incident to
## a diagonal window [first, last] run from b_{first-1} to b_{last}: one more than
## there are diagonal entries, which is the alignment linop_jacobi() documents.
##
## The limits may be periodic: a_j = a_inf[(j mod P) + 1] outside the window, and
## likewise b_j, anchored to the absolute index so that the window and the tail
## agree about which period they are in.
jacobi_a <- function(A, j) {
  lo <- A$args$first
  lim <- A$args$a_inf
  out <- lim[(j %% length(lim)) + 1L]
  inw <- j >= lo & j <= lo + length(A$args$a) - 1L
  out[inw] <- A$args$a[j[inw] - lo + 1L]
  out
}

jacobi_b <- function(A, j) {
  lo <- A$args$first - 1L
  lim <- A$args$b_inf
  out <- lim[(j %% length(lim)) + 1L]
  inw <- j >= lo & j <= lo + length(A$args$b) - 1L
  out[inw] <- A$args$b[j[inw] - lo + 1L]
  out
}

## The essential spectrum, exact for constant limits. For period-P limits it is
## the union of P bands, and by Floquet theory band k runs between the k-th
## eigenvalues of the P x P Bloch matrices at theta = 0 and theta = pi (the band
## functions are monotone on [0, pi]), so the edges are eigenvalues of two small
## matrices: computed, to rounding. One row per band.
jacobi_bloch <- function(A, phase) {
  a <- A$args$a_inf; b <- A$args$b_inf
  P <- length(a)
  H <- diag(a, P)
  for (j in seq_len(P - 1L)) {
    H[j, j + 1L] <- b[j]
    H[j + 1L, j] <- b[j]
  }
  H[P, 1L] <- H[P, 1L] + phase * b[P]
  H[1L, P] <- H[1L, P] + phase * b[P]
  H
}

jacobi_bands <- function(A) {
  a <- A$args$a_inf; b <- A$args$b_inf
  P <- length(a)
  if (P == 1L) return(matrix(a + c(-2, 2) * abs(b), 1L))
  bloch <- function(theta) {
    H <- jacobi_bloch(A, theta)
    sort(Re(eigen(H, symmetric = TRUE, only.values = TRUE)$values))
  }
  e0 <- bloch(1); epi <- bloch(-1)
  cbind(pmin(e0, epi), pmax(e0, epi))
}

## The hull of the essential spectrum, its lowest and highest point.
jacobi_band <- function(A) {
  B <- jacobi_bands(A)
  c(min(B[, 1L]), max(B[, 2L]))
}

## The signed distance from each value to the essential spectrum: positive
## outside every band, negative inside one.
jacobi_band_gap <- function(A, q) {
  B <- jacobi_bands(A)
  vapply(q, function(x) min(pmax(x - B[, 2L], B[, 1L] - x)), numeric(1))
}

## An upper bound on ||H||, closed form. H is self-adjoint, so its spectral radius
## is bounded by the largest absolute row sum, and every row is |a_j| + |b_{j-1}|
## + |b_j|. Outside the window that is the limiting row, so a finite maximum over
## the window and the limit covers every row there is.
##
## This is the one quantity in the package that wants an UPPER bound on a norm
## rather than a lower one, and it is not an exception to what norm2() promises:
## it multiplies eps to make an arithmetic floor, where overstating the norm
## widens the floor, and a floor that is too wide reports qualified where a
## sharper one reports pass. norm2() sits in a Rigal-Gaches denominator, where the
## safe direction is the other one.
jacobi_norm_bound <- function(A) {
  lo <- A$args$first - 1L
  P <- length(A$args$a_inf) * length(A$args$b_inf)
  ## the window, and one whole period of the limits beyond it on each side
  j <- seq.int(lo - P, lo + length(A$args$b) + P)
  max(abs(jacobi_a(A, j)) + abs(jacobi_b(A, j - 1L)) + abs(jacobi_b(A, j)))
}

## The operator has no block to act on, so nothing routes here through
## linop_apply(): require_finite_dim() refuses first, with the verb the caller
## used. Reaching the handler directly is the one way in, and it gets the same
## sentence.
jacobi_apply <- function(op, X, mode) {
  stopf(paste0("a Jacobi operator on l^2(Z) has no block to act on.\n",
               "  Truncate it with finite_section(A, n) and apply that."))
}

#' A Jacobi operator on the sequence space
#'
#' A self-adjoint operator on `l^2(Z)`, tridiagonal in the standard basis, given
#' by its diagonal and off-diagonal sequences. Both are real and eventually
#' constant: they take the values you supply on a finite window and their limiting
#' values everywhere else.
#'
#' The operator has infinite dimensions, so it has no matrix and nothing numeric
#' runs on it directly. [finite_section()] truncates it to an operator that does,
#' and the certificate that comes back describes the operator here rather than the
#' truncation.
#'
#' The window is `from` to `from + length(diagonal) - 1`. `offdiagonal[i]` couples
#' index `from + i - 2` to `from + i - 1`, so the couplings incident to the
#' diagonal window run from the one on its left edge to the one on its right, and
#' there are one more of them than there are diagonal entries.
#'
#' A constructor knows the operator is self-adjoint the way [linop_eye()] knows it
#' about the identity, so `hermitian` and `symmetric` are recorded with
#' `source = "construction"` rather than as declarations.
#'
#' @param diagonal The diagonal entries on the window, as a real vector. A
#'   discrete Schrodinger operator is this argument alone: it is the potential.
#' @param offdiagonal The couplings, as a real vector of length 1 for a constant
#'   off-diagonal, or `length(diagonal) + 1` for a perturbed one. Zero is not
#'   allowed as the limiting value.
#' @param from The index of `diagonal[1]`. The default centres the window on the
#'   origin.
#' @param diagonal_limit The value the diagonal takes outside the window, or one
#'   period of values: outside the window `a_j = diagonal_limit[(j mod P) + 1]`.
#'   With periodic limits the essential spectrum is a union of bands with gaps
#'   between them, and the certificate measures isolation against the nearest
#'   band.
#' @param offdiagonal_limit The value the off-diagonal takes outside the window,
#'   or one period of values, anchored the same way. Defaults to `offdiagonal`
#'   when that is a single number.
#' @return A `linop` of dimension `Inf x Inf`.
#' @examples
#' ## the potential v at the origin, whose top eigenvalue is sign(v) sqrt(v^2 + 4)
#' H <- linop_jacobi(diagonal = 1)
#' dim(H)
#' F <- finite_section(H, n = 50)
#' fit <- eigs(F, k = 1, which = "largest_algebraic")
#' fit$values - sqrt(5)
#' fit$certificate
#' @export
linop_jacobi <- function(diagonal, offdiagonal = 1, from = NULL,
                         diagonal_limit = 0, offdiagonal_limit = NULL) {
  if (is.complex(diagonal) || is.complex(offdiagonal)) {
    stopf(paste0("the coefficients must be real; a complex diagonal or off-diagonal makes\n",
                 "  the operator non-self-adjoint, and every bound this class certifies is a\n",
                 "  self-adjointness argument."))
  }
  if (!is.numeric(diagonal) || !length(diagonal) || anyNA(diagonal)) {
    stopf("diagonal must be a real vector of at least one entry")
  }
  L <- length(diagonal)
  if (!is.numeric(offdiagonal) || anyNA(offdiagonal) ||
      !length(offdiagonal) %in% c(1L, L + 1L)) {
    stopf(paste0("offdiagonal is one number, or %d of them for a window of %d diagonal entries.\n",
                 "  b_j couples j and j+1, so the couplings touching the window are the one on\n",
                 "  each edge and the %d between them."),
          L + 1L, L, L - 1L)
  }
  offdiagonal_limit <- offdiagonal_limit %||%
    (if (length(offdiagonal) == 1L) offdiagonal else 1)
  if (!is.numeric(offdiagonal_limit) || !length(offdiagonal_limit) ||
      !all(is.finite(offdiagonal_limit)) || any(offdiagonal_limit == 0)) {
    stopf(paste0("offdiagonal_limit must be nonzero real numbers, one period of them.\n",
                 "  At zero the operator has no coupling outside the window at all, so it is a\n",
                 "  direct sum of a finite matrix and a multiplication operator rather than a\n",
                 "  Jacobi operator, and its essential spectrum is a set of points."))
  }
  if (!is.numeric(diagonal_limit) || !length(diagonal_limit) ||
      !all(is.finite(diagonal_limit))) {
    stopf("diagonal_limit must be finite real numbers, one period of them")
  }
  ## one period for both, of the length both divide
  per <- length(diagonal_limit) * length(offdiagonal_limit) /
    gcd_int(length(diagonal_limit), length(offdiagonal_limit))
  diagonal_limit <- rep_len(as.numeric(diagonal_limit), per)
  offdiagonal_limit <- rep_len(as.numeric(offdiagonal_limit), per)
  from <- as.integer(from %||% -((L - 1L) %/% 2L))
  if (is.na(from)) stopf("from must be a single whole number")
  b <- if (length(offdiagonal) == 1L) rep(offdiagonal, L + 1L) else offdiagonal

  ## How far the perturbation reaches from the origin, counting the coupling on
  ## each edge of the window. finite_section() needs a window strictly larger
  ## than this, since the closed forms all read the tail as free.
  radius <- max(abs(from - 1L), abs(from + L - 1L))

  caps <- list(hermitian = capability(TRUE, ev_construction()),
               symmetric = capability(TRUE, ev_construction()),
               real = capability(TRUE, ev_construction()))
  new_linop("jacobi", c(Inf, Inf), "double", do.call(new_caps, caps),
            list(a = as.numeric(diagonal), b = as.numeric(b), first = from,
                 a_inf = as.numeric(diagonal_limit),
                 b_inf = as.numeric(offdiagonal_limit), radius = radius),
            cost = Inf)
}

#' How fast an eigenvector decays
#'
#' Outside the window the eigenvalue equation is the free recurrence
#' `b (u[j-1] + u[j+1]) + a u[j] = lambda u[j]`, whose solutions are `z^j` for
#' the two roots of `z^2 - mu z + 1` with `mu = (lambda - a_inf) / b_inf`. The
#' roots are reciprocal, so for `|mu| > 2` exactly one lies inside the unit disc, and
#' membership in `l^2` forces it on each side. An eigenvector is therefore
#' *exactly* geometric outside the window, at the rate this returns, rather than
#' asymptotically so.
#'
#' `NA` for a value inside the essential spectrum: there both roots lie on the
#' unit circle, nothing decays, and there is no eigenvector to decay.
#'
#' @param A A `linop` from [linop_jacobi()], or a [finite_section()] of one.
#' @param lambda One or more real values.
#' @return A numeric vector in `(0, 1)`, with `NA` where `lambda` lies in the
#'   essential spectrum.
#' @examples
#' H <- linop_jacobi(diagonal = 1)
#' ## the golden-ratio conjugate, for the eigenvalue sqrt(5)
#' decay_rate(H, sqrt(5))
#' decay_rate(H, 0)
#' @export
decay_rate <- function(A, lambda) {
  P <- jacobi_of(A, "decay_rate")
  if (!is.numeric(lambda) || !length(lambda)) {
    stopf("lambda must be a real vector of at least one value")
  }
  if (length(P$args$a_inf) > 1L) return(periodic_decay_rate(P, lambda))
  ## Written as the reciprocal of the growing root rather than as the difference
  ## of two nearly equal numbers: |mu|/2 - sqrt((mu/2)^2 - 1) cancels to nothing
  ## for a well separated eigenvalue, which is the case this is most used on.
  h <- abs(lambda - P$args$a_inf) / (2 * abs(P$args$b_inf))
  s <- sqrt(pmax(h^2 - 1, 0))
  ifelse(h > 1, 1 / (h + s), NA_real_)
}

## Periodic limits: the free recurrence b_j u_{j+1} = (lambda - a_j) u_j -
## b_{j-1} u_{j-1} over one period is a 2 x 2 transfer matrix M of determinant 1,
## so its eigenvalues are reciprocal; outside the bands |tr M| > 2 and the
## smaller one, per period, is the decay. Per site it is that to the power 1/P,
## which the width search divides by; within a period the tail oscillates about
## that rate rather than following it exactly.
periodic_decay_rate <- function(P, lambda) {
  a <- P$args$a_inf; b <- P$args$b_inf
  per <- length(a)
  vapply(lambda, function(l) {
    M <- diag(2)
    for (j in seq_len(per)) {
      bprev <- b[((j - 2L) %% per) + 1L]
      Tj <- matrix(c((l - a[j]) / b[j], 1, -bprev / b[j], 0), 2L)
      M <- Tj %*% M
    }
    tr <- M[1, 1] + M[2, 2]
    if (abs(tr) <= 2) return(NA_real_)
    h <- abs(tr) / 2
    (1 / (h + sqrt(h^2 - 1)))^(1 / per)
  }, numeric(1))
}

gcd_int <- function(x, y) if (y == 0L) x else gcd_int(y, x %% y)

## The Jacobi operator behind whatever was handed in, so decay_rate() and the
## certificate read the same argument off a section as off the operator itself.
jacobi_of <- function(A, verb) {
  if (is_linop(A)) {
    if (A$node == "jacobi") return(A)
    if (A$node == "section") return(A$args$A)
  }
  stopf(paste0("%s() expects an operator from linop_jacobi(), or a finite_section() of one.\n",
               "  Both the decay rate and the essential spectrum are read off the limiting\n",
               "  coefficients, which only that class carries."), verb)
}

register_jacobi_node <- function() {
  linop_register_node("jacobi", jacobi_apply, NULL, function(op) Inf,
                      spectrum = section_eigs, section = jacobi_section,
                      methods = jacobi_gap_hooks(), overwrite = TRUE)
}
