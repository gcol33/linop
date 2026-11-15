## Phase 6's first function space: operators on functions on the circle, held as
## operators on their Fourier coefficients, a sequence on Z.
## dev_notes/phase6-function-spaces-design.md has the design.
##
##     (H u)_k = s(k) u_k + sum_{|j| <= d} vhat_j u_{k-j}
##
## is a constant-coefficient differential operator sum_m c_m (d/dx)^m, whose
## symbol s(k) = sum_m c_m (i k)^m is diagonal in the Fourier basis, plus
## multiplication by a trigonometric polynomial V(x) = sum_j vhat_j e^{ijx} of
## degree d, which is a banded Laurent operator. So H is the Jacobi node of step
## 3 with a band of d instead of 1 and an unbounded diagonal, and the
## unboundedness is the new mathematics: for an operator of order p >= 1 with a
## real symbol and a real V, H is self-adjoint on the domain where s(k) u_k is in
## l^2 (V is bounded, so Kato-Rellich), and |s(k)| -> Inf makes its resolvent
## compact and its spectrum discrete. A residual bound is then a distance to an
## eigenvalue, and there is no essential spectrum for a finite section to pollute.

## The Fourier coefficient vhat_j of the multiplier, for any integer j.
fourier_vhat <- function(P, j) {
  d <- P$args$d
  out <- rep(if (is.complex(P$args$vhat)) 0i else 0, length(j))
  inb <- abs(j) <= d
  out[inb] <- P$args$vhat[j[inb] + d + 1L]
  out
}

fourier_symbol <- function(P, k) {
  cm <- P$args$coef
  s <- rep(0i, length(k))
  for (m in seq_along(cm)) s <- s + cm[m] * (1i * k)^(m - 1L)
  if (P$args$real_symbol) Re(s) else s
}

fourier_apply <- function(op, X, mode) {
  stopf(paste0("an operator on the Fourier coefficients of a function has no block to act on.\n",
               "  Truncate it with finite_section(A, n) and apply that."))
}

#' An operator on a function space
#'
#' The one constructor for operators on functions, held as operators on their
#' coefficients in the space's basis. `kind = "fourier"` is the circle: a
#' function is its Fourier coefficients, a sequence on the integers, and the
#' operator is a constant-coefficient differential operator plus multiplication
#' by a trigonometric polynomial,
#'
#' `H u = sum_m derivative[m + 1] (d/dx)^m u + V u`, with
#' `V(x) = sum_j multiplier[j + d + 1] exp(i j x)`.
#'
#' In coefficients that is a banded operator on `l^2(Z)` with an unbounded
#' diagonal, so it has infinite dimensions: [finite_section()] truncates it to
#' the frequencies `-n` to `n`, and [eigs()] runs on it directly and chooses the
#' truncation. With a real symbol (the derivative coefficients times `i^m` real)
#' and a real `V` (the multiplier conjugate-symmetric) it is self-adjoint by
#' construction, and of order at least one its spectrum is discrete, so the
#' certificate's residual bound is a distance to an eigenvalue of the operator on
#' the circle itself.
#'
#' `kind = "chebyshev"` is the interval `[-1, 1]`: a function is its Chebyshev
#' coefficients, the operator `sum_m c_m u^(m) + a(x) u` is held in the
#' ultraspherical method (Olver and Townsend), mapping Chebyshev coefficients to
#' those of the ultraspherical basis `C^(p)` where it is banded, and `boundary`
#' gives one condition `u^(d)(+-1) = value` per order. [solve()] takes the
#' Chebyshev coefficients of the right-hand side and returns those of the
#' solution, and [eigs()] with `which = "smallest"` the eigenvalues of smallest
#' magnitude; both choose the number of coefficients and certify the exact
#' residual on the untruncated system.
#'
#' `kind = "hermite"` is the line, in the Hermite functions, where `x` and
#' `d/dx` are tridiagonal: `sum_m c_m (d/dx)^m + V(x)` with `V` a polynomial given
#' by its monomial coefficients. `kind = "laguerre"` is the half line, in the
#' Laguerre functions `exp(-x/2) L_n(x)`, where `x` is tridiagonal and the
#' Laguerre operator `-(x u')' + x u / 4` is diagonal: `c` times that operator
#' plus a polynomial `V(x)`, with `derivative = c`. Both are banded operators on
#' coefficient sequences, truncated and certified as on the circle.
#'
#' @param kind The function space: `"fourier"`, the circle; `"chebyshev"`, the
#'   interval `[-1, 1]`; `"hermite"`, the line; `"laguerre"`, the half line.
#' @param derivative Coefficients `c_0, c_1, ..., c_p` of `sum c_m (d/dx)^m`.
#'   The default, `c(0, 0, -1)`, is `-d^2/dx^2`.
#' @param multiplier The coefficients of the multiplying function: Fourier
#'   coefficients from frequency `-d` to `d` (an odd number of them) on the
#'   circle, Chebyshev coefficients `a_0, a_1, ...` on the interval, monomial
#'   coefficients `v_0, v_1, ...` of `V(x)` on the line and half line. A single
#'   number is a constant.
#' @param boundary For `"chebyshev"`, a list with `at` (each `-1` or `1`),
#'   `derivative` and `value`, one entry per order of the operator.
#' @return A `linop` of dimension `Inf x Inf`.
#' @examples
#' ## -u'' + V u with V = a^2 sin^2 x - a cos x has the ground state
#' ## exp(a cos x) at eigenvalue 0
#' a <- 2
#' H <- linop_space("fourier", multiplier = c(-a^2 / 4, -a / 2, a^2 / 2, -a / 2, -a^2 / 4))
#' fit <- eigs(H, k = 1, which = "smallest_algebraic")
#' fit$values
#' fit$n
#'
#' ## u'' = -pi^2 sin(pi x) on [-1, 1] with u(-1) = u(1) = 0, whose solution is
#' ## sin(pi x); the right-hand side's Chebyshev coefficients are Bessel values
#' k <- 0:40
#' f <- -pi^2 * ifelse(k %% 2 == 1, 2 * (-1)^((k - 1) / 2) * besselJ(pi, k), 0)
#' L <- linop_space("chebyshev", derivative = c(0, 0, 1),
#'                  boundary = list(at = c(-1, 1), derivative = c(0, 0), value = c(0, 0)))
#' u <- solve(L, f)
#' x <- 0.3
#' sum(u * cos((seq_along(u) - 1) * acos(x))) - sin(pi * x)
#' @export
linop_space <- function(kind = "fourier", derivative = c(0, 0, -1), multiplier = 0,
                        boundary = NULL) {
  kinds <- c("fourier", "chebyshev", "hermite", "laguerre")
  if (!is_scalar_string(kind) || !kind %in% kinds) {
    stopf("kind is one of %s; got %s", paste0("\"", kinds, "\"", collapse = ", "), deparse(kind))
  }
  if (kind == "chebyshev") return(new_chebyshev(derivative, multiplier, boundary))
  if (kind %in% c("hermite", "laguerre")) {
    if (!is.numeric(multiplier) || !length(multiplier) || anyNA(multiplier)) {
      stopf("multiplier is the coefficients v_0, v_1, ... of the polynomial V(x) = sum v_m x^m")
    }
    if (!is.numeric(derivative) || !length(derivative) || anyNA(derivative)) {
      stopf("derivative must be real coefficients")
    }
    return(new_basis_operator(kind, derivative, multiplier))
  }
  if (!(is.numeric(derivative) || is.complex(derivative)) || !length(derivative) ||
      anyNA(derivative)) {
    stopf("derivative must be the coefficients c_0, ..., c_p, at least one of them")
  }
  if (!(is.numeric(multiplier) || is.complex(multiplier)) || length(multiplier) %% 2L != 1L ||
      anyNA(multiplier)) {
    stopf("multiplier must be the Fourier coefficients from -d to d, an odd number of them")
  }
  coef <- derivative
  while (length(coef) > 1L && coef[length(coef)] == 0) coef <- coef[-length(coef)]
  order <- length(coef) - 1L
  d <- (length(multiplier) - 1L) %/% 2L
  ## i^m c_m real for every m is exactly a real symbol on every integer k.
  real_symbol <- all(Mod(Im(coef * 1i^(seq_along(coef) - 1L))) == 0)
  real_v <- all(multiplier == Conj(rev(multiplier)))
  self_adjoint <- real_symbol && real_v
  ## The coefficient matrix is real exactly when V is even as well as real.
  real_matrix <- real_symbol && all(Im(multiplier) == 0) && all(multiplier == rev(multiplier))
  caps <- list()
  if (self_adjoint) caps$hermitian <- capability(TRUE, ev_construction())
  if (real_matrix && self_adjoint) caps$symmetric <- capability(TRUE, ev_construction())
  new_linop("fourier", c(Inf, Inf), if (real_matrix) "double" else "complex",
            do.call(new_caps, caps),
            list(coef = coef, order = order, vhat = multiplier, d = d,
                 real_symbol = real_symbol, self_adjoint = self_adjoint,
                 discrete = self_adjoint && order >= 1L,
                 vnorm = sum(Mod(multiplier))),
            cost = Inf)
}

## ------------------------------------------------------------- the section --

fourier_section <- function(P, n) {
  if (n < max(1L, P$args$d)) {
    stopf("n = %d is narrower than the multiplier's degree %d; the section has to hold a whole band",
          n, P$args$d)
  }
  k <- seq.int(-n, n)
  m <- 2L * n + 1L
  diagv <- fourier_symbol(P, k) + fourier_vhat(P, 0L)
  real <- P$dtype == "double"
  caps <- list()
  if (P$args$self_adjoint) caps$hermitian <- capability(TRUE, ev_construction())
  if (real && P$args$self_adjoint) caps$symmetric <- capability(TRUE, ev_construction())
  new_linop("fourier_section", c(m, m), P$dtype, do.call(new_caps, caps),
            list(A = P, n = n, diag = if (real) Re(diagv) else diagv),
            cost = (2 * P$args$d + 1) * m, finite = TRUE)
}

## (H_n u)_k = diag_k u_k + sum_{j != 0} vhat_j u_{k-j} on |k| <= n. Its transpose
## has vhat_{-j} at offset j, its adjoint conj(vhat_{-j}), its conjugate
## conj(vhat_j).
fourier_section_apply <- function(op, X, mode) {
  P <- op$args$A
  d <- P$args$d
  m <- nrow(X)
  dg <- op$args$diag
  if (mode %in% c("C", "R") && is.complex(dg)) dg <- Conj(dg)
  Y <- dg * X
  for (j in setdiff(seq.int(-d, d), 0L)) {
    v <- fourier_vhat(P, if (mode %in% c("T", "C")) -j else j)
    if (mode %in% c("C", "R")) v <- Conj(v)
    if (v == 0 || abs(j) >= m) next
    if (j > 0) {
      Y[(j + 1L):m, ] <- Y[(j + 1L):m, ] + v * X[1:(m - j), , drop = FALSE]
    } else {
      a <- -j
      Y[1:(m - a), ] <- Y[1:(m - a), ] + v * X[(a + 1L):m, , drop = FALSE]
    }
  }
  Y
}

fourier_section_materialize <- function(op) {
  m <- op$dim[1L]
  E <- diag(1, m)
  if (op$dtype == "complex") storage.mode(E) <- "complex"
  fourier_section_apply(op, E, "N")
}

## ------------------------------------------------------- the certificate ----

## The section identity of step 3 with a band of d. With u solving the section
## and u~ its extension by zero, (H - q) u~ agrees with (H_n - q) u on |k| <= n,
## and outside it is nonzero only on the d frequencies each side the band reaches:
##
##     r_{n+t} = sum_{s=t}^{d} vhat_s u_{n+t-s},  r_{-(n+t)} = sum_{s=t}^{d} vhat_{-s} u_{-(n+t)+s},
##
## for t = 1..d, exactly. So ||(H - q) u~||^2 = ||(H_n - q) u||^2 + ||r_out||^2,
## both terms reported, and for a self-adjoint H
##
##     dist(q, sigma(H)) <= sqrt(residual^2 + truncation^2),
##
## a distance to an eigenvalue when the spectrum is discrete. There is no ||H|| to
## scale a tolerance by, so a pair is held to tol times its own scale,
## |q| + ||V|| + |c_p|: the eigenvalue, a bound on the bounded part, and the
## leading coefficient, which sets the spacing of the spectrum near zero.
fourier_edge <- function(P, X) {
  d <- P$args$d
  m <- nrow(X)
  n <- (m - 1L) %/% 2L
  if (d == 0L) return(rep(0, ncol(X)))
  idx <- function(l) l + n + 1L
  tot <- rep(0, ncol(X))
  for (t in seq_len(d)) {
    rp <- 0; rn <- 0
    for (s in t:d) {
      rp <- rp + fourier_vhat(P, s) * X[idx(n + t - s), ]
      rn <- rn + fourier_vhat(P, -s) * X[idx(-(n + t) + s), ]
    }
    tot <- tot + Mod(rp)^2 + Mod(rn)^2
  }
  sqrt(tot)
}

fourier_scale <- function(P, q) Mod(q) + P$args$vnorm + Mod(P$args$coef[length(P$args$coef)])

fourier_certificate <- function(A, values, vectors, tol, iterations, maxit,
                                floor_const = SOLVE_FLOOR_CONST,
                                requested = NULL, stop_reason = NULL) {
  P <- A$args$A
  band_certificate(A, values, vectors, tol, iterations, maxit, floor_const,
                   requested, stop_reason, P,
                   hn = max(Mod(A$args$diag)) + P$args$vnorm,
                   scale = function(q) fourier_scale(P, q),
                   edge = function(X) fourier_edge(P, X),
                   self_adjoint = P$args$self_adjoint, discrete = P$args$discrete,
                   reach = sprintf("the %d frequencies each side the band reaches", P$args$d),
                   space = "the circle")
}

fourier_certify <- function(ctx) {
  fourier_certificate(ctx$A, ctx$values, ctx$vectors, tol = ctx$tol,
                      iterations = ctx$iterations, maxit = ctx$maxit,
                      floor_const = ctx$floor_const, requested = ctx$requested,
                      stop_reason = ctx$stop_reason)
}

fourier_eigs <- function(A, k, control, args) {
  band_eigs(A, k, control, args,
            n_min = max(A$args$d, as.integer(ceiling(k / 2)) + 1L))
}

register_fourier_nodes <- function() {
  linop_register_node("fourier", fourier_apply, NULL, function(op) Inf,
                      spectrum = fourier_eigs, section = fourier_section, overwrite = TRUE)
  linop_register_node("fourier_section", fourier_section_apply, fourier_section_materialize,
                      function(op) op$cost, certify = fourier_certify, overwrite = TRUE)
}
