## Phase 6's first function space: operators on functions on the circle, held as
## operators on their Fourier coefficients, a sequence on Z.
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
#' gives one condition `u^(d)(+-1) = value` per order. [`solve()`][solve.linop] takes the
#' Chebyshev coefficients of the right-hand side and returns those of the
#' solution, and [eigs()] with `which = "smallest"` the eigenvalues of smallest
#' magnitude; both choose the number of coefficients and certify the exact
#' residual on the untruncated system. The boundary conditions define the
#' operator's domain: [adjoint()] returns the adjoint operator with the adjoint
#' domain, by integration by parts (mixed conditions where the boundary form
#' requires them), and the operator records whether it is self-adjoint in L2.
#' For a self-adjoint `-(k u')' + q u` with `k > 0` and homogeneous Dirichlet,
#' Neumann or Robin conditions, `eigs(..., method = "galerkin")` discretises the
#' weak form (natural conditions enter it as boundary terms, Dirichlet ones as
#' a basis that vanishes there) and returns upper bounds of the eigenvalues
#' (Courant-Fischer), lifted back to Chebyshev coefficients.
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
#'   The default, `c(0, 0, -1)`, is `-d^2/dx^2`. For `"chebyshev"` it may also be
#'   a list of polynomial coefficients `a_0(x), ..., a_p(x)`, each given by its
#'   Chebyshev coefficients, for `sum a_m(x) u^(m)`.
#'   Each entry of that list may also be a function of `x`.
#' @param multiplier The coefficients of the multiplying function: Fourier
#'   coefficients from frequency `-d` to `d` (an odd number of them) on the
#'   circle, Chebyshev coefficients `a_0, a_1, ...` on the interval, monomial
#'   coefficients `v_0, v_1, ...` of `V(x)` on the line and half line. A single
#'   number is a constant.
#'   A function of `x` is sampled adaptively in the selected space. It must
#'   return one finite value per point. Scalar functions that error on vector
#'   input are evaluated pointwise after one detected fallback.
#'   Certificate construction may make one cached complex vector call to estimate
#'   an analytic tail on an ellipse or strip. The sampled maximum and analyticity
#'   are estimates. Real-only functions retain coefficient and off-grid estimates.
#' @param boundary For `"chebyshev"`, a list with `at` (each `-1` or `1`),
#'   `value`, and either `derivative` (one per condition) or `weights` (a list of
#'   weight vectors over `u, u', ...`, for mixed conditions), one condition per
#'   order of the operator. A term of a sum may omit it; the problem needs it.
#'   `value` may be a function of the endpoint, or a list with a function or
#'   scalar per endpoint. Functions are approximated on the interval first.
#' @param tol Relative tolerance for adaptive function coefficients. Defaults
#'   to machine epsilon. Chopping is an estimate recorded in the certificate.
#'   Solve and eigen requests may shorten stored multipliers when the estimated
#'   sup-norm error fits one tenth of the request and the arithmetic floor permits it.
#' @param nmax Maximum sample grid length for each function. An unresolved
#'   function is refused by argument name. Defaults to `65536`.
#' @references Aurentz, J. L. and Trefethen, L. N. (2017). Chopping a Chebyshev
#'   series. ACM Transactions on Mathematical Software 43(4).
#'   \doi{10.1145/2998442}.
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
                        boundary = NULL, tol = .Machine$double.eps, nmax = 65536L) {
  records <- list()
  convert <- function(x, name) {
    if (!is.function(x)) return(x)
    z <- function_coefficients(x, kind, name, tol, nmax)
    records[[length(records) + 1L]] <<- z
    z$coefficients
  }
  if (is_scalar_string(kind) && kind %in% names(coefficient_registry())) {
    function_multiplier <- is.function(multiplier)
    multiplier <- convert(multiplier, "multiplier")
    if (kind == "chebyshev") {
      if (is.function(derivative)) derivative <- list(derivative)
      if (is.list(derivative)) derivative <- lapply(seq_along(derivative), function(i)
        convert(derivative[[i]], sprintf("derivative[[%d]]", i)))
      if (is.list(boundary) && !is.null(boundary$value)) {
        bv <- boundary$value
        if (is.function(bv)) {
          cf <- convert(bv, "boundary$value")
          boundary$value <- vapply(boundary$at, function(x) cheb_value_at(cf, x), numeric(1))
        } else if (is.list(bv)) {
          boundary$value <- vapply(seq_along(bv), function(i) {
            cf <- convert(bv[[i]], sprintf("boundary$value[[%d]]", i))
            if (is.function(bv[[i]])) cheb_value_at(cf, boundary$at[i]) else cf
          }, numeric(1))
        }
      }
    }
    if (kind %in% c("hermite", "laguerre") && function_multiplier) {
      if (any(Im(multiplier) != 0)) stopf("multiplier function must be real in %s", kind)
      A <- new_basis_operator(kind, derivative, Re(multiplier), "orthogonal")
    } else A <- new_space_coefficients(kind, derivative, multiplier, boundary)
  } else A <- new_space_coefficients(kind, derivative, multiplier, boundary)
  carry_function_records(A, records)
}

new_space_coefficients <- function(kind, derivative, multiplier, boundary) {
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

fourier_section_norm <- function(A) max(Mod(A$args$diag)) + A$args$A$args$vnorm

fourier_pair_converged <- function(A, rn, theta, tol, floor_const) {
  target <- tol * fourier_scale(A$args$A, theta)
  floor <- floor_const * .Machine$double.eps * fourier_section_norm(A)
  rn <= pmax(target / 2 - floor, floor)
}

fourier_certificate <- function(A, values, vectors, tol, iterations, maxit,
                                floor_const = SOLVE_FLOOR_CONST,
                                requested = NULL, stop_reason = NULL) {
  P <- A$args$A
  band_certificate(A, values, vectors, tol, iterations, maxit, floor_const,
                   requested, stop_reason, P,
                   hn = fourier_section_norm(A),
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
                      function(op) op$cost, certify = fourier_certify, overwrite = TRUE,
                      converged = fourier_pair_converged)
}
