## Phase 6 step 2: functions on [-1, 1] held as Chebyshev coefficients, and the
## ultraspherical method for differential operators on them (Olver and Townsend,
## SIAM Review 2013).
##
## u = sum_k c_k T_k. Differentiation is banded only between different spaces:
## D_m maps Chebyshev coefficients to coefficients in the ultraspherical basis
## C^(m), column k to row k - m with weight 2^(m-1) (m-1)! k, and the conversions
## S_lambda: C^(lambda) -> C^(lambda+1) have two diagonals. So
##
##     L u = sum_m c_m u^(m) + a(x) u
##
## of order p is S_{p-1} ... S_m D_m on each derivative term and the conversion of
## a Chebyshev-basis product on the multiplier, and lands in C^(p): its domain and
## range are different spaces, which is what keeps it banded. Boundary conditions
## u^(d)(+-1) = value are dense rows. Every one of these maps a finitely supported
## coefficient vector to a finitely supported one, so the residual of a truncated
## solution on the original infinite system is computable exactly: the same
## identity the Fourier and Jacobi sections certify with.

## T_k^(d)(1) = prod_{j<d} (k^2 - j^2) / (2j + 1), and T_k^(d)(-1) = (-1)^(k+d) times it.
cheb_boundary_row <- function(at, d, N) {
  k <- seq.int(0L, N - 1L)
  v <- rep(1, N)
  for (j in seq_len(d) - 1L) v <- v * (k^2 - j^2) / (2 * j + 1)
  if (at < 0) v <- v * (-1)^(k + d)
  v
}

## A condition sum_d w_d u^(d)(at) as a row on the first N coefficients. A
## condition on one derivative is a unit weight; the adjoint of a problem has
## mixed (Robin) conditions in general, which is why the weights are kept whole.
bc_row <- function(at, w, N) {
  row <- numeric(N)
  for (d in seq_along(w) - 1L) if (w[d + 1L] != 0) row <- row + w[d + 1L] * cheb_boundary_row(at, d, N)
  row
}

## D_m c, in C^(m) coefficients.
cheb_diff <- function(c, m) {
  if (m == 0L) return(c)
  N <- length(c)
  if (N <= m) return(numeric(0))
  k <- seq.int(m, N - 1L)
  2^(m - 1) * factorial(m - 1) * k * c[k + 1L]
}

## S_lambda v: C^(lambda) -> C^(lambda + 1). T_0 = C^(1)_0, T_k = (C^(1)_k -
## C^(1)_{k-2}) / 2; for lambda >= 1, C^(lambda)_k = lambda / (lambda + k) (C^(lambda+1)_k
## - C^(lambda+1)_{k-2}). Row k reads columns k and k + 2.
cheb_convert <- function(v, lambda) {
  N <- length(v)
  if (!N) return(v)
  k <- seq.int(0L, N - 1L)
  w <- if (lambda == 0L) ifelse(k == 0L, 1, 0.5) else lambda / (lambda + k)
  out <- w * v
  if (N > 2L) out[1:(N - 2L)] <- out[1:(N - 2L)] - w[3:N] * v[3:N]
  out
}

cheb_convert_chain <- function(v, from, to) {
  for (lam in seq_len(to - from) + from - 1L) v <- cheb_convert(v, lam)
  v
}

## (a u) in Chebyshev coefficients, exactly: T_i T_k = (T_{i+k} + T_{|i-k|}) / 2.
cheb_mult <- function(a, c) {
  na <- length(a); nc <- length(c)
  if (!na || !nc) return(numeric(0))
  out <- rep(if (is.complex(a) || is.complex(c)) 0i else 0, na + nc - 1L)
  for (i in seq_len(na) - 1L) {
    if (a[i + 1L] == 0) next
    for (k in seq_len(nc) - 1L) {
      h <- a[i + 1L] * c[k + 1L] / 2
      out[i + k + 1L] <- out[i + k + 1L] + h
      out[abs(i - k) + 1L] <- out[abs(i - k) + 1L] + h
    }
  }
  out
}

pad_to <- function(v, n) if (length(v) >= n) v else c(v, rep(0, n - length(v)))

## The derivative of a Chebyshev series, in Chebyshev coefficients:
## c'_{k-1} = c'_{k+1} + 2 k c_k, halved at k = 1.
cheb_dcoef <- function(cf) {
  n <- length(cf)
  if (n <= 1L) return(0)
  d <- numeric(n + 1L)
  for (k in (n - 1L):1L) d[k] <- d[k + 2L] + 2 * k * cf[k + 1L]
  d[1L] <- d[1L] / 2
  d[seq_len(n - 1L)]
}

cheb_dcoef_n <- function(cf, m) {
  for (i in seq_len(m)) cf <- cheb_dcoef(cf)
  cf
}

## The d-th derivative of a Chebyshev series at +-1.
cheb_value_at <- function(cf, at, d = 0L) sum(cf * cheb_boundary_row(at, d, length(cf)))

## L c in C^(p), exact and of finite length. A constant coefficient takes the
## banded route, D_m then the conversions; a polynomial coefficient a_m(x) takes
## the Chebyshev derivative, the Chebyshev product and the conversions from
## C^(0), which is no longer banded but still maps a finitely supported vector to
## one, so the residual on the untruncated system stays exact.
cheb_apply_L <- function(P, c) {
  p <- P$args$order
  out <- numeric(0)
  for (m in seq_along(P$args$terms) - 1L) {
    am <- P$args$terms[[m + 1L]]
    if (all(am == 0)) next
    term <- if (length(am) == 1L) am * cheb_convert_chain(cheb_diff(c, m), m, p)
            else cheb_convert_chain(cheb_mult(am, cheb_dcoef_n(c, m)), 0L, p)
    n <- max(length(out), length(term))
    out <- pad_to(out, n) + pad_to(term, n)
  }
  out
}

trim_poly <- function(v) {
  v <- as.numeric(v)
  while (length(v) > 1L && v[length(v)] == 0) v <- v[-length(v)]
  if (!length(v)) 0 else v
}

## The N x N system on the first N coefficients: the boundary rows, then the
## first N - p rows of L. Columns are L applied to unit vectors, so the band is
## read off the operator rather than restated.
cheb_system <- function(P, N) {
  p <- P$args$order
  bc <- P$args$boundary
  A <- matrix(0, N, N)
  for (i in seq_len(p)) A[i, ] <- bc_row(bc$at[i], bc$weights[[i]], N)
  for (k in seq_len(N)) {
    e <- numeric(N); e[k] <- 1
    col <- pad_to(cheb_apply_L(P, e), N - p)
    A[p + seq_len(N - p), k] <- col[seq_len(N - p)]
  }
  A
}

cheb_apply <- function(op, X, mode) {
  stopf(paste0("an operator on Chebyshev coefficients has no block to act on.\n",
               "  solve() and eigs() run on it directly and choose the number of coefficients."))
}

## derivative is constants c_0, c_1, ... or a list of Chebyshev coefficient
## vectors a_0(x), a_1(x), ..., one per order; multiplier adds to a_0.
new_chebyshev <- function(derivative, multiplier, boundary) {
  terms <- if (is.list(derivative)) lapply(derivative, trim_poly)
           else lapply(as.numeric(derivative), trim_poly)
  if (!length(terms)) terms <- list(0)
  if (length(multiplier) && any(multiplier != 0)) {
    n <- max(length(terms[[1L]]), length(multiplier))
    terms[[1L]] <- trim_poly(pad_to(terms[[1L]], n) + pad_to(as.numeric(multiplier), n))
  }
  while (length(terms) > 1L && all(terms[[length(terms)]] == 0)) terms <- terms[-length(terms)]
  p <- length(terms) - 1L
  coef <- if (all(lengths(terms) == 1L)) unlist(terms) else NULL
  ## A term of a sum may come without boundary conditions, which belong to the
  ## problem; solve() and eigs() require them, one per order.
  if (!is.null(boundary)) boundary <- normalise_boundary(boundary)
  P <- new_linop("chebyshev", c(Inf, Inf), "double", new_caps(),
                 list(terms = terms, coef = coef, order = p, boundary = boundary,
                      domain = "chebyshev", range = sprintf("ultraspherical(%d)", p)),
                 cost = Inf)
  ## self-adjointness in L2 on the interval, NA without a full set of conditions
  P$args$self_adjoint <- cheb_self_adjoint(P)
  P
}

## Conditions as given: at, value, and either derivative (one per condition) or
## weights (a list of weight vectors over u, u', u'', ...).
normalise_boundary <- function(boundary) {
  if (!is.list(boundary) || !all(c("at", "value") %in% names(boundary)) ||
      !any(c("derivative", "weights") %in% names(boundary))) {
    stopf("boundary is a list of at (each -1 or 1), value, and derivative or weights")
  }
  nb <- length(boundary$at)
  weights <- boundary$weights
  if (is.null(weights)) {
    d <- boundary$derivative
    if (length(d) != nb || any(d < 0) || any(d != trunc(d))) {
      stopf("boundary derivatives are whole numbers from 0, one per condition")
    }
    weights <- lapply(d, function(k) c(numeric(k), 1))
  }
  if (!is.list(weights) || length(weights) != nb || length(boundary$value) != nb) {
    stopf("boundary needs at, value and derivative (or weights) of the same length")
  }
  if (!all(boundary$at %in% c(-1, 1))) stopf("boundary conditions sit at -1 or 1")
  list(at = as.numeric(boundary$at), weights = lapply(weights, as.numeric),
       value = as.numeric(boundary$value))
}

## What solve() and eigs() need before a system exists: a differential operator
## with one condition per order, each on a derivative below the order.
cheb_check_problem <- function(P, verb) {
  p <- P$args$order
  bc <- P$args$boundary
  if (p < 1L) stopf("%s() needs a differential operator; this one is of order 0", verb)
  nb <- if (is.null(bc)) 0L else length(bc$at)
  if (nb != p) {
    stopf("%s() on an operator of order %d needs %d boundary conditions; it has %d",
          verb, p, p, nb)
  }
  if (any(lengths(bc$weights) > p)) {
    stopf("boundary derivatives run from 0 to %d for an operator of order %d", p - 1L, p)
  }
  invisible(TRUE)
}

## ------------------------------------------------------------------ solve --

CHEB_N_START <- 32L
CHEB_N_MAX <- 4096L

## The computed coefficients c on N terms solve the original infinite system
## exactly with the right-hand side moved by r: the boundary residual, the rows of
## L c - F the section kept (rounding), and the rows it dropped (computed exactly,
## since L c has finite support). That makes ||r|| / ||F|| a backward error in the
## data, as Rigal-Gaches' is for a matrix, with nothing estimated. A forward error
## would need a bound on L^-1, which is not computed.
cheb_measure <- function(P, c, F, values) {
  p <- P$args$order
  N <- length(c)
  bc <- P$args$boundary
  rb <- vapply(seq_len(p), function(i)
    sum(bc_row(bc$at[i], bc$weights[[i]], N) * c) - values[i], numeric(1))
  Lc <- cheb_apply_L(P, c)
  n <- max(length(Lc), length(F))
  r <- pad_to(Lc, n) - pad_to(F, n)
  kept <- seq_len(min(N - p, n))
  list(finite = sqrt(sum(rb^2) + sum(r[kept]^2)),
       truncation = sqrt(sum(r[-kept]^2)),
       scale = sqrt(sum(F^2) + sum(values^2)))
}

cheb_solve <- function(A, b, args) {
  A <- function_request_operator(A, args$tol)
  P <- A
  cheb_check_problem(P, "solve")
  records <- function_records(A)
  if (is.function(b)) {
    z <- function_coefficients(b, "chebyshev", "b", args$coefficient_tol %||% .Machine$double.eps,
                              args$nmax %||% 65536L)
    b <- z$coefficients
    records <- c(records, list(z))
  }
  f <- as.numeric(b)
  if (!length(f) || anyNA(f)) stopf("b is the Chebyshev coefficients of the right-hand side")
  tol <- args$tol %||% 1e-10
  p <- P$args$order
  F <- cheb_convert_chain(f, 0L, p)
  values <- P$args$boundary$value
  N <- max(CHEB_N_START, length(f) + p + 2L)
  widths <- integer(0)
  repeat {
    widths <- c(widths, N)
    S <- cheb_system(P, N)
    rhs <- c(values, pad_to(F, N - p)[seq_len(N - p)])
    c <- tryCatch(solve(S, rhs), error = function(e)
      stopf("the %d-coefficient system is singular: %s", N, conditionMessage(e)))
    m <- cheb_measure(P, c, F, values)
    if (m$truncation <= tol * m$scale / 2 || N >= CHEB_N_MAX) break
    N <- min(2L * N, CHEB_N_MAX)
  }
  cert <- cheb_solve_certificate(m, tol, N, widths)
  cert <- function_certificate(cert, records)
  structure(c, certificate = cert)
}

cheb_solve_certificate <- function(m, tol, N, widths) {
  r <- cert_rows()
  eps <- .Machine$double.eps
  rel_f <- m$finite / max(m$scale, .Machine$double.xmin)
  rel_t <- m$truncation / max(m$scale, .Machine$double.xmin)
  flr <- SOLVE_FLOOR_CONST * eps * N
  r$add("finite residual", cert_level(rel_f, tol, flr),
        sprintf("||boundary and kept rows of L c - F|| / ||F|| = %.3e on %d coefficients", rel_f, N))
  r$add("truncation residual", cert_level(rel_t, tol, flr),
        sprintf("the rows the truncation dropped, computed exactly: %.3e relative; widths %s",
                rel_t, paste(widths, collapse = " -> ")),
        source = "theorem", guarantee = "deterministic_bound", confidence = 1)
  bw <- sqrt(rel_f^2 + rel_t^2)
  r$add("backward error", cert_level(bw, tol, flr),
        sprintf("the coefficients solve the infinite system with the right-hand side moved by %.3e relative", bw))
  r$add("forward error", "not_checked",
        "a distance to the true solution needs a bound on the inverse of the operator, which is not computed",
        source = "computation", guarantee = "identity", confidence = NA_real_)
  build_certificate(r$collect(), subject = "chebyshev solve",
                    values = list(finite = rel_f, truncation = rel_t, backward_error = bw,
                                  n = N, widths = widths),
                    evidence = r$collect_evidence())
}

## ------------------------------------------------------------------- eigs --

## L u = lambda u with the boundary conditions is [B; L] c = lambda [0; S c], with
## S the conversion to C^(p) and the boundary rows zero on the right: eigenvalues
## at infinity, one per condition. Swapped, S c = mu [B; L] c with mu = 1 / lambda,
## they sit at mu = 0, out of the way of the largest |mu|, which are the smallest
## |lambda|. The residual reported is the original operator's,
## [B c; L c - lambda S c] over every row, exact for the same reason as the solve's,
## and it is measured in ultraspherical coefficients rather than in L2, so it bounds
## no distance to the spectrum: the forward row says so.
cheb_eigs <- function(A, k, control, args) {
  A <- function_request_operator(A, args$tol)
  P <- A
  cheb_check_problem(P, "eigs")
  if (identical(args$method, "galerkin")) return(cheb_galerkin_eigs(P, k, args))
  p <- P$args$order
  tol <- args$tol %||% 1e-10
  which <- args$which %||% "smallest"
  if (!which %in% c("smallest", "SM")) {
    stopf("eigs() on a Chebyshev operator finds the eigenvalues of smallest magnitude; which = \"smallest\"")
  }
  N <- max(CHEB_N_START, 4L * k + p)
  widths <- integer(0)
  prev <- NULL
  repeat {
    widths <- c(widths, N)
    S <- cheb_system(P, N)
    Bm <- matrix(0, N, N)
    for (j in seq_len(N)) {
      e <- numeric(N); e[j] <- 1
      Bm[p + seq_len(N - p), j] <- pad_to(cheb_convert_chain(e, 0L, p), N - p)[seq_len(N - p)]
    }
    K <- solve(S, Bm)
    fit <- eigs(linop(K), k = k, which = "largest", tol = min(tol, 1e-12), method = "arnoldi")
    mu <- fit$values
    lam <- 1 / mu
    X <- fit$vectors
    meas <- lapply(seq_along(lam), function(i) {
      c <- X[, i]
      Lc <- cheb_apply_L(P, c)
      Sc <- cheb_convert_chain(c, 0L, p)
      n <- max(length(Lc), length(Sc))
      rr <- pad_to(Lc, n) - lam[i] * pad_to(Sc, n)
      bc <- P$args$boundary
      rb <- vapply(seq_len(p), function(j) sum(bc_row(bc$at[j], bc$weights[[j]], N) * c),
                   numeric(1))
      kept <- seq_len(min(N - p, n))
      sc <- sqrt(sum(Mod(Sc)^2))
      c(finite = sqrt(sum(Mod(rb)^2) + sum(Mod(rr[kept])^2)) / sc,
        truncation = sqrt(sum(Mod(rr[-kept])^2)) / sc)
    })
    meas <- do.call(rbind, meas)
    ok <- meas[, "truncation"] <= tol * pmax(Mod(lam), 1) / 2
    stable <- !is.null(prev) && length(prev) == length(lam) &&
              max(Mod(sort_by_mod(lam) - sort_by_mod(prev)) / pmax(Mod(lam), 1)) <= tol
    if ((all(ok) && stable) || N >= CHEB_N_MAX) break
    prev <- lam
    N <- min(2L * N, CHEB_N_MAX)
  }
  ord <- order(Mod(lam))
  lam <- lam[ord]; X <- X[, ord, drop = FALSE]; meas <- meas[ord, , drop = FALSE]
  if (all(Im(lam) == 0)) lam <- Re(lam)
  cert <- cheb_eigen_certificate(lam, meas, tol, N, widths, k)
  cert <- function_certificate(cert, function_records(A), P)
  new_eigen_result(values = lam, vectors = X, certificate = cert, k = k, which = which,
                   method = "arnoldi on the ultraspherical pencil", iterations = fit$iterations,
                   maxit = fit$maxit, rounds = fit$restarts,
                   converged = cert$values$converged, residuals = meas[, "finite"],
                   dispatch = list(requested = args$method %||% "auto",
                                   chosen = "arnoldi on the ultraspherical pencil",
                                   reason = sprintf("Chebyshev operator of order %d; %s coefficients",
                                                    P$args$order, paste(widths, collapse = " -> ")),
                                   candidates = NULL),
                   engine = NULL)
}

sort_by_mod <- function(z) z[order(Mod(z), Re(z), Im(z))]

cheb_eigen_certificate <- function(lam, meas, tol, N, widths, k) {
  r <- cert_rows()
  scale <- pmax(Mod(lam), 1)
  flr <- SOLVE_FLOOR_CONST * .Machine$double.eps * N
  fin <- meas[, "finite"] / scale; trn <- meas[, "truncation"] / scale
  met <- sqrt(fin^2 + trn^2) <= tol + flr
  r$add("finite residual", cert_level(max(fin), tol, flr),
        sprintf("worst ||[B c; L c - lambda S c]|| on the kept rows / (max(|lambda|, 1) ||S c||) = %.3e", max(fin)))
  r$add("truncation residual", cert_level(max(trn), tol, flr),
        sprintf("the dropped rows, computed exactly: %.3e; %s coefficients", max(trn),
                paste(widths, collapse = " -> ")),
        source = "theorem", guarantee = "deterministic_bound", confidence = 1)
  r$add("target identity", "not_checked",
        "a residual does not say which eigenvalue was found",
        source = "computation", guarantee = "identity", confidence = NA_real_)
  cert_add_convergence(r, met, k, N, CHEB_N_MAX, NULL, "the coefficients stopped converging")
  r$add("forward error", "not_checked",
        "the residual is measured in ultraspherical coefficients, not in L2, so it bounds no distance to the spectrum",
        source = "computation", guarantee = "identity", confidence = NA_real_)
  build_certificate(r$collect(), subject = "chebyshev eigenproblem",
                    values = list(finite = fin, truncation = trn, n = N, widths = widths,
                                  converged = met),
                    evidence = r$collect_evidence())
}

register_chebyshev_node <- function() {
  linop_register_node("chebyshev", cheb_apply, NULL, function(op) Inf,
                      spectrum = cheb_eigs, solve = cheb_solve, overwrite = TRUE)
}
