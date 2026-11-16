## The Chebyshev operator's its domain, its adjoint
## with the adjoint domain, and a Galerkin discretisation of its weak form.
##
## The domain of L = sum_m c_m D^m + a(x) on [-1, 1] is the set of u satisfying
## its homogeneous boundary conditions. For constant derivative coefficients and
## a polynomial multiplier, integration by parts gives the adjoint exactly
## (the tractable class):
##
##     int (L u) conj(v) = int u conj(L* v) + [B(u, v)]_{-1}^{1},
##     L* = sum_m (-1)^m conj(c_m) D^m + conj(a),
##     B(u, v) = sum_m c_m sum_{j<m} (-1)^j u^(m-1-j) conj(v)^(j),
##
## and the adjoint conditions at each endpoint are those that make B vanish for
## every boundary value u can take there: with W a basis of the admissible
## (u, u', ...) at the endpoint and M the matrix of B, they are the rows of W^T M.
## They are mixed (Robin) conditions in general. Self-adjoint means L* = L as
## coefficients and the two sets of conditions span the same rows at each end.
## That is self-adjointness in L2 on the interval, which is not the hermitian
## capability: the Chebyshev coefficients are not orthonormal in L2, so the
## coefficient matrix of a self-adjoint operator is not hermitian.

## B(u, v) at one endpoint as a matrix on (u, ..., u^(p-1)) x (v, ..., v^(p-1)).
## With variable coefficients the boundary form carries their derivatives,
##     B(u, v) = sum_m sum_{j<m} (-1)^j u^(m-1-j) (a_m v)^(j),
##     (a_m v)^(j) = sum_{i<=j} choose(j, i) a_m^(j-i) v^(i),
## so the matrix differs between the two ends.
cheb_boundary_form <- function(terms, p, at) {
  M <- matrix(0, p, p)
  for (m in seq_along(terms) - 1L) {
    am <- terms[[m + 1L]]
    if (m == 0L || all(am == 0)) next
    for (j in 0:(m - 1L)) {
      r <- m - 1L - j
      for (i in 0:j) {
        M[r + 1L, i + 1L] <- M[r + 1L, i + 1L] +
          (-1)^j * choose(j, i) * cheb_value_at(am, at, j - i)
      }
    }
  }
  M
}

## The formal adjoint's coefficients, L* v = sum_m (-1)^m (a_m v)^(m) expanded:
## b_j = sum_{m >= j} (-1)^m choose(m, j) a_m^(m-j).
cheb_adjoint_terms <- function(terms) {
  p <- length(terms) - 1L
  lapply(0:p, function(j) {
    out <- 0
    for (m in j:p) {
      d <- cheb_dcoef_n(terms[[m + 1L]], m - j)
      n <- max(length(out), length(d))
      out <- pad_to(out, n) + (-1)^m * choose(m, j) * pad_to(d, n)
    }
    trim_poly(out)
  })
}

null_basis <- function(R, p) {
  if (!length(R) || !nrow(R)) return(diag(p))
  s <- svd(R, nu = 0L, nv = p)
  rank <- sum(s$d > 1e-12 * max(s$d, 1))
  if (rank >= p) return(matrix(0, p, 0L))
  s$v[, (rank + 1L):p, drop = FALSE]
}

cheb_endpoint_rows <- function(bc, at, p) {
  i <- which(bc$at == at)
  if (!length(i)) return(matrix(0, 0L, p))
  do.call(rbind, lapply(bc$weights[i], function(w) c(w, numeric(p - length(w)))[seq_len(p)]))
}

## The adjoint conditions, at -1 and 1, as weight vectors with value 0.
cheb_adjoint_boundary <- function(P) {
  p <- P$args$order
  bc <- P$args$boundary
  at <- numeric(0); weights <- list()
  for (e in c(-1, 1)) {
    M <- cheb_boundary_form(P$args$terms, p, e)
    W <- null_basis(cheb_endpoint_rows(bc, e, p), p)
    if (!ncol(W)) next
    rows <- t(W) %*% M
    keep <- apply(rows, 1L, function(r) max(abs(r)) > 1e-12)
    rows <- rows[keep, , drop = FALSE]
    for (i in seq_len(nrow(rows))) {
      w <- rows[i, ]
      w <- w / w[which(abs(w) > 1e-12 * max(abs(w)))[1L]]
      at <- c(at, e); weights <- c(weights, list(w))
    }
  }
  list(at = at, weights = weights, value = numeric(length(at)))
}

cheb_adjoint <- function(P) {
  if (is.null(P$args$boundary)) {
    stopf(paste0("the adjoint of a Chebyshev operator needs its domain, which its boundary\n",
                 "  conditions define; this one has none"))
  }
  new_chebyshev(cheb_adjoint_terms(P$args$terms), 0, cheb_adjoint_boundary(P))
}

same_rowspace <- function(A, B) {
  if (!nrow(A) && !nrow(B)) return(TRUE)
  if (nrow(A) != nrow(B)) return(FALSE)
  rk <- function(X) { d <- svd(X, nu = 0L, nv = 0L)$d; sum(d > 1e-10 * max(d, 1)) }
  rk(rbind(A, B)) == rk(A) && rk(A) == rk(B)
}

## Self-adjoint in L2: the formal adjoint has the same coefficients and the
## adjoint conditions span the operator's own at each endpoint.
cheb_self_adjoint <- function(P) {
  bc <- P$args$boundary
  p <- P$args$order
  if (is.null(bc) || p < 1L || length(bc$at) != p) return(NA)
  adj <- cheb_adjoint_terms(P$args$terms)
  same_poly <- function(x, y) {
    n <- max(length(x), length(y))
    max(abs(pad_to(x, n) - pad_to(y, n))) <= 1e-12 * max(1, abs(x), abs(y))
  }
  if (length(adj) != length(P$args$terms) ||
      !all(mapply(same_poly, adj, P$args$terms))) return(FALSE)
  adj <- cheb_adjoint_boundary(P)
  all(vapply(c(-1, 1), function(e)
    same_rowspace(cheb_endpoint_rows(bc, e, p), cheb_endpoint_rows(adj, e, p)), logical(1)))
}

cheb_domain_text <- function(P) {
  bc <- P$args$boundary
  if (is.null(bc)) return("no boundary conditions")
  term <- function(w, at) {
    parts <- character(0)
    for (d in seq_along(w) - 1L) {
      if (abs(w[d + 1L]) < 1e-14) next
      c1 <- w[d + 1L]
      lead <- if (abs(c1 - 1) < 1e-14) "" else if (abs(c1 + 1) < 1e-14) "-"
              else paste0(format(signif(c1, 4)), " ")
      parts <- c(parts, sprintf("%su%s(%d)", lead, strrep("'", d), as.integer(at)))
    }
    gsub("+ -", "- ", paste(parts, collapse = " + "), fixed = TRUE)
  }
  paste0("{u : ", paste(vapply(seq_along(bc$at), function(i)
    paste(term(bc$weights[[i]], bc$at[i]), "= 0"), character(1)), collapse = ", "), "}")
}

## ---------------------------------------------------------------- Galerkin --

## The weak form of the self-adjoint L u = -(k u')' + q u, with k = -a_2 and
## q = a_0 read off the coefficients (self-adjointness makes a_1 = a_2'), is
##
##     a(u, v) = int k u' v' + int q u v + boundary terms,   m(u, v) = int u v,
##
## the sesquilinear form, treated as a peer of the operator. A
## Dirichlet end is an essential condition and the basis vanishes there; a
## Robin or Neumann end, alpha u + beta u' = 0 with beta nonzero, is natural:
## integrating by parts leaves -[k u' v] at the ends, which the condition turns
## into + k(1) (alpha/beta) u(1) v(1) - k(-1) (alpha/beta) u(-1) v(-1), and the
## basis is unconstrained there. The basis is T_{j+2} - T_j for two Dirichlet
## ends, T_{j+1} - T_j or T_{j+1} + T_j for one at 1 or at -1, T_j for none, so the
## trial space lies in the form domain. Gauss-Legendre quadrature on enough nodes
## integrates every entry exactly, the integrands being polynomials. The pencil
## is definite (M_n is the Gram matrix of independent functions) and goes to
## eigs(B = ); the vectors lift back to Chebyshev coefficients, and the values
## are upper bounds of the true eigenvalues by Courant-Fischer.
gauss_legendre <- function(Q) {
  k <- seq_len(Q - 1L)
  b <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, Q, Q)
  J[cbind(k, k + 1L)] <- b; J[cbind(k + 1L, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1L, ]^2)
}

## Each end's condition, as "dirichlet" or the Robin ratio alpha / beta.
cheb_galerkin_ends <- function(P) {
  bc <- P$args$boundary
  lapply(c(-1, 1), function(e) {
    i <- which(bc$at == e)
    w <- c(bc$weights[[i]], 0, 0)
    if (w[2L] == 0) list(kind = "dirichlet") else list(kind = "natural", ratio = w[1L] / w[2L])
  })
}

cheb_galerkin_ok <- function(P) {
  bc <- P$args$boundary
  if (P$args$order != 2L) return("it is for second-order operators")
  if (is.null(bc) || length(bc$at) != 2L || !setequal(bc$at, c(-1, 1))) {
    return("it needs one condition at each end")
  }
  if (any(bc$value != 0)) return("it needs homogeneous conditions")
  if (!isTRUE(P$args$self_adjoint)) {
    return("it is for -(k u')' + q u with self-adjoint conditions, the operators with a symmetric weak form")
  }
  g <- gauss_legendre(40L)
  kx <- -vapply(g$x, function(x) sum(P$args$terms[[3L]] * cos((seq_along(P$args$terms[[3L]]) - 1L) * acos(x))), numeric(1))
  if (any(kx <= 0)) return("it needs k = -a_2 positive on [-1, 1], for the form to be bounded below")
  TRUE
}

cheb_galerkin_basis <- function(ends, j) {
  d <- vapply(ends, function(e) e$kind == "dirichlet", logical(1))
  if (all(d)) list(hi = j + 2L, lo = j, sign = -1)
  else if (d[2L]) list(hi = j + 1L, lo = j, sign = -1)
  else if (d[1L]) list(hi = j + 1L, lo = j, sign = 1)
  else list(hi = j, lo = NA_integer_, sign = 0)
}

cheb_galerkin_matrices <- function(P, n) {
  ends <- cheb_galerkin_ends(P)
  kcf <- -P$args$terms[[3L]]
  qcf <- P$args$terms[[1L]]
  Q <- n + 6L + ceiling(max(length(kcf), length(qcf)) / 2)
  g <- gauss_legendre(Q)
  th <- acos(pmin(1, pmax(-1, g$x)))
  Tk <- function(k) cos(k * th)
  dTk <- function(k) if (k == 0L) rep(0, Q) else k * sin(k * th) / sin(th)
  at_end <- function(k, e) if (e > 0) 1 else (-1)^k
  phi <- function(j, f) {
    b <- cheb_galerkin_basis(ends, j)
    if (b$sign == 0) f(b$hi) else f(b$hi) + b$sign * f(b$lo)
  }
  Phi <- sapply(0:(n - 1L), phi, f = Tk)
  dPhi <- sapply(0:(n - 1L), phi, f = dTk)
  series <- function(cf) vapply(g$x, function(x) sum(cf * cos((seq_along(cf) - 1L) * acos(x))), numeric(1))
  K <- crossprod(dPhi, (g$w * series(kcf)) * dPhi) + crossprod(Phi, (g$w * series(qcf)) * Phi)
  for (e in 1:2) {
    if (ends[[e]]$kind != "natural") next
    x0 <- c(-1, 1)[e]
    pv <- vapply(0:(n - 1L), phi, numeric(1), f = function(k) at_end(k, x0))
    ke <- cheb_value_at(kcf, x0)
    K <- K + (if (x0 > 0) 1 else -1) * ke * ends[[e]]$ratio * outer(pv, pv)
  }
  Mm <- crossprod(Phi, g$w * Phi)
  list(K = (K + t(K)) / 2, M = (Mm + t(Mm)) / 2, ends = ends)
}

cheb_lift <- function(x, ends) {
  n <- length(x)
  cvec <- numeric(n + 2L)
  for (j in seq_len(n) - 1L) {
    b <- cheb_galerkin_basis(ends, j)
    cvec[b$hi + 1L] <- cvec[b$hi + 1L] + x[j + 1L]
    if (b$sign != 0) cvec[b$lo + 1L] <- cvec[b$lo + 1L] + b$sign * x[j + 1L]
  }
  cvec
}

cheb_galerkin_eigs <- function(P, k, args) {
  ok <- cheb_galerkin_ok(P)
  if (!isTRUE(ok)) stopf("eigs(method = \"galerkin\") cannot run: %s", ok)
  tol <- args$tol %||% 1e-10
  p <- 2L
  n <- max(16L, 2L * k + 4L)
  widths <- integer(0); prev <- NULL
  repeat {
    widths <- c(widths, n)
    G <- cheb_galerkin_matrices(P, n)
    B <- linop(G$M)
    ## the Gram matrix of independent functions is positive definite by
    ## construction, at any size
    B$caps$positive_definite <- capability(TRUE, ev_construction())
    fit <- eigs(linop(G$K), k = k, which = "smallest_algebraic", B = B,
                tol = min(tol, 1e-12), method = "lanczos")
    lam <- fit$values
    U <- apply(fit$vectors, 2L, cheb_lift, ends = G$ends)
    strong <- vapply(seq_len(k), function(i) {
      cvec <- U[, i]
      Lc <- cheb_apply_L(P, cvec); Sc <- cheb_convert_chain(cvec, 0L, p)
      m <- max(length(Lc), length(Sc))
      sqrt(sum((pad_to(Lc, m) - lam[i] * pad_to(Sc, m))^2)) / sqrt(sum(Sc^2)) / max(abs(lam[i]), 1)
    }, numeric(1))
    stable <- !is.null(prev) && max(abs(lam - prev) / pmax(abs(lam), 1)) <= tol
    if ((all(strong <= tol) && stable) || n >= CHEB_N_MAX) break
    prev <- lam
    n <- min(2L * n, CHEB_N_MAX)
  }
  r <- cert_rows()
  flr <- SOLVE_FLOOR_CONST * .Machine$double.eps * n
  r$add("galerkin residual", cert_level(max(fit$residuals / pmax(abs(lam), 1)), tol, flr),
        sprintf("worst ||A_n x - theta M_n x|| relative = %.3e on %d basis functions", max(fit$residuals), n))
  r$add("strong residual", cert_level(max(strong), tol, flr),
        sprintf("the lifted function in the strong form, ||L u - theta u|| in C^(2) coefficients, relative: %.3e; widths %s",
                max(strong), paste(widths, collapse = " -> ")))
  r$add("upper bound", "pass",
        "each theta_i is at least the i-th eigenvalue of the operator (Courant-Fischer: the trial space lies in the form domain), to the arithmetic floor",
        source = "theorem", guarantee = "deterministic_bound", confidence = 1)
  r$add("target identity", "not_checked",
        "the upper bound pairs theta_i with lambda_i by rank; a lower bound is not computed",
        source = "computation", guarantee = "identity", confidence = NA_real_)
  met <- strong <= tol + flr
  cert_add_convergence(r, met, k, n, CHEB_N_MAX, NULL, "the basis stopped converging")
  r$add("forward error", "not_checked",
        "the strong residual is in ultraspherical coefficients, not in L2, so it bounds no distance to the spectrum; the upper bound row is the one-sided statement",
        source = "computation", guarantee = "identity", confidence = NA_real_)
  cert <- build_certificate(r$collect(), subject = "galerkin eigenproblem",
                            values = list(strong = strong, n = n, widths = widths, converged = met),
                            evidence = r$collect_evidence())
  cert <- function_certificate(cert, function_records(P), P)
  new_eigen_result(values = lam, vectors = U, certificate = cert, k = k,
                   which = "smallest_algebraic", method = "galerkin",
                   iterations = fit$iterations, maxit = fit$maxit, rounds = fit$restarts,
                   converged = met, residuals = strong,
                   dispatch = list(requested = "galerkin", chosen = "galerkin",
                                   reason = sprintf("weak form on a Chebyshev basis vanishing at the Dirichlet ends, Lanczos in the mass inner product; %s basis functions",
                                                    paste(widths, collapse = " -> ")),
                                   candidates = NULL),
                   engine = NULL)
}

## ------------------------------------------------- adjoints of the other kinds

## Fourier, Hermite and Laguerre coefficients are orthonormal in L2, so the
## adjoint is the coefficients': for the circle the symbol conjugated (c_m to
## (-1)^m conj(c_m)) and V to conj(V); in Hermite functions the odd derivatives
## change sign; the Laguerre operator and real multipliers are their own.
space_adjoint <- function(A) {
  out <- switch(A$node,
    chebyshev = cheb_adjoint(A),
    fourier = linop_space("fourier",
                          derivative = Conj(A$args$coef) * (-1)^(seq_along(A$args$coef) - 1L),
                          multiplier = Conj(rev(A$args$vhat))),
    basis = new_basis_operator(A$args$kind,
                               if (A$args$kind == "hermite")
                                 A$args$coef * (-1)^(seq_along(A$args$coef) - 1L)
                               else A$args$coef,
                               A$args$v, A$args$multiplier_basis %||% "monomial"))
  carry_function_records(out, function_records(A))
}
