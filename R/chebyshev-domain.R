## Plan sections 8.4 to 8.6 for the Chebyshev operator: its domain, its adjoint
## with the adjoint domain, and a Galerkin discretisation of its weak form.
##
## The domain of L = sum_m c_m D^m + a(x) on [-1, 1] is the set of u satisfying
## its homogeneous boundary conditions. For constant derivative coefficients and
## a polynomial multiplier, integration by parts gives the adjoint exactly
## (plan 8.4 names this class as the tractable one):
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

## B(u, v) as a matrix on (u, u', ..., u^(p-1)) x (v, ..., v^(p-1)).
cheb_boundary_form <- function(coef, p) {
  M <- matrix(0, p, p)
  for (m in seq_along(coef) - 1L) {
    if (m == 0L || coef[m + 1L] == 0) next
    for (j in 0:(m - 1L)) {
      i <- m - 1L - j
      M[i + 1L, j + 1L] <- M[i + 1L, j + 1L] + coef[m + 1L] * (-1)^j
    }
  }
  M
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
  M <- cheb_boundary_form(P$args$coef, p)
  at <- numeric(0); weights <- list()
  for (e in c(-1, 1)) {
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
  coef <- P$args$coef * (-1)^(seq_along(P$args$coef) - 1L)
  new_chebyshev(coef, P$args$a, cheb_adjoint_boundary(P))
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
  odd <- seq_along(P$args$coef) %% 2L == 0L
  if (any(P$args$coef[odd] != 0)) return(FALSE)
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

## The weak form of L u = -c u'' + (c0 + a(x)) u with homogeneous Dirichlet
## conditions, a(u, v) = c int u' v' + int (c0 + a) u v against m(u, v) = int u v,
## is the sesquilinear form plan 8.5 asks to treat as a peer of the operator. It
## is discretised on phi_j = T_{j+2} - T_j, which vanish at +-1 exactly, so the
## trial space lies inside the form domain H^1_0. Gauss-Legendre quadrature on
## enough nodes integrates every entry exactly, the integrands being polynomials.
## The pencil (A_n, M_n) is definite, M_n being the Gram matrix of independent
## functions, and goes to eigs(B = ) in the B inner product. The Ritz vectors
## lift back to Chebyshev coefficients, and the Galerkin values are upper bounds
## of the true eigenvalues by Courant-Fischer, since the trial space is a
## subspace of the form domain.
gauss_legendre <- function(Q) {
  k <- seq_len(Q - 1L)
  b <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, Q, Q)
  J[cbind(k, k + 1L)] <- b; J[cbind(k + 1L, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1L, ]^2)
}

cheb_galerkin_ok <- function(P) {
  bc <- P$args$boundary
  cf <- c(P$args$coef, 0, 0, 0)
  if (P$args$order != 2L || cf[2L] != 0) {
    return("it is for -c u'' + (c0 + a) u, a second-order operator with no first-order term")
  }
  if (is.null(bc) || length(bc$at) != 2L || !setequal(bc$at, c(-1, 1)) ||
      !all(vapply(bc$weights, function(w) length(w) >= 1L && w[1L] != 0 && all(w[-1L] == 0), logical(1))) ||
      any(bc$value != 0)) {
    return("it needs homogeneous Dirichlet conditions, u(-1) = u(1) = 0")
  }
  TRUE
}

cheb_galerkin_matrices <- function(P, n) {
  c2 <- P$args$coef[3L]; c0 <- P$args$coef[1L]
  a <- P$args$a
  Q <- n + 4L + ceiling((length(a) - 1L) / 2)
  g <- gauss_legendre(Q)
  th <- acos(pmin(1, pmax(-1, g$x)))
  Tk <- function(k) cos(k * th)
  dTk <- function(k) if (k == 0L) rep(0, Q) else k * sin(k * th) / sin(th)
  Phi <- sapply(0:(n - 1L), function(j) Tk(j + 2L) - Tk(j))
  dPhi <- sapply(0:(n - 1L), function(j) dTk(j + 2L) - dTk(j))
  ax <- c0 + vapply(g$x, function(x) sum(a * cos((seq_along(a) - 1L) * acos(x))), numeric(1))
  K <- (-c2) * crossprod(dPhi, g$w * dPhi) + crossprod(Phi, (g$w * ax) * Phi)
  Mm <- crossprod(Phi, g$w * Phi)
  list(K = (K + t(K)) / 2, M = (Mm + t(Mm)) / 2)
}

cheb_lift <- function(x) {
  n <- length(x)
  cvec <- numeric(n + 2L)
  cvec[seq_len(n) + 2L] <- cvec[seq_len(n) + 2L] + x
  cvec[seq_len(n)] <- cvec[seq_len(n)] - x
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
    U <- apply(fit$vectors, 2L, cheb_lift)
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
        "each theta_i is at least the i-th eigenvalue of the operator (Courant-Fischer: the trial space lies in the form domain H^1_0), to the arithmetic floor",
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
  new_eigen_result(values = lam, vectors = U, certificate = cert, k = k,
                   which = "smallest_algebraic", method = "galerkin",
                   iterations = fit$iterations, maxit = fit$maxit, rounds = fit$restarts,
                   converged = met, residuals = strong,
                   dispatch = list(requested = "galerkin", chosen = "galerkin",
                                   reason = sprintf("weak form on phi_j = T_{j+2} - T_j, Lanczos in the mass inner product; %s basis functions",
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
  switch(A$node,
    chebyshev = cheb_adjoint(A),
    fourier = linop_space("fourier",
                          derivative = Conj(A$args$coef) * (-1)^(seq_along(A$args$coef) - 1L),
                          multiplier = Conj(rev(A$args$vhat))),
    basis = new_basis_operator(A$args$kind,
                               if (A$args$kind == "hermite")
                                 A$args$coef * (-1)^(seq_along(A$args$coef) - 1L)
                               else A$args$coef,
                               A$args$v))
}
