## Functions on the line and the half line, held as coefficients
## in the Hermite and Laguerre functions.
##
## Both bases make multiplication by x tridiagonal, and each has a differential
## operator that is diagonal in it. For the Hermite functions psi_n on R,
##     x psi_n = sqrt(n/2) psi_{n-1} + sqrt((n+1)/2) psi_{n+1},
##     psi_n'  = sqrt(n/2) psi_{n-1} - sqrt((n+1)/2) psi_{n+1},
## so d/dx is tridiagonal too and -d^2 + x^2 is diagonal, 2n + 1. For the Laguerre
## functions phi_n = exp(-x/2) L_n(x) on [0, Inf),
##     x phi_n = -n phi_{n-1} + (2n + 1) phi_n - (n + 1) phi_{n+1},
## and T u = -(x u')' + x u / 4 is diagonal, n + 1/2. So -u'' + V u with V a
## polynomial (Hermite) and c T u + V u (Laguerre) are banded operators on
## sequences indexed by N, and their truncations certify by the same edge
## identity as the circle's: banded-space.R.

## X and D (Hermite) or X and T (Laguerre) on the first M coefficients.
basis_matrices <- function(kind, M) {
  n <- seq_len(M - 1L)
  rec <- basis_recurrence(kind, M)
  X <- matrix(0, M, M)
  diag(X) <- rec$b
  X[cbind(n, n + 1L)] <- rec$a[n]; X[cbind(n + 1L, n)] <- rec$a[n]
  if (kind == "hermite") {
    D <- matrix(0, M, M)
    D[cbind(n, n + 1L)] <- sqrt(n / 2); D[cbind(n + 1L, n)] <- -sqrt(n / 2)
    list(X = X, D = D)
  } else {
    list(X = X, T = diag(seq.int(0L, M - 1L) + 0.5, M))
  }
}

mat_power <- function(A, m) {
  out <- diag(1, nrow(A))
  for (i in seq_len(m)) out <- out %*% A
  out
}

## The operator's matrix on the first M coefficients, exact: products of
## tridiagonal matrices of degree m truncated at M + m agree with the infinite
## products on the first M rows and columns, so the block is computed on a
## padded matrix and cut.
basis_block <- function(P, M) {
  pad <- P$args$band + 2L
  Mb <- M + pad
  B <- basis_matrices(P$args$kind, Mb)
  L <- matrix(0, Mb, Mb)
  if (P$args$kind == "hermite") {
    for (m in seq_along(P$args$coef) - 1L) {
      if (P$args$coef[m + 1L] != 0) L <- L + P$args$coef[m + 1L] * mat_power(B$D, m)
    }
  } else {
    L <- L + P$args$coef[1L] * B$T
  }
  L <- basis_polynomial_matrix(P$args$v, B$X, P$args$kind,
                                P$args$multiplier_basis %||% "monomial", L)
  L[seq_len(M), seq_len(M), drop = FALSE]
}

basis_apply <- function(op, X, mode) {
  stopf(paste0("an operator on %s coefficients has no block to act on.\n",
               "  Truncate it with finite_section(A, n) or call eigs() on it."), op$args$kind)
}

new_basis_operator <- function(kind, derivative, multiplier, multiplier_basis = "monomial") {
  v <- as.numeric(multiplier)
  while (length(v) > 1L && v[length(v)] == 0) v <- v[-length(v)]
  deg <- length(v) - 1L
  leading <- if (multiplier_basis == "orthogonal" && deg > 0L)
    v[deg + 1L] / prod(basis_recurrence(kind, deg)$a) else v[deg + 1L]
  if (kind == "hermite") {
    coef <- as.numeric(derivative)
    while (length(coef) > 1L && coef[length(coef)] == 0) coef <- coef[-length(coef)]
    order <- length(coef) - 1L
    odd <- seq_along(coef) %% 2L == 0L
    self_adjoint <- all(coef[odd] == 0)
    ## -c u'' + V u with c > 0 and V of even degree with a positive leading
    ## coefficient confines: compact resolvent, discrete spectrum. A free particle
    ## (V constant) has continuous spectrum.
    discrete <- self_adjoint && order == 2L && coef[3L] < 0 && deg >= 2L &&
                deg %% 2L == 0L && leading > 0
    band <- max(order, deg)
    lead <- abs(coef[length(coef)]) + abs(leading)
  } else {
    if (length(derivative) != 1L || !is.numeric(derivative)) {
      stopf(paste0("for kind = \"laguerre\", derivative is the single coefficient of the\n",
                   "  Laguerre operator -(x u')' + x u / 4, which is diagonal in that basis"))
    }
    coef <- as.numeric(derivative)
    order <- 2L
    self_adjoint <- TRUE
    discrete <- coef > 0 && (deg == 0L || leading > 0)
    band <- max(1L, deg)
    lead <- abs(coef) + abs(leading)
  }
  caps <- list()
  if (self_adjoint) {
    caps$hermitian <- capability(TRUE, ev_construction())
    caps$symmetric <- capability(TRUE, ev_construction())
  }
  new_linop("basis", c(Inf, Inf), "double", do.call(new_caps, caps),
            list(kind = kind, coef = coef, order = order, v = v, degree = deg,
                 band = band, self_adjoint = self_adjoint, discrete = discrete,
                 lead = lead, multiplier_basis = multiplier_basis),
            cost = Inf)
}

## ------------------------------------------------------------- the section --

## The indices 0 to n, n + 1 coefficients.
basis_section <- function(P, n) {
  if (n < P$args$band) {
    stopf("n = %d is narrower than the band %d; the section has to hold a whole band",
          n, P$args$band)
  }
  M <- basis_block(P, n + 1L)
  caps <- list()
  if (P$args$self_adjoint) {
    caps$hermitian <- capability(TRUE, ev_construction())
    caps$symmetric <- capability(TRUE, ev_construction())
  }
  new_linop("basis_section", c(n + 1L, n + 1L), "double", do.call(new_caps, caps),
            list(A = P, n = n, M = M), cost = 2 * (n + 1) * (2 * P$args$band + 1),
            finite = TRUE)
}

basis_section_apply <- function(op, X, mode) {
  M <- op$args$M
  if (mode_transposes(mode)) crossprod(M, X) else M %*% X
}

basis_section_materialize <- function(op) op$args$M

## The rows past the cut that the band reaches, from an exact block one band wider.
basis_edge <- function(P, X) {
  N <- nrow(X)
  r <- P$args$band
  big <- basis_block(P, N + r)
  out <- big[N + seq_len(r), seq_len(N), drop = FALSE] %*% X
  col_norms(out)
}

basis_certify <- function(ctx) {
  A <- ctx$A
  P <- A$args$A
  band_certificate(A, ctx$values, ctx$vectors, ctx$tol, ctx$iterations, ctx$maxit,
                   ctx$floor_const, ctx$requested, ctx$stop_reason, P,
                   hn = max(rowSums(abs(A$args$M))),
                   scale = function(q) Mod(q) + P$args$lead,
                   edge = function(X) basis_edge(P, X),
                   self_adjoint = P$args$self_adjoint, discrete = P$args$discrete,
                   reach = sprintf("the %d coefficients past the cut the band reaches", P$args$band),
                   space = if (P$args$kind == "hermite") "the line" else "the half line")
}

basis_eigs <- function(A, k, control, args) {
  band_eigs(A, k, control, args, n_min = max(A$args$band, k + 1L), space = "l^2(N)")
}

register_basis_nodes <- function() {
  linop_register_node("basis", basis_apply, NULL, function(op) Inf,
                      spectrum = basis_eigs, section = basis_section, overwrite = TRUE)
  linop_register_node("basis_section", basis_section_apply, basis_section_materialize,
                      function(op) op$cost, certify = basis_certify, overwrite = TRUE)
}
