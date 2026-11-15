## Phase 6 step 2: functions on [-1, 1] held as Chebyshev coefficients, and the
## ultraspherical method for differential operators on them (Olver and Townsend,
## SIAM Review 2013). dev_notes/phase6-function-spaces-design.md, section 4.
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

## L c in C^(p), exact and of finite length.
cheb_apply_L <- function(P, c) {
  p <- P$args$order
  out <- numeric(0)
  for (m in seq_along(P$args$coef) - 1L) {
    cm <- P$args$coef[m + 1L]
    if (cm == 0) next
    term <- cm * cheb_convert_chain(cheb_diff(c, m), m, p)
    n <- max(length(out), length(term))
    out <- pad_to(out, n) + pad_to(term, n)
  }
  if (any(P$args$a != 0)) {
    term <- cheb_convert_chain(cheb_mult(P$args$a, c), 0L, p)
    n <- max(length(out), length(term))
    out <- pad_to(out, n) + pad_to(term, n)
  }
  out
}

## The N x N system on the first N coefficients: the boundary rows, then the
## first N - p rows of L. Columns are L applied to unit vectors, so the band is
## read off the operator rather than restated.
cheb_system <- function(P, N) {
  p <- P$args$order
  bc <- P$args$boundary
  A <- matrix(0, N, N)
  for (i in seq_len(p)) A[i, ] <- cheb_boundary_row(bc$at[i], bc$derivative[i], N)
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

new_chebyshev <- function(derivative, multiplier, boundary) {
  coef <- derivative
  while (length(coef) > 1L && coef[length(coef)] == 0) coef <- coef[-length(coef)]
  p <- length(coef) - 1L
  if (p < 1L) stopf("a Chebyshev operator needs a derivative term; derivative = %s is order 0",
                    deparse(derivative))
  if (is.null(boundary) || !is.list(boundary) ||
      !all(c("at", "derivative", "value") %in% names(boundary))) {
    stopf("boundary is a list of at (each -1 or 1), derivative (whole numbers) and value")
  }
  nb <- length(boundary$at)
  if (nb != p || length(boundary$derivative) != nb || length(boundary$value) != nb) {
    stopf("an operator of order %d needs %d boundary conditions; got %d", p, p, nb)
  }
  if (!all(boundary$at %in% c(-1, 1))) stopf("boundary conditions sit at -1 or 1")
  if (any(boundary$derivative < 0) || any(boundary$derivative >= p)) {
    stopf("boundary derivatives run from 0 to %d for an operator of order %d", p - 1L, p)
  }
  new_linop("chebyshev", c(Inf, Inf), "double", new_caps(),
            list(coef = as.numeric(coef), order = p, a = as.numeric(multiplier),
                 boundary = list(at = as.numeric(boundary$at),
                                 derivative = as.integer(boundary$derivative),
                                 value = as.numeric(boundary$value)),
                 domain = "chebyshev", range = sprintf("ultraspherical(%d)", p)),
            cost = Inf)
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
    sum(cheb_boundary_row(bc$at[i], bc$derivative[i], N) * c) - values[i], numeric(1))
  Lc <- cheb_apply_L(P, c)
  n <- max(length(Lc), length(F))
  r <- pad_to(Lc, n) - pad_to(F, n)
  kept <- seq_len(min(N - p, n))
  list(finite = sqrt(sum(rb^2) + sum(r[kept]^2)),
       truncation = sqrt(sum(r[-kept]^2)),
       scale = sqrt(sum(F^2) + sum(values^2)))
}

cheb_solve <- function(A, b, args) {
  P <- A
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
  P <- A
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
      rb <- vapply(seq_len(p), function(j) sum(cheb_boundary_row(bc$at[j], bc$derivative[j], N) * c),
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
