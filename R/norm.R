## Section 5.3 lists norm_bound(type) with evidence among the optional operator
## capabilities. That protocol is Phase 3. The arithmetic floor S0.6 requires is
## Phase 2 and needs ||A|| now, so this is where it comes from, and it carries
## its own evidence for the reason a capability does: a norm read off the
## structure of a diagonal leaf and a norm reached by power iteration are not the
## same kind of number, and a certificate resting on one has to say which.

new_norm_estimate <- function(value, evidence, method, detail = "") {
  structure(list(value = as.numeric(value)[1L], evidence = evidence,
                 method = method, detail = detail),
            class = "linop_norm_estimate")
}

#' @export
print.linop_norm_estimate <- function(x, ...) {
  cat(sprintf("||A||_2 ~ %.6g   (%s: %s)\n", x$value, x$method, x$detail))
  cat(sprintf("  evidence: %s\n", format_evidence(x$evidence)))
  invisible(x)
}

## An expression supplies the adjoint modes only if every fun leaf inside it
## does. Nodes registered by a third party report no children, so an operator
## built from one is assumed to have an adjoint; if it does not, linop_apply()
## raises its own message naming the missing callback.
expr_has_adjoint <- function(A) {
  if (A$node == "fun") return(!is.null(A$args$adjoint))
  kids <- node_children(A)
  if (!length(kids)) return(TRUE)
  all(vapply(kids, expr_has_adjoint, logical(1)))
}

#' The spectral norm of an operator, with the evidence for it
#'
#' Exact where the structure of the node gives it exactly, exact where the
#' operator is a stored matrix small enough to factor, and a power-iteration
#' estimate otherwise.
#'
#' @param A A `linop`.
#' @param tol Relative change at which the power iteration stops. The arithmetic
#'   floor of [solve_certificate()] scales linearly in the norm, so a value good
#'   to a few percent is ample and the default is loose accordingly.
#' @param maxit Power-iteration budget.
#' @param seed Seed for the starting vector, so the estimate is reproducible.
#'   The caller's random stream is restored afterwards.
#' @param exact_max Largest stored matrix, in entries, to factor exactly.
#' @param n_probe Probes used when the operator has no adjoint.
#' @return An object of class `linop_norm_estimate`.
#' @noRd
norm2 <- function(A, tol = 1e-2, maxit = 20L, seed = 1L, exact_max = 1e5,
                  n_probe = 8L) {
  if (!is_linop(A)) stopf("norm2() expects a linop")
  if (prod(as.numeric(A$dim)) == 0) {
    return(new_norm_estimate(0, ev_construction(), "structure", "empty operator"))
  }
  args <- list(tol = tol, maxit = maxit, seed = seed, exact_max = exact_max,
               n_probe = n_probe)

  ## A view or a scale inherits its child's norm and its child's evidence, so an
  ## estimate underneath stays visible to evidence_satisfies() at the top.
  inherit <- function(child, value, detail) {
    inner <- do.call(norm2, c(list(A = child), args))
    new_norm_estimate(value(inner$value),
                      ev_construction(depends_on = list(inner$evidence)),
                      inner$method, sprintf(detail, inner$detail))
  }

  if (A$node == "identity") {
    return(new_norm_estimate(1, ev_construction(), "structure",
                             "every singular value of the identity is 1"))
  }
  if (A$node == "diag") {
    return(new_norm_estimate(max(Mod(A$args$d)), evidence("computation", "identity", 1),
                             "structure", "largest |d_i|"))
  }
  if (A$node == "scale") {
    a <- Mod(A$args$a)
    return(inherit(A$args$A, function(v) a * v, "|a| times [%s]"))
  }
  if (A$node %in% c("transpose", "adjoint", "conjugate")) {
    ## ||A^T||_2 = ||A^H||_2 = ||conj(A)||_2 = ||A||_2, since all four share a
    ## singular value set.
    return(inherit(A$args$A, function(v) v,
                   paste0(A$node, " preserves the spectral norm; [%s]")))
  }

  if (A$node %in% c("dense", "sparse") && prod(as.numeric(A$dim)) <= exact_max) {
    M <- as.matrix(linop_materialize(A))
    s <- svd(M, nu = 0L, nv = 0L)$d
    return(new_norm_estimate(if (length(s)) s[1L] else 0,
                             evidence("computation", "identity", 1), "svd",
                             sprintf("largest singular value of the stored %d x %d matrix",
                                     A$dim[1L], A$dim[2L])))
  }

  if (expr_has_adjoint(A)) norm2_power(A, tol, maxit, seed)
  else norm2_probe(A, n_probe, seed)
}

## Power iteration on A^H A. For a unit v, ||A v|| <= ||A||_2 for every v, so
## each iterate is a lower bound and the sequence climbs towards the norm. It is
## recorded as an estimate rather than a deterministic_bound: "bound" on a norm
## reads as an upper bound, and this is the other side.
norm2_power <- function(A, tol, maxit, seed) {
  n <- A$dim[2L]
  v <- with_preserved_seed(seed, {
    z <- stats::rnorm(n)
    if (A$dtype == "complex") z <- z + 1i * stats::rnorm(n)
    matrix(z, n, 1L)
  })
  nv <- col_norms(v)
  if (nv == 0) nv <- 1
  v <- v / nv

  sigma <- 0
  gap <- NA_real_
  used <- 0L
  for (it in seq_len(maxit)) {
    used <- it
    w <- linop_apply(A, v, "N")
    s_new <- col_norms(w)
    gap <- if (s_new > 0) abs(s_new - sigma) / s_new else 0
    sigma <- s_new
    u <- linop_apply(A, w, "C")
    nu <- col_norms(u)
    if (nu == 0 || gap <= tol) break
    v <- u / nu
  }
  new_norm_estimate(sigma, evidence("computation", "estimate", NA_real_), "power",
                    sprintf("power iteration on A^H A, %d steps, last relative change %.2e",
                            used, gap))
}

## Without an adjoint there is no A^H A to iterate on, so the estimate is the
## largest ||A v|| / ||v|| over random probes. Still a lower bound, and a much
## weaker one, which is what the detail says.
norm2_probe <- function(A, n_probe, seed) {
  n <- A$dim[2L]
  V <- with_preserved_seed(seed, {
    z <- matrix(stats::rnorm(n * n_probe), n, n_probe)
    if (A$dtype == "complex") z <- z + 1i * matrix(stats::rnorm(n * n_probe), n, n_probe)
    z
  })
  nv <- col_norms(V)
  nv[nv == 0] <- 1
  V <- scale_cols(V, 1 / nv)
  sigma <- max(col_norms(linop_apply(A, V, "N")))
  new_norm_estimate(sigma, evidence("computation", "estimate", NA_real_), "probe",
                    sprintf("largest ||A v|| over %d random unit probes; the operator declares no adjoint, so no A^H A iteration is available",
                            n_probe))
}
