## What every function-space node shares once its operator is banded in the
## space's basis: the certificate of a truncation, and the search for its width.
## The Fourier node (on Z) and the Hermite and Laguerre nodes (on N) differ only
## in how the band is built, how the residual on the cut rows is formed, and
## what scale a pair is held to; those come in as arguments.
##
## The certificate rests on the identity the Jacobi section introduced: with u
## solving the section and u~ its extension by zero, (H - q) u~ agrees with
## (H_n - q) u inside the section and is nonzero outside it only on the rows the
## band reaches past the cut, which edge() computes exactly. So
## ||(H - q) u~||^2 = ||(H_n - q) u||^2 + ||r_out||^2, both terms reported, and for
## a self-adjoint H, dist(q, sigma(H)) <= sqrt(residual^2 + truncation^2), a
## distance to an eigenvalue when the spectrum is discrete.
band_residuals <- function(A, values, X, edge) {
  xn <- col_norms(X)
  sx <- ifelse(xn > 0, xn, 1)
  res <- col_norms(linop_apply(A, X, "N") - scale_cols(X, values)) / sx
  trunc <- edge(X) / sx
  list(residual = res, truncation = trunc, total = sqrt(res^2 + trunc^2),
       norm = xn)
}

band_certificate <- function(A, values, vectors, tol, iterations, maxit, floor_const,
                             requested, stop_reason, P, hn, scale, edge,
                             self_adjoint, discrete, reach, space) {
  r <- cert_rows()
  eps <- .Machine$double.eps
  n <- A$args$n
  X <- as_block(vectors)
  q <- if (is.complex(values)) as.vector(values) else as.numeric(values)
  k <- ncol(X)
  requested <- requested %||% k
  if (!k) {
    return(eigen_certificate_empty(r, requested, iterations, maxit, stop_reason,
                                   "finite section"))
  }
  flr <- floor_const * eps * hn
  want <- tol * scale(q)
  measured <- band_residuals(A, q, X, edge)
  res <- measured$residual
  trunc <- measured$truncation
  total <- measured$total
  bound <- total + flr
  met <- bound <= want

  r$add("arithmetic floor", "pass",
        sprintf("||H_n|| <= %.4g; c = %g, eps = %.3g, floor %.3e", hn, floor_const, eps, flr),
        source = "construction", guarantee = "deterministic_bound", confidence = 1)
  r$add("finite residual", cert_level(max(res / want), 1, flr / min(want)),
        sprintf("worst ||H_n u - q u|| / ||u|| = %.3e against %.3e; %d of %d iterations",
                max(res), min(want), iterations, maxit))
  if (self_adjoint) {
    cert_add_orthogonality(r, X, floor_const)
  } else {
    r$add("orthogonality", "not_checked",
          "the operator is not self-adjoint, and its eigenvectors need not be orthogonal",
          source = "computation", guarantee = "identity", confidence = NA_real_)
  }
  r$add("truncation bound", cert_level(max(trunc / want), 1, flr / min(want)),
        sprintf("worst ||r_out|| / ||u|| = %.3e against %.3e at n = %d, from %s%s",
                max(trunc), min(want), n, reach,
                if (max(trunc) <= flr) "; the tail is at the arithmetic floor" else ""),
        source = "theorem", guarantee = "deterministic_bound", confidence = 1)
  r$add("target identity", "not_checked",
        "a residual bound says the value is near some eigenvalue, not which one",
        source = "computation", guarantee = "identity", confidence = NA_real_)
  cert_add_convergence(r, met, requested, iterations, maxit,
                       if (!all(met) && all(res + flr <= want))
                         sprintf("the eigensolve met its target and the section did not; the truncation term is %.3e at n = %d",
                                 max(trunc), n)
                       else stop_reason,
                       "the subspace had stopped improving")
  if (self_adjoint) {
    r$add("forward error", if (all(met)) "pass" else "fail",
          sprintf("min_j |q - lambda_j(H)| <= %.3e = sqrt(%.3e^2 + %.3e^2) + %.3e against %.3e%s",
                  max(bound), max(res), max(trunc), flr, min(want),
                  if (discrete) "; the spectrum is discrete, so lambda_j is an eigenvalue"
                  else "; the spectrum is not established discrete, so lambda_j is a point of it, not necessarily an eigenvalue"),
          evidence = evidence("theorem", "deterministic_bound", 1))
  } else {
    r$add("forward error", "not_checked",
          sprintf("q lies in the %.3e-pseudospectrum of the operator on %s; it is not self-adjoint, so no distance to an eigenvalue follows",
                  max(total), space),
          source = "computation", guarantee = "identity", confidence = NA_real_)
  }
  build_certificate(r$collect(), subject = "finite section",
                    values = list(residual = res, truncation = trunc, bound = bound,
                                  target = want, norm = hn, floor = flr, n = n,
                                  converged = met),
                    evidence = r$collect_evidence())
}

## The width search for a node whose eigenvectors decay faster than any
## geometric rate: two measurements give the local rate, and extrapolating it is
## conservative, since the true decay accelerates.
BAND_N_START <- 16L

band_eigs <- function(A, k, control, args, n_min, space = "l^2(Z)") {
  if (!is.null(args$v0)) {
    stopf("eigs() on the operator itself does not take v0; the truncation width is what it chooses")
  }
  ctl <- section_width_control(control, BAND_N_START, n_min, free_tail = FALSE)
  width_search(A, ctl$n, ctl$n_max, k, args,
               accept = function(v) min(v$target) / 2,
               predict = function(v, history) {
                 h <- nrow(history)
                 if (h < 2L || history$eta[h - 1L] <= 0) return(2L * history$n[h])
                 rate <- (history$eta[h] / history$eta[h - 1L])^(1 / (history$n[h] - history$n[h - 1L]))
                 if (!is.finite(rate) || rate >= 1 || rate <= 0) return(2L * history$n[h])
                 target <- min(v$target) * SECTION_N_AIM
                 as.integer(ceiling(history$n[h] + log(target / history$eta[h]) / log(rate)))
               }, space = space)
}
