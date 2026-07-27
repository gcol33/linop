## Section 6. A status says whether a check passed; evidence says what kind of
## argument produced it. The two are reported separately because the classes are
## not comparable on one axis.

build_certificate <- function(df, subject, probes = NULL, values = NULL) {
  overall <- if (any(df$status == "fail")) "fail"
             else if (any(df$status %in% c("qualified", "not_checked"))) "qualified"
             else "pass"
  weak <- df$check[df$guarantee != "identity" | df$status == "not_checked"]
  structure(list(subject = subject, checks = df, overall = overall,
                 without_deterministic_bound = weak, probes = probes,
                 values = values),
            class = "linop_certificate")
}

cert_status <- function(cert, check) cert$checks$status[cert$checks$check == check]

#' @export
print.linop_certificate <- function(x, ...) {
  df <- x$checks
  cat(sprintf("%-30s %-12s %-16s %-21s %s\n", "check", "status", "source", "guarantee", "conf"))
  cat(strrep("-", 90), "\n", sep = "")
  for (i in seq_len(nrow(df))) {
    conf <- if (is.na(df$confidence[i])) "-" else format(df$confidence[i])
    cat(sprintf("%-30s %-12s %-16s %-21s %s\n",
                df$check[i], df$status[i], df$source[i], df$guarantee[i], conf))
  }
  cat(strrep("-", 90), "\n", sep = "")
  cat(sprintf("%-30s %s", "overall", x$overall))
  if (length(x$without_deterministic_bound)) {
    cat(sprintf("   no deterministic bound on: %s",
                paste(x$without_deterministic_bound, collapse = ", ")))
  }
  cat("\n")
  failed <- df[df$status == "fail", ]
  if (nrow(failed)) {
    cat("\nfailures:\n")
    for (i in seq_len(nrow(failed))) cat(sprintf("  %s: %s\n", failed$check[i], failed$detail[i]))
  }
  invisible(x)
}

## Collects rows in the order they are added, so a certificate reads top to
## bottom in the order the checks were made.
cert_rows <- function() {
  rows <- list()
  list(
    add = function(check, status, detail = "", source = "computation",
                   guarantee = "identity", confidence = 1) {
      rows[[length(rows) + 1L]] <<- data.frame(
        check = check, status = status, source = source, guarantee = guarantee,
        confidence = confidence, detail = detail, stringsAsFactors = FALSE)
    },
    collect = function() do.call(rbind, rows))
}

## ------------------------------------------------------- the arithmetic floor

## S0.6 measured the correction this exists for. The a posteriori bound keeps
## decaying past machine epsilon while the true error plateaus at about 1e-15, so
## a fully converged result certifies as fail. Plan section 6's table carries no
## roundoff term anywhere, and both the residual and the backward-error line need
## one.
##
## For a linear solve the floor is Higham's bound on the computed residual: for
## any x, fl(b - A x) differs from b - A x by at most a small multiple of
## eps * (||A|| ||x|| + ||b||), so no x, however exactly it solves the system,
## can be shown to have a residual below that level. That same quantity is the
## denominator of the Rigal-Gaches identity for the normwise relative backward
## error,
##
##     omega(x) = ||b - A x|| / (||A|| ||x|| + ||b||),
##
## so once the backward error is expressed relatively the floor is a plain
## c * eps.
##
## c counts terms in an inner product a matrix-free operator never exposes, so it
## is not knowable here. It is an argument with a default rather than a constant,
## and any line that meets its tolerance only because of it is reported as
## qualified with an estimate guarantee, never as a clean identity. The default
## comes from S0.6's measurement: a plateau of about 2.2e-15 against ||H|| of 4
## to 6 puts c between 2 and 4.
SOLVE_FLOOR_CONST <- 4

#' Certify the result of a linear solve
#'
#' @param A The operator solved against.
#' @param B The right-hand side block.
#' @param X The computed solution block.
#' @param tol The tolerance that was requested.
#' @param norm_estimate `||A||`, from `norm2()`.
#' @param iterations Iterations spent.
#' @param maxit The iteration budget.
#' @param floor_const `c` in the arithmetic floor.
#' @param norm_lower_bound A lower bound on `||A||` observed for free during the
#'   iteration. Used only when the norm estimate is not exact.
#' @return An object of class `linop_certificate`.
#' @noRd
solve_certificate <- function(A, B, X, tol, norm_estimate, iterations, maxit,
                              floor_const = SOLVE_FLOOR_CONST,
                              norm_lower_bound = 0) {
  r <- cert_rows()
  eps <- .Machine$double.eps

  ## The true residual, not the recurrence residual. One extra apply, and it is
  ## the whole difference between a certificate and a progress report: the
  ## recurrence drifts, and it drifts exactly where the answer is most converged.
  R <- B - linop_apply(A, X, "N")
  rn <- col_norms(R)
  bn <- col_norms(B)
  xn <- col_norms(X)

  ## A zero right-hand side has the exact solution 0, where a relative test is
  ## undefined; fall back to an absolute one for that column.
  scale <- ifelse(bn > 0, bn, 1)
  rel <- rn / scale

  ## Every route to ||A|| returns a lower bound, and so does the ratio the
  ## iteration observed, so the larger of the two is the tighter one and is still
  ## a lower bound. An exact norm is left alone.
  norm_is_exact <- evidence_satisfies(norm_estimate$evidence,
                                      requirement(guarantees = "identity"))
  nA <- if (norm_is_exact) norm_estimate$value
        else max(norm_estimate$value, norm_lower_bound)
  norm_guarantee <- if (norm_is_exact) "identity" else "estimate"
  norm_conf <- if (norm_is_exact) 1 else NA_real_

  denom <- nA * xn + bn
  floor_abs <- floor_const * eps * denom
  floor_rel <- floor_abs / scale

  ## ------------------------------------------------------------------ floor --
  r$add("arithmetic floor", "pass",
        sprintf("||A|| ~ %.4g by %s; c = %g, eps = %.3g, worst relative floor %.3e",
                nA, norm_estimate$method, floor_const, eps, max(floor_rel)),
        source = norm_estimate$evidence$source,
        guarantee = norm_guarantee, confidence = norm_conf)

  ## --------------------------------------------------------------- residual --
  clean <- all(rel <= tol)
  floored <- all(rel <= tol + floor_rel)
  r$add("residual",
        if (clean) "pass" else if (floored) "qualified" else "fail",
        sprintf("worst ||b - A x|| / ||b|| = %.3e against tol %.3e%s",
                max(rel), tol,
                if (clean) "" else sprintf(", floor %.3e", max(floor_rel))),
        guarantee = if (clean) "identity" else "estimate",
        confidence = if (clean) 1 else NA_real_)

  ## --------------------------------------------------------- backward error --
  ## Rigal-Gaches: omega is exactly the smallest normwise relative perturbation
  ## of A and b for which x is the exact solution. An identity, not a bound,
  ## which is why this line can carry an identity guarantee whenever ||A|| does.
  omega <- ifelse(denom > 0, rn / denom, 0)
  bclean <- all(omega <= tol)
  bfloored <- all(omega <= tol + floor_const * eps)
  r$add("backward error",
        if (bclean) "pass" else if (bfloored) "qualified" else "fail",
        sprintf("worst ||r|| / (||A|| ||x|| + ||b||) = %.3e against tol %.3e",
                max(omega), tol),
        guarantee = if (bclean) norm_guarantee else "estimate",
        confidence = if (bclean) norm_conf else NA_real_)

  ## ------------------------------------------------------------ convergence --
  ## Whether the iteration reached what was asked, and if not, why it stopped.
  ## An iteration that gave up early is not a converged one merely because it had
  ## budget left over.
  spent <- iterations >= maxit
  r$add("convergence", if (floored) "pass" else "fail",
        sprintf("%d of at most %d iterations; %s", iterations, maxit,
                if (floored) "target reached"
                else if (spent) "budget exhausted"
                else "stopped early, the residual had stopped decreasing"))

  ## --------------------------------------------------------- forward error --
  ## ||x - x*|| / ||x*|| <= kappa(A) * omega needs ||A^-1||, which no matrix-free
  ## operator supplies and no residual implies. Section 6.1's restraint: the
  ## honest answer is that this was not checked.
  r$add("forward error", "not_checked",
        "needs a condition estimate; ||A^-1|| is not available matrix-free",
        source = "computation", guarantee = "identity", confidence = NA_real_)

  build_certificate(r$collect(), subject = "solve",
                    values = list(residual = rel, backward_error = omega,
                                  norm = nA, floor = floor_rel))
}
