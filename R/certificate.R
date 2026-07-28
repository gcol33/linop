## Section 6. A status says whether a check passed; evidence says what kind of
## argument produced it. The two are reported separately because the classes are
## not comparable on one axis.

## The guarantee kinds that are deterministic. Both of them: an identity holds
## exactly and a deterministic_bound holds always, and the summary line names the
## checks that have neither. Testing `!= "identity"` would list a theorem's bound
## among the checks with no deterministic bound behind them, which is the one
## reading of that line that cannot be right.
DETERMINISTIC_GUARANTEES <- c("identity", "deterministic_bound")

build_certificate <- function(df, subject, probes = NULL, values = NULL,
                              evidence = NULL) {
  overall <- if (any(df$status == "fail")) "fail"
             else if (any(df$status %in% c("qualified", "not_checked"))) "qualified"
             else "pass"
  weak <- df$check[!df$guarantee %in% DETERMINISTIC_GUARANTEES |
                   df$status == "not_checked"]
  structure(list(subject = subject, checks = df, overall = overall,
                 without_deterministic_bound = weak, probes = probes,
                 values = values, evidence = evidence),
            class = "linop_certificate")
}

cert_status <- function(cert, check) cert$checks$status[cert$checks$check == check]

#' @export
print.linop_certificate <- function(x, ...) {
  df <- x$checks
  ## What ran and why, where a dispatcher chose it. Section 1.1 asks that the
  ## result record which route was taken, and this is where it is recorded: the
  ## certificate is the one object a caller always gets back.
  if (!is.null(x$dispatch)) {
    cat(sprintf("solved by %s (%s)\n", x$dispatch$chosen, x$dispatch$reason))
    cat(strrep("-", 90), "\n", sep = "")
  }
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
##
## A row may carry an evidence() object instead of the three fields, in which case
## the fields are read off it rather than typed a second time. That is what a
## conditional check needs: a bound resting on a declared capability records the
## dependency, so evidence_satisfies() sees the declaration underneath the row the
## way it already sees it underneath a propagated capability (section 5.3). The
## table keeps the flat fields, since it is what the printed certificate shows.
cert_rows <- function() {
  rows <- list()
  evs <- list()
  list(
    add = function(check, status, detail = "", source = "computation",
                   guarantee = "identity", confidence = 1, evidence = NULL) {
      if (!is.null(evidence)) {
        source <- evidence$source
        guarantee <- evidence$guarantee
        confidence <- evidence$confidence
        evs[[check]] <<- evidence
      }
      rows[[length(rows) + 1L]] <<- data.frame(
        check = check, status = status, source = source, guarantee = guarantee,
        confidence = confidence, detail = detail, stringsAsFactors = FALSE)
    },
    collect = function() do.call(rbind, rows),
    collect_evidence = function() if (length(evs)) evs else NULL)
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
#' @param least_squares Read the backward error in the least-squares sense where
#'   the residual test does not hold. Costs one apply in mode `"C"`, so it is on
#'   for the methods whose problem is a minimisation and off for the ones whose
#'   problem is an equation.
#' @param stop_reason Why the iteration stopped, where the solver knows something
#'   the iteration count does not say.
#' @return An object of class `linop_certificate`.
#' @noRd
solve_certificate <- function(A, B, X, tol, norm_estimate, iterations, maxit,
                              floor_const = SOLVE_FLOOR_CONST,
                              norm_lower_bound = 0, least_squares = FALSE,
                              stop_reason = NULL) {
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

  ## Whether the compatible reading holds, per column. Everything below branches
  ## on this rather than on which method was run, so a least-squares method
  ## handed a compatible system certifies exactly as the square methods do.
  met_r <- rel <= tol + floor_rel

  ## ------------------------------------------------- the least-squares reading
  ## For an incompatible problem b - A x does not go to zero, so its size is a
  ## property of the problem rather than a measure of the solve, and testing it
  ## against tol would report a converged least-squares solution as a failure.
  ## What does go to zero is A^H r, and Stewart's perturbation makes that a
  ## backward error rather than a heuristic. With
  ##
  ##     dA = -(r r^H A) / ||r||^2,      A + dA = (I - P_r) A,
  ##
  ## the perturbed residual is r + P_r A x and (A + dA)^H (r + P_r A x) = 0
  ## identically, because (I - P_r) r = 0 and (I - P_r) P_r = 0. So x is the
  ## exact least-squares solution of min ||b - (A + dA) x||, and the size of the
  ## perturbation is ||dA||_2 = ||A^H r|| / ||r|| exactly. Relative to ||A||,
  ##
  ##     ||A^H r|| / (||A|| ||r||),
  ##
  ## which is Paige and Saunders' second stopping rule. The perturbation is
  ## exhibited rather than bounded, so this is an upper bound on the smallest one
  ## that works, and ||A|| entering through a lower bound can only enlarge it:
  ## both are the direction a certificate needs.
  ##
  ## Its floor is c eps again. The computed A^H r inherits the residual's own
  ## error through ||A|| floor_abs and adds its own c eps ||A|| ||r||, so against
  ## the denominator ||A|| ||r|| the floor is c eps + floor_abs / ||r||. The
  ## second term matters only where ||r|| has itself fallen to the residual
  ## floor, which is the regime where the compatible reading is the one in use.
  if (least_squares) {
    atn <- col_norms(linop_apply(A, R, "C"))
    optimality <- ifelse(rn > 0, atn / (nA * rn), 0)
    ls_floor <- floor_const * eps + ifelse(rn > 0, floor_abs / rn, 0)
    met_ls <- optimality <= tol + ls_floor
  } else {
    optimality <- NULL
    ls_floor <- NULL
    met_ls <- rep(FALSE, length(rel))
  }
  ## A column has been solved when either reading holds. They are not ranked:
  ## one says x solves a nearby compatible system, the other says x is the exact
  ## minimiser of a nearby least-squares problem, and which one applies is a fact
  ## about b and the range of A rather than about the iteration.
  met <- met_r | met_ls

  ## ------------------------------------------------------------------ floor --
  r$add("arithmetic floor", "pass",
        sprintf("||A|| ~ %.4g by %s; c = %g, eps = %.3g, worst relative floor %.3e%s",
                nA, norm_estimate$method, floor_const, eps, max(floor_rel),
                if (least_squares)
                  sprintf("; least-squares floor %.3e", max(ls_floor)) else ""),
        source = norm_estimate$evidence$source,
        guarantee = norm_guarantee, confidence = norm_conf)

  ## --------------------------------------------------------------- residual --
  ## A residual that cannot reach tol because b is not in the range of A is not
  ## checked rather than failed: nothing about the solve produced it. The
  ## distinction is measured here and not taken from the solver, which is why a
  ## column that misses both readings still reports fail.
  clean <- all(rel <= tol)
  floored <- all(met_r)
  incompatible <- !met_r & met_ls
  r$add("residual",
        if (clean) "pass"
        else if (floored) "qualified"
        else if (all(met)) "not_checked"
        else "fail",
        sprintf("worst ||b - A x|| / ||b|| = %.3e against tol %.3e%s",
                max(rel), tol,
                if (clean) ""
                else if (all(met) && !floored)
                  sprintf("; b is not in the range of A for %d of %d columns, and the distance to it is not a measure of the solve",
                          sum(incompatible), length(rel))
                else sprintf(", floor %.3e", max(floor_rel))),
        guarantee = if (clean) "identity" else "estimate",
        confidence = if (clean) 1 else NA_real_)

  ## --------------------------------------------------------- backward error --
  ## Rigal-Gaches: omega is exactly the smallest normwise relative perturbation
  ## of A and b for which x is the exact solution. An identity, not a bound,
  ## which is why this line can carry an identity guarantee whenever ||A|| does.
  ##
  ## Both readings land on this one row because both are backward errors. The
  ## column decides which of the two it is entitled to, and the detail says so;
  ## a least-squares column has no compatible reading and a compatible one has
  ## no use for Stewart's, since ||A^H r|| / (||A|| ||r||) stays O(1) as r goes
  ## to zero.
  omega <- ifelse(denom > 0, rn / denom, 0)
  if (least_squares) {
    bw <- ifelse(met_r, omega, optimality)
    bw_floor <- ifelse(met_r, floor_const * eps, ls_floor)
    reading <- if (all(met_r)) "||r|| / (||A|| ||x|| + ||b||)"
               else if (any(met_r)) "||r|| / (||A|| ||x|| + ||b||) and ||A^H r|| / (||A|| ||r||)"
               else "||A^H r|| / (||A|| ||r||)"
  } else {
    bw <- omega
    bw_floor <- rep(floor_const * eps, length(omega))
    reading <- "||r|| / (||A|| ||x|| + ||b||)"
  }
  bclean <- all(bw <= tol)
  bfloored <- all(bw <= tol + bw_floor)
  r$add("backward error",
        if (bclean) "pass" else if (bfloored) "qualified" else "fail",
        sprintf("worst %s = %.3e against tol %.3e", reading, max(bw), tol),
        guarantee = if (bclean) norm_guarantee else "estimate",
        confidence = if (bclean) norm_conf else NA_real_)

  ## ------------------------------------------------------------ convergence --
  ## Whether the iteration reached what was asked, and if not, why it stopped.
  ## An iteration that gave up early is not a converged one merely because it had
  ## budget left over.
  spent <- iterations >= maxit
  r$add("convergence", if (all(met)) "pass" else "fail",
        sprintf("%d of at most %d iterations; %s", iterations, maxit,
                if (all(met)) "target reached"
                else if (!is.null(stop_reason)) stop_reason
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
                    values = list(residual = rel, backward_error = bw,
                                  norm = nA, floor = floor_rel,
                                  optimality = optimality, converged = met),
                    evidence = r$collect_evidence())
}

## ------------------------------------------------------ the eigen certificate

## Section 6's table and 6.1. This is the third shape rather than a flag on the
## second. Both readings of solve_certificate() are backward errors on a linear
## system, where b is data; here theta and x are both outputs and A is the only
## datum, so the Rigal-Gaches denominator ||A|| ||x|| + ||b|| has no second term
## and only A may be perturbed.
##
## The perturbation is exhibited, as Stewart's is for the least-squares reading.
## A Ritz value is the Rayleigh quotient theta = x^H A x / (x^H x), so the
## residual r = A x - theta x satisfies x^H r = 0, and for a unit x
##
##     E = -r x^H - x r^H
##
## is hermitian with E x = -r, hence (A + E) x = theta x exactly. In the
## orthonormal basis {x, r/||r||} it is [[0, -||r||], [-||r||, 0]], so
## ||E||_2 = ||r||. The pair is exact for an operator ||r|| away from A and in the
## same class as A rather than merely near it, which is what the forward line
## below needs.
##
## Weyl: a hermitian perturbation moves each eigenvalue by at most its norm, and
## theta is an exact eigenvalue of A + E, so
##
##     min_j |theta - lambda_j(A)| <= ||r|| / ||x||.
##
## That is a deterministic bound and the first one in the package. It is one only
## as far as the hermitian declaration goes, so the row records that dependency
## and evidence_satisfies() sees a bare user_declaration underneath exactly as it
## does underneath a propagated capability. And it bounds the distance to SOME
## eigenvalue: which one is target identity, and that stays not_checked.
##
## @param A The operator. Hermitian for the forward line to be available.
## @param values The Ritz values, real.
## @param vectors The Ritz vectors, one per column.
## @param hermitian_evidence Evidence for A being hermitian, or `NULL`. The
##   forward line rests on it and records it.
## @noRd
eigen_certificate <- function(A, values, vectors, tol, norm_estimate,
                              iterations, maxit,
                              floor_const = SOLVE_FLOOR_CONST,
                              norm_lower_bound = 0,
                              hermitian_evidence = NULL, requested = NULL,
                              stop_reason = NULL, subject = "eigen") {
  r <- cert_rows()
  eps <- .Machine$double.eps
  X <- as_block(vectors)
  theta <- as.numeric(values)
  k <- ncol(X)
  ## A run that ran out of budget before its subspace was even k wide comes back
  ## with fewer pairs than were asked for, and is never padded to k with anything.
  ## The convergence line counts against what was asked rather than against what
  ## survived, or a run that produced one pair and converged it would read as
  ## complete.
  requested <- requested %||% k

  ## The true residual, recomputed here and not taken from the iteration, for the
  ## reason the solve certificate recomputes b - A x: what a recurrence believes
  ## about itself drifts, and it drifts where the answer is most converged.
  R <- linop_apply(A, X, "N") - scale_cols(X, theta)
  rn <- col_norms(R)
  xn <- col_norms(X)
  scale_x <- ifelse(xn > 0, xn, 1)

  norm_is_exact <- evidence_satisfies(norm_estimate$evidence,
                                      requirement(guarantees = "identity"))
  nA <- if (norm_is_exact) norm_estimate$value
        else max(norm_estimate$value, norm_lower_bound)
  norm_guarantee <- if (norm_is_exact) "identity" else "estimate"
  norm_conf <- if (norm_is_exact) 1 else NA_real_
  ## A zero operator has every eigenvalue at zero and no scale to measure
  ## against; an absolute test is the only one left.
  if (!is.finite(nA) || nA <= 0) nA <- 1

  ## ------------------------------------------------------------ floors ------
  ## fl(A x - theta x) differs from the exact difference by at most a small
  ## multiple of eps * (||A|| ||x|| + |theta| ||x||): Higham's bound again, with
  ## theta x in the place b held. Both readings below divide it by their own
  ## denominator.
  floor_abs <- floor_const * eps * (nA + abs(theta)) * scale_x
  bw_floor <- floor_abs / (nA * scale_x)

  bw <- rn / (nA * scale_x)
  met <- bw <= tol + bw_floor

  ## The eigenvalue equation read relatively against the term it has to cancel.
  ## Where |theta| has itself fallen to the level of eps ||A|| that ratio stops
  ## being a measure of the solve, in the way ||b - A x|| stops being one for a b
  ## outside the range of A, and it is reported the same way.
  meaningful <- abs(theta) > floor_const * eps * nA
  rel <- ifelse(meaningful, rn / (abs(theta) * scale_x), NA_real_)
  rel_floor <- ifelse(meaningful, floor_abs / (abs(theta) * scale_x), NA_real_)
  met_r <- !meaningful | rel <= tol + rel_floor

  ## ------------------------------------------------------------- floor row --
  r$add("arithmetic floor", "pass",
        sprintf("||A|| ~ %.4g by %s; c = %g, eps = %.3g, worst backward floor %.3e",
                nA, norm_estimate$method, floor_const, eps, max(bw_floor)),
        source = norm_estimate$evidence$source,
        guarantee = norm_guarantee, confidence = norm_conf)

  ## ---------------------------------------------------------- residual row --
  clean <- all(meaningful) && all(rel <= tol)
  floored <- all(meaningful) && all(met_r)
  r$add("residual",
        if (clean) "pass"
        else if (floored) "qualified"
        else if (all(met_r) && all(met)) "not_checked"
        else "fail",
        sprintf("worst ||A x - theta x|| / (|theta| ||x||) = %.3e against tol %.3e%s",
                if (any(meaningful)) max(rel[meaningful]) else NA_real_, tol,
                if (clean) ""
                else if (!all(meaningful))
                  sprintf("; |theta| is at the level of eps ||A|| for %d of %d pairs, where that ratio is not a measure of the solve",
                          sum(!meaningful), k)
                else sprintf(", floor %.3e", max(rel_floor, na.rm = TRUE))),
        guarantee = if (clean) "identity" else "estimate",
        confidence = if (clean) 1 else NA_real_)

  ## ----------------------------------------------------- orthogonality row --
  ## Measured on what is returned, not on the basis the iteration held. A dot
  ## product of length n carries a rounding error of about n eps for unit
  ## arguments, which is the level a reorthogonalised basis is entitled to reach
  ## and no better.
  G <- crossprod(Conj(X), X) - diag(1, k)
  orth <- max(Mod(G))
  orth_floor <- floor_const * eps * nrow(X)
  r$add("orthogonality",
        if (orth <= orth_floor) "pass"
        else if (orth <= sqrt(eps)) "qualified" else "fail",
        sprintf("||X^H X - I||_max = %.3e against the dot-product floor %.3e",
                orth, orth_floor))

  ## ---------------------------------------------------- backward error row --
  ## E = -r x^H - x r^H is hermitian, exhibited, and of norm exactly ||r||, and
  ## no smaller hermitian perturbation makes the pair exact. So this row is an
  ## identity whenever ||A|| is one, in the way the Rigal-Gaches row is.
  bclean <- all(bw <= tol)
  bfloored <- all(met)
  r$add("backward error",
        if (bclean) "pass" else if (bfloored) "qualified" else "fail",
        sprintf("worst ||A x - theta x|| / (||A|| ||x||) = %.3e against tol %.3e",
                max(bw), tol),
        guarantee = if (bclean) norm_guarantee else "estimate",
        confidence = if (bclean) norm_conf else NA_real_)

  ## ---------------------------------------------------- target identity row --
  ## Section 6.1, and the reason this package's certificate is not a diagnostics
  ## object. A small residual proves the pair is an approximate eigenpair. That it
  ## is the largest, the smallest or the nearest one requires inertia counts, an
  ## interval enclosure, a complete decomposition or a separation bound, and
  ## matrix-free there is generally none of those. Convergence history is not
  ## evidence for it.
  r$add("target identity", "not_checked",
        "needs inertia counts, an enclosure or a separation bound; a residual does not imply which eigenvalue was found",
        source = "computation", guarantee = "identity", confidence = NA_real_)

  ## ------------------------------------------------------- convergence row --
  spent <- iterations >= maxit
  complete <- all(met) && k >= requested
  r$add("convergence", if (complete) "pass" else "fail",
        sprintf("%d of %d requested pairs converged in %d of at most %d iterations; %s",
                sum(met), requested, iterations, maxit,
                if (complete) "target reached"
                else if (!is.null(stop_reason)) stop_reason
                else if (spent) "budget exhausted"
                else "stopped early, the subspace had stopped improving"))

  ## ----------------------------------------------------- forward error row --
  ## The one line in the package that carries a theorem rather than a
  ## measurement, and it carries the declaration it rests on with it.
  fwd_abs <- rn / scale_x
  if (is.null(hermitian_evidence)) {
    r$add("forward error", "not_checked",
          paste0("Weyl's bound needs the operator to be hermitian, and nothing establishes it; ",
                 "without that, a residual constrains no eigenvalue"),
          source = "computation", guarantee = "identity", confidence = NA_real_)
  } else {
    r$add("forward error", if (all(fwd_abs <= tol * nA)) "pass" else "fail",
          sprintf(paste0("min_j |theta - lambda_j(A)| <= %.3e, that is %.3e relative to ||A||; ",
                         "the eigenvectors are not covered, which needs a certified separation and not a Ritz gap"),
                  max(fwd_abs), max(fwd_abs) / nA),
          evidence = evidence("theorem", "deterministic_bound", 1,
                              depends_on = list(hermitian_evidence)))
  }

  build_certificate(r$collect(), subject = subject,
                    values = list(residual = rel, backward_error = bw,
                                  forward_bound = fwd_abs, orthogonality = orth,
                                  norm = nA, floor = bw_floor, converged = met),
                    evidence = r$collect_evidence())
}
