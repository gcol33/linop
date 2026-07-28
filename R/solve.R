## Sections 1.1 and 7.1. One verb over the seven methods.
##
## solve() rather than solve_iterative(). This reverses an earlier decision in the
## plan, which had solve() error without a direct capability and point at a second
## verb: refusing was the wrong way to honour "nothing silent", because reporting
## achieves it without spending a verb, and for a linop there is usually no
## factorization available at all, so iterative is the only meaning solve() can
## have.
##
## It costs no export. solve() is a base generic, so an S3 method reaches it the
## way t(), crossprod() and %*% already do, and the API budget of section 1.1 does
## not move. test-api-budget.R asserts that.

## Every method of section 7.1. The two least-squares methods are the only ones a
## rectangular operator admits, and the only ones auto will not choose.
SOLVE_METHODS <- c("cg", "minres", "gmres", "fgmres", "bicgstab", "lsqr", "lsmr")
SOLVE_LEAST_SQUARES <- c("lsqr", "lsmr")

## The arguments solve() supplies itself, which control may therefore not carry:
## a knob with two places to set it has one place too many.
SOLVE_OWN_ARGS <- c("A", "b", "tol", "maxit", "x0", "preconditioner")

check_control <- function(control, known, what) {
  bad <- setdiff(names(control), known)
  if (!length(bad)) return(invisible(TRUE))
  stopf("%s does not take: %s\n  Through control it takes: %s",
        what, paste(bad, collapse = ", "),
        if (length(known)) paste(known, collapse = ", ") else "nothing")
}

## Resolved at call time rather than held in a list at load time, because the
## files defining these load after this one.
solver_for <- function(method) {
  switch(method,
         cg = cg_solve, minres = minres_solve, gmres = gmres_solve,
         fgmres = fgmres_solve, bicgstab = bicgstab_solve,
         lsqr = lsqr_solve, lsmr = lsmr_solve)
}

## What auto asks of a capability before it will choose a method on the caller's
## behalf, where the method itself has no requirement constant of its own. The
## same reading of the evidence fields as CG's and MINRES's: an exact check on
## data the operator already holds is an identity rather than a probe, and probes
## are excluded by the guarantee field rather than by the source list.
AUTO_DIRECT_REQUIREMENT <- requirement(
  sources = c("construction", "adapter_contract", "theorem", "computation"),
  guarantees = "identity", min_confidence = 1)

## A capability that is TRUE and whose evidence clears the minimum. Unknown is
## not false, and a value behind evidence auto does not accept is not a value
## auto may act on: that is the section 5.3 laundering case reaching dispatch.
auto_has <- function(A, name, req) {
  isTRUE(capv(A, name)) && evidence_satisfies(cape(A, name), req)
}

## Section 1.1: auto picks from declared capabilities and shape, and records the
## choice and the reason. It never picks a method whose precondition is NA or
## whose evidence is below that method's declared minimum.
##
## The order is by how much the method is allowed to assume. A direct route
## first, then the two that require something of the operator, then the fallback
## that requires nothing, which is what makes auto total rather than partial.
##
## GMRES rather than BiCGSTAB as the fallback. Both require nothing and both run
## without an adjoint, so both would answer; GMRES minimises the residual over
## the whole Krylov space, so its residual cannot rise and it has no breakdown,
## and those are the properties to prefer when the package is choosing rather
## than the caller. BiCGSTAB stays reachable by name.
auto_method <- function(A) {
  if (auto_has(A, "diagonal", AUTO_DIRECT_REQUIREMENT) &&
      auto_has(A, "full_rank", AUTO_DIRECT_REQUIREMENT)) {
    return(list(method = "direct",
                reason = "the operator is diagonal and of full rank, so the solve is exact"))
  }
  if (auto_has(A, "hermitian", CG_PD_REQUIREMENT) &&
      auto_has(A, "positive_definite", CG_PD_REQUIREMENT)) {
    return(list(method = "cg",
                reason = "the operator is hermitian positive definite"))
  }
  if (auto_has(A, "hermitian", MINRES_HERMITIAN_REQUIREMENT)) {
    return(list(method = "minres",
                reason = "the operator is hermitian, and definiteness is not established"))
  }
  list(method = "gmres",
       reason = "no capability is established, and gmres requires none")
}

## The exact solve of a diagonal operator, which is the whole of the direct half
## of section 1.1 in v0.1. The entries are recovered with one apply of a block of
## ones rather than read off a leaf, so the route works for any operator that can
## establish diagonality, not only for linop_scaling().
diag_solve <- function(A, b, tol, floor_const, norm_control) {
  was_vector <- is.null(dim(b))
  B <- as_block(b)
  n <- A$dim[1L]
  if (nrow(B) != n) {
    stopf("non-conformable: operator is %d x %d, right-hand side has %d rows",
          n, A$dim[2L], nrow(B))
  }
  d <- linop_apply(A, matrix(1, n, 1L), "N")[, 1L]
  if (any(d == 0)) {
    stopf(paste0("the operator declares diagonal and full_rank, and a diagonal entry is zero.\n",
                 "  The declaration is contradicted by the operator's own action."))
  }
  X <- B / d
  norm_est <- do.call(norm2, c(list(A = A), norm_control))
  cert <- solve_certificate(A, B, X, tol = tol, norm_estimate = norm_est,
                            iterations = 0L, maxit = 0L,
                            floor_const = floor_const)
  list(x = undo_block(X, was_vector), certificate = cert, method = "direct",
       iterations = 0L, restarts = 0L,
       converged = cert_status(cert, "residual") != "fail",
       residual = cert$values$residual, history = NULL)
}

#' Solve a linear system
#'
#' Solves `a x = b`. If the operator declares a capability that makes a direct
#' solve available, that is used; otherwise one of the seven Krylov methods of
#' the package runs. Either way the result records which happened, what residual
#' was reached, and what argument establishes it.
#'
#' @param a A `linop`.
#' @param b A vector or block of right-hand sides.
#' @param method `"auto"`, or one of `"cg"`, `"minres"`, `"gmres"`, `"fgmres"`,
#'   `"bicgstab"`, `"lsqr"`, `"lsmr"`. Naming a method is the caller asserting
#'   their own declaration and is not filtered by evidence; `"auto"` is the
#'   package choosing, and applies each method's declared evidence minimum.
#' @param preconditioner A `preconditioner`, or `NULL`.
#' @param tol Relative tolerance.
#' @param maxit Iteration budget. Defaults to the method's own default.
#' @param x0 Starting iterate, or `NULL` for zero.
#' @param details Return the full solve object rather than the solution with the
#'   certificate attached.
#' @param control A list of arguments for the chosen method, such as `restart`
#'   for GMRES or `conlim` for the bidiagonal methods. A name the method does not
#'   take is an error rather than being ignored.
#' @param ... Unused.
#' @return By default the solution, a vector for a vector `b` and a matrix
#'   otherwise, carrying the certificate in its `"certificate"` attribute.
#'   Arithmetic on the result drops that attribute, which is the documented cost
#'   of behaving like a matrix. With `details = TRUE`, the full object.
#' @export
solve.linop <- function(a, b, method = "auto", preconditioner = NULL,
                        tol = 1e-8, maxit = NULL, x0 = NULL, details = FALSE,
                        control = list(), ...) {
  if (missing(b)) {
    stopf(paste0("solve() on a linop needs a right-hand side.\n",
                 "  solve(a) would be the inverse, which is not formed here: a linop is defined by\n",
                 "  its action, and an explicit inverse is the one thing that is never available."))
  }
  if (!is_scalar_string(method)) stopf("method must be a single string")
  if (!method %in% c("auto", SOLVE_METHODS)) {
    stopf("unknown method '%s'; it is one of: auto, %s",
          method, paste(SOLVE_METHODS, collapse = ", "))
  }
  if (!is.list(control)) stopf("control must be a list")

  square <- a$dim[1L] == a$dim[2L]
  requested <- method

  if (method == "auto") {
    ## Section 1.1: base R's solve() requires square, and least squares is a
    ## different mathematical request rather than a fallback for a shape that
    ## did not fit. The caller states the intent exactly where the meaning
    ## changes, and both methods are reachable through this same verb.
    if (!square) {
      stopf(paste0("solve() needs a square operator; this one is %d x %d.\n",
                   "  A rectangular system is a least-squares problem, min ||b - a x||, which is a\n",
                   "  different request rather than the same one on a different shape. Name the method\n",
                   "  to make it: solve(a, b, method = \"lsqr\") or method = \"lsmr\"."),
            a$dim[1L], a$dim[2L])
    }
    picked <- auto_method(a)
    method <- picked$method
    reason <- picked$reason
  } else {
    if (!square && !method %in% SOLVE_LEAST_SQUARES) {
      stopf(paste0("method '%s' needs a square operator; this one is %d x %d.\n",
                   "  The least-squares methods are the ones defined for it: %s."),
            method, a$dim[1L], a$dim[2L],
            paste(SOLVE_LEAST_SQUARES, collapse = " and "))
    }
    reason <- "named by the caller, so no evidence minimum was applied"
  }

  args <- list(A = a, b = b, tol = tol, x0 = x0)
  if (!is.null(maxit)) args$maxit <- maxit
  if (method == "direct") {
    if (!is.null(preconditioner)) {
      stopf(paste0("a direct solve takes no preconditioner.\n",
                   "  The operator is diagonal, so the solve is exact and there is nothing to accelerate.\n",
                   "  Name an iterative method if the preconditioner is the point."))
    }
    if (!is.null(x0)) {
      stopf("a direct solve takes no starting iterate; the operator is diagonal and the solve is exact")
    }
    known <- c("floor_const", "norm_control")
    check_control(control, known, "a direct solve")
    fit <- do.call(diag_solve, utils::modifyList(
      list(A = a, b = b, tol = tol, floor_const = SOLVE_FLOOR_CONST,
           norm_control = list()), control))
  } else {
    fn <- solver_for(method)
    ## control carries the method's own knobs and nothing that solve() already
    ## has an argument for, so there is one place each of those is set.
    known <- setdiff(names(formals(fn)), SOLVE_OWN_ARGS)
    check_control(control, known, sprintf("method '%s'", method))
    args$preconditioner <- preconditioner
    fit <- do.call(fn, c(args, control))
  }

  fit$requested <- requested
  fit$reason <- reason
  fit$certificate$dispatch <- list(requested = requested, chosen = method,
                                   reason = reason)
  if (details) return(fit)
  structure(fit$x, certificate = fit$certificate)
}
