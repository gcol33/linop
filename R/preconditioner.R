## Section 4.3. Two representations are common and mixing them is a known bug
## source: M ~ A applied by solving M z = r, and P ~ A^-1 applied as z = P r.
## Internally only one direction ever exists: apply_inverse.

PRECOND_SIDES <- c("left", "right", "split")

#' Build a preconditioner
#'
#' @param apply_inverse `function(R)` returning `z = M^-1 r`, always this
#'   direction.
#' @param fixed Is the map fixed across calls? A flexible preconditioner is
#'   rejected by solvers that require a fixed one.
#' @param hermitian,positive_definite Declared properties, three-valued.
#' @param side `"left"`, `"right"` or `"split"`.
#' @param dim `c(n, n)`.
#' @param apply_inverse_adjoint `function(R)` returning `M^-H r`. Needed only by
#'   a method that runs on `A M^-1` and has to apply its adjoint, which in v0.1
#'   is LSQR and LSMR. A hermitian `M` supplies it for nothing, since
#'   `M^-H = M^-1` there, so this is for the case where the two differ.
#' @return An object of class `preconditioner`.
#' @export
preconditioner <- function(apply_inverse, fixed = TRUE, hermitian = NA,
                           positive_definite = NA, side = "left", dim = NULL,
                           apply_inverse_adjoint = NULL) {
  if (!is.function(apply_inverse)) stopf("apply_inverse must be a function of one block")
  if (!is.null(apply_inverse_adjoint) && !is.function(apply_inverse_adjoint)) {
    stopf("apply_inverse_adjoint must be a function of one block, or NULL")
  }
  if (!isTRUE(fixed) && !isFALSE(fixed)) stopf("fixed must be TRUE or FALSE")
  if (!side %in% PRECOND_SIDES) {
    stopf("side must be one of: %s", paste(PRECOND_SIDES, collapse = ", "))
  }
  structure(list(apply_inverse = apply_inverse, fixed = fixed,
                 hermitian = hermitian, positive_definite = positive_definite,
                 side = side, dim = if (is.null(dim)) NULL else as.integer(dim),
                 apply_inverse_adjoint = apply_inverse_adjoint),
            class = "preconditioner")
}

is_preconditioner <- function(x) inherits(x, "preconditioner")

#' Treat an approximate inverse as a preconditioner
#'
#' A `linop` `M ~ A` is applied by solving `M z = r`, so it needs a solver; the
#' inverse is never formed here. A `linsolve` already applies `A^-1`, so it needs
#' nothing further, and its contract (section 4.2) determines whether the
#' resulting preconditioner is fixed.
#'
#' @param x A `linop` approximating `A`, or a `linsolve`.
#' @param solver `function(M, R)` performing the solve. `linop` method only.
#' @param ... Passed to [preconditioner()].
#' @return A `preconditioner`.
#' @export
as_preconditioner <- function(x, ...) UseMethod("as_preconditioner")

#' @rdname as_preconditioner
#' @export
as_preconditioner.linop <- function(x, solver, ...) {
  if (missing(solver)) stopf("as_preconditioner() needs a solver; M^-1 is not formed here")
  ## M^-H r is the solve of M^H z = r, and the solver was declared as a function
  ## of an operator and a block, so handing it the adjoint node is within what it
  ## promised. adjoint(x) costs nothing to form.
  preconditioner(apply_inverse = function(R) solver(x, R),
                 apply_inverse_adjoint = function(R) solver(adjoint(x), R),
                 hermitian = capv(x, "hermitian"),
                 positive_definite = capv(x, "positive_definite"),
                 dim = x$dim, ...)
}

#' @rdname as_preconditioner
#' @param side `"left"`, `"right"` or `"split"`. Defaults by contract; see
#'   details.
#' @export
as_preconditioner.linsolve <- function(x, side = NULL, ...) {
  ## A solve is a fixed linear map only when its contract says so on all three
  ## axes. This mirrors verify.linsolve(), which holds a solve to repeatability
  ## under exactly this condition.
  fixed <- x$contract$fidelity %in% c("exact", "linear_approximation") &&
           x$contract$determinacy == "fixed" &&
           x$contract$randomness == "deterministic"

  ## FGMRES is the only method that accepts a flexible preconditioner and it is
  ## defined for right preconditioning only, so "right" is the default there
  ## rather than the "left" default the constructor carries.
  if (is.null(side)) side <- if (fixed) "left" else "right"

  ## hermitian and positive_definite stay NA: a linsolve declares a fidelity and
  ## a determinacy, never a symmetry. An exact solve of an SPD operator is not
  ## the same claim as a preconditioner declaring itself SPD, and CG refusing it
  ## by name is the three-valued rule working rather than a gap.
  ## The linsolve contract already carries an adjoint field and reports whether
  ## one is available, so nothing has to be inferred here.
  preconditioner(apply_inverse = x$apply_inverse,
                 apply_inverse_adjoint = x$adjoint, fixed = fixed,
                 side = side, dim = x$dim, ...)
}

#' @rdname as_preconditioner
#' @export
as_preconditioner.default <- function(x, ...) {
  stopf("as_preconditioner() expects a linop or a linsolve, not %s", class(x)[1L])
}

#' Treat an operator P ~ A^-1 as a preconditioner
#'
#' Applied by computing `z = P r`.
#' @param P A `linop` approximating `A^-1`.
#' @param ... Passed to [preconditioner()].
#' @return A `preconditioner`.
#' @export
as_preconditioner_inverse <- function(P, ...) {
  if (!is_linop(P)) stopf("as_preconditioner_inverse() expects a linop")
  preconditioner(apply_inverse = function(R) linop_apply(P, R, "N"),
                 apply_inverse_adjoint =
                   if (expr_has_adjoint(P)) function(R) linop_apply(P, R, "C"),
                 hermitian = capv(P, "hermitian"),
                 positive_definite = capv(P, "positive_definite"),
                 dim = P$dim, ...)
}

## Why the two least-squares rows below are right-only. LSQR runs on A M^-1 and
## recovers x = M^-1 y, so M^-1 is applied to vectors of the domain and M is
## n x n. A left preconditioner acts on the residual, which lives in the
## codomain, and for a rectangular A that is not the same space at all: an m x m
## map is not conformable with an n x n one. Where the two spaces do coincide it
## is still the wrong object, because min ||M^-1 (b - A x)|| is a weighted
## least-squares problem whose minimiser is a different x, so the side would be
## changing the answer rather than the path to it. Split has the same defect,
## since L^-1 A L^-H puts L^-1 on the codomain again.
LSQ_SIDE_REASON <- paste0(
  "A least-squares preconditioner acts on the domain of A: the method runs on A M^-1 and\n",
  "  recovers x = M^-1 y, so M is n x n. Left and split preconditioning act on the residual,\n",
  "  which lives in the codomain, and they change which x minimises rather than how fast it\n",
  "  is reached. Build the same map with side = 'right'.")

## Section 4.3 requirement table: the properties each method requires of a
## preconditioner, and the sides it can consume. A fixed = FALSE preconditioner
## passed to CG is an error naming the flag, not a warning.
##
## Two rows are restricted by side, and the plan's own instruction is that a row
## stays unrestricted until a check narrows it, so each restriction carries the
## reason it was narrowed on. The rows that remain unrestricted admit more than
## one standard formulation, left in the M-inner product, split through
## M = L L^H, or right, and all three are implemented for them.
PRECOND_REQUIREMENTS <- list(
  cg       = list(flags = c("fixed", "hermitian", "positive_definite"), sides = PRECOND_SIDES),
  minres   = list(flags = c("fixed", "hermitian", "positive_definite"), sides = PRECOND_SIDES),
  gmres    = list(flags = "fixed",     sides = PRECOND_SIDES),
  fgmres   = list(flags = character(), sides = "right",
                  reason = paste0("FGMRES forms its update from the right-preconditioned basis, ",
                                  "so 'right' is the only side it is defined for.")),
  bicgstab = list(flags = "fixed",     sides = PRECOND_SIDES),
  lsqr     = list(flags = "fixed",     sides = "right", reason = LSQ_SIDE_REASON),
  lsmr     = list(flags = "fixed",     sides = "right", reason = LSQ_SIDE_REASON)
)

check_preconditioner <- function(P, method) {
  if (is.null(P)) return(invisible(TRUE))
  if (!is_preconditioner(P)) stopf("preconditioner must come from preconditioner()")
  req <- PRECOND_REQUIREMENTS[[method]]
  if (is.null(req)) stopf("unknown method '%s'", method)
  for (flag in req$flags) {
    v <- P[[flag]]
    if (!isTRUE(v)) {
      stopf(paste0("method '%s' requires a preconditioner with %s = TRUE; this one has %s = %s.\n",
                   "  %s"),
            method, flag, flag, format(v),
            if (flag == "fixed")
              "A flexible preconditioner changes between applications; use fgmres, which takes a right-side preconditioner."
            else
              sprintf("Declare it with preconditioner(..., %s = TRUE) if you can justify it.", flag))
    }
  }
  if (!P$side %in% req$sides) {
    stopf("method '%s' does not accept a '%s' preconditioner; it takes %s.%s",
          method, P$side, paste(sprintf("'%s'", req$sides), collapse = " or "),
          if (is.null(req$reason)) "" else paste0("\n  ", req$reason))
  }
  invisible(TRUE)
}

## M^-1. A NULL preconditioner is the identity, and returning the block
## untouched costs nothing.
precond_applier <- function(P) {
  if (is.null(P)) return(function(R) R)
  guarded(P$apply_inverse)
}

## M^-H, for a method that runs on A M^-1 and has to apply its adjoint. A
## hermitian M supplies it for nothing, because M^-H = M^-1 there, and that is
## the common case: the declaration is what makes the shortcut legitimate rather
## than a guess, so an undeclared M is refused by name instead.
precond_adjoint_applier <- function(P, method) {
  if (is.null(P)) return(function(R) R)
  if (!is.null(P$apply_inverse_adjoint)) return(guarded(P$apply_inverse_adjoint))
  if (isTRUE(P$hermitian)) return(guarded(P$apply_inverse))
  stopf(paste0("method '%s' applies M^-H as well as M^-1, and this preconditioner supplies neither.\n",
               "  It runs on A M^-1, whose adjoint is M^-H A^H, so one apply of the iteration needs it.\n",
               "  Either pass preconditioner(..., apply_inverse_adjoint = ) or declare hermitian = TRUE,\n",
               "  which makes M^-H = M^-1 and costs nothing further."),
        method)
}

## Every solver applies M^-1 and M^-H the same way and has to defend against the
## same failure, an apply that returns the wrong shape, so the check lives once
## here rather than once per method and once per direction.
guarded <- function(f) {
  function(R) {
    Z <- as_block(f(R))
    if (!identical(dim(Z), dim(R))) {
      stopf("the preconditioner returned a %d x %d block for a %d x %d residual",
            nrow(Z), ncol(Z), nrow(R), ncol(R))
    }
    Z
  }
}

#' @export
print.preconditioner <- function(x, ...) {
  cat("<preconditioner>\n")
  cat(sprintf("  side:              %s\n", x$side))
  cat(sprintf("  fixed:             %s\n", x$fixed))
  cat(sprintf("  hermitian:         %s\n", x$hermitian))
  cat(sprintf("  positive_definite: %s\n", x$positive_definite))
  cat(sprintf("  M^-H:              %s\n",
              if (!is.null(x$apply_inverse_adjoint)) "supplied"
              else if (isTRUE(x$hermitian)) "equal to M^-1, by the hermitian declaration"
              else "absent"))
  invisible(x)
}
