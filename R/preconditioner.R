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
#' @return An object of class `preconditioner`.
#' @export
preconditioner <- function(apply_inverse, fixed = TRUE, hermitian = NA,
                           positive_definite = NA, side = "left", dim = NULL) {
  if (!is.function(apply_inverse)) stopf("apply_inverse must be a function of one block")
  if (!isTRUE(fixed) && !isFALSE(fixed)) stopf("fixed must be TRUE or FALSE")
  if (!side %in% PRECOND_SIDES) {
    stopf("side must be one of: %s", paste(PRECOND_SIDES, collapse = ", "))
  }
  structure(list(apply_inverse = apply_inverse, fixed = fixed,
                 hermitian = hermitian, positive_definite = positive_definite,
                 side = side, dim = if (is.null(dim)) NULL else as.integer(dim)),
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
  preconditioner(apply_inverse = function(R) solver(x, R),
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
  preconditioner(apply_inverse = x$apply_inverse, fixed = fixed,
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
                 hermitian = capv(P, "hermitian"),
                 positive_definite = capv(P, "positive_definite"),
                 dim = P$dim, ...)
}

## Section 4.3 requirement table: the properties each method requires of a
## preconditioner, and the sides it can consume. A fixed = FALSE preconditioner
## passed to CG is an error naming the flag, not a warning.
##
## FGMRES is the one row restricted by side: it forms its update from the
## right-preconditioned basis, so right is the only side it is defined for.
## Every other method admits more than one standard formulation, left in the
## M-inner product, split through M = L L^H, or right, so all three sides are
## accepted there.
PRECOND_REQUIREMENTS <- list(
  cg       = list(flags = c("fixed", "hermitian", "positive_definite"), sides = PRECOND_SIDES),
  minres   = list(flags = c("fixed", "hermitian", "positive_definite"), sides = PRECOND_SIDES),
  gmres    = list(flags = "fixed",     sides = PRECOND_SIDES),
  fgmres   = list(flags = character(), sides = "right"),
  bicgstab = list(flags = "fixed",     sides = PRECOND_SIDES),
  lsqr     = list(flags = "fixed",     sides = PRECOND_SIDES),
  lsmr     = list(flags = "fixed",     sides = PRECOND_SIDES)
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
          if (identical(method, "fgmres"))
            paste0("\n  FGMRES forms its update from the right-preconditioned basis, ",
                   "so 'right' is the only side it is defined for.")
          else "")
  }
  invisible(TRUE)
}

## Every solver applies M^-1 the same way and has to defend against the same
## failure, an apply_inverse() that returns the wrong shape, so the guard lives
## once here rather than once per method. A NULL preconditioner is the identity,
## and returning the block untouched costs nothing.
precond_applier <- function(P) {
  if (is.null(P)) return(function(R) R)
  function(R) {
    Z <- as_block(P$apply_inverse(R))
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
  invisible(x)
}
