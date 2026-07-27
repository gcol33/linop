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

#' Treat an operator M ~ A as a preconditioner
#'
#' Applied by solving `M z = r`.
#' @param M A `linop` approximating `A`.
#' @param solver `function(M, R)` performing the solve.
#' @param ... Passed to [preconditioner()].
#' @return A `preconditioner`.
#' @export
as_preconditioner <- function(M, solver, ...) {
  if (!is_linop(M)) stopf("as_preconditioner() expects a linop")
  if (missing(solver)) stopf("as_preconditioner() needs a solver; M^-1 is not formed here")
  preconditioner(apply_inverse = function(R) solver(M, R),
                 hermitian = capv(M, "hermitian"),
                 positive_definite = capv(M, "positive_definite"),
                 dim = M$dim, ...)
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

## Section 4.3 requirement table. A fixed = FALSE preconditioner passed to CG is
## an error naming the flag, not a warning.
PRECOND_REQUIREMENTS <- list(
  cg       = c("fixed", "hermitian", "positive_definite"),
  minres   = c("fixed", "hermitian", "positive_definite"),
  gmres    = c("fixed"),
  fgmres   = character(),
  bicgstab = c("fixed"),
  lsqr     = c("fixed"),
  lsmr     = c("fixed")
)

check_preconditioner <- function(P, method) {
  if (is.null(P)) return(invisible(TRUE))
  if (!is_preconditioner(P)) stopf("preconditioner must come from preconditioner()")
  req <- PRECOND_REQUIREMENTS[[method]]
  if (is.null(req)) stopf("unknown method '%s'", method)
  for (flag in req) {
    v <- P[[flag]]
    if (!isTRUE(v)) {
      stopf(paste0("method '%s' requires a preconditioner with %s = TRUE; this one has %s = %s.\n",
                   "  %s"),
            method, flag, flag, format(v),
            if (flag == "fixed")
              "A flexible preconditioner changes between applications; use fgmres instead."
            else
              sprintf("Declare it with preconditioner(..., %s = TRUE) if you can justify it.", flag))
    }
  }
  invisible(TRUE)
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
