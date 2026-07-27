## Section 1.1. A linop behaves like a matrix under R's normal arithmetic.
##
## One matrixOps group method covers %*%, crossprod and tcrossprod. S0.1 verified
## this against a namespace-installed package: three generics, one method, and
## `.Generic` tells them apart. R >= 4.4.0 is required because crossprod and
## tcrossprod only became S3 generic then.

#' @export
matrixOps.linop <- function(x, y) {
  if (missing(y)) y <- NULL
  switch(.Generic,
    "%*%"        = linop_matmul(x, y),
    "crossprod"  = linop_matmul(make_view(as_linop_operand(x), "transpose"), y %||% x),
    "tcrossprod" = linop_matmul(x, make_view(as_linop_operand(y %||% x), "transpose")),
    stopf("unsupported matrix operation '%s' on a linop", .Generic))
}

as_linop_operand <- function(x) if (is_linop(x)) x else linop(x)

## Dispatch on the second argument, exactly as the plan specifies. Returning a
## lazy node where a matrix would have returned a vector is the fastest way to
## break the mental model, so it does not happen.
linop_matmul <- function(x, y) {
  if (is_linop(y)) {
    return(make_product(list(as_linop_operand(x), y)))
  }
  if (!is_linop(x)) stopf("linop_matmul() needs at least one linop")
  was_vector <- is.null(dim(y))
  Y <- linop_apply(x, y, "N")
  undo_block(Y, was_vector)
}

#' @export
Ops.linop <- function(e1, e2) {
  if (missing(e2)) {
    return(switch(.Generic,
      "+" = e1,
      "-" = make_scale(e1, -1),
      stopf("unary '%s' is not defined for a linop", .Generic)))
  }
  switch(.Generic,
    "+" = linop_add(e1, e2, 1),
    "-" = linop_add(e1, e2, -1),
    "*" = linop_mul(e1, e2),
    "/" = {
      if (is_linop(e2)) stopf("division by a linop is not defined; use solve()")
      linop_mul(e1, 1 / e2)
    },
    stopf("'%s' is not defined for a linop", .Generic))
}

linop_add <- function(e1, e2, sign) {
  if (!is_linop(e1) || !is_linop(e2)) {
    stopf(paste0("adding a scalar to an operator is ambiguous.\n",
                 "  For A + c*I write  A + c * linop_eye(nrow(A))"))
  }
  make_sum(list(e1, if (sign == 1) e2 else make_scale(e2, -1)))
}

linop_mul <- function(e1, e2) {
  if (is_linop(e1) && is_linop(e2)) {
    stopf(paste0("'*' between two operators is elementwise for matrices and is not defined here.\n",
                 "  For the matrix product write  A %%*%% B"))
  }
  if (is_linop(e1)) {
    if (length(e2) != 1L) stopf("only scalar multiplication is defined; got a value of length %d", length(e2))
    make_scale(e1, e2)
  } else {
    if (length(e1) != 1L) stopf("only scalar multiplication is defined; got a value of length %d", length(e1))
    make_scale(e2, e1)
  }
}

#' @export
t.linop <- function(x) make_view(x, "transpose")

#' @export
Complex.linop <- function(z) {
  if (.Generic != "Conj") stopf("'%s' is not defined for a linop", .Generic)
  make_view(z, "conjugate")
}

#' The conjugate transpose
#'
#' Distinct from [t()], which is the transpose. For a real operator the two have
#' the same action; they remain separate nodes so the printed tree keeps showing
#' what was written.
#'
#' @param x A `linop`.
#' @param ... Unused.
#' @return A `linop`.
#' @export
adjoint <- function(x, ...) UseMethod("adjoint")

#' @rdname adjoint
#' @export
adjoint.linop <- function(x, ...) make_view(x, "adjoint")

#' @rdname adjoint
#' @export
adjoint.default <- function(x, ...) Conj(t(x))
