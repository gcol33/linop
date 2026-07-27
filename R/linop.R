#' Build a linear operator
#'
#' One generic, dispatching on its first argument. A matrix or `Matrix` becomes
#' an operator wrapping that storage; a function becomes a callback operator.
#' Third-party packages implement `linop.myclass()` and nothing else.
#'
#' A `linop` behaves like a matrix: use ordinary matrix operations.
#'
#' @param x A matrix, a `Matrix`, a function implementing the forward action, or
#'   any object with a `linop()` method.
#' @param ... Passed to methods.
#' @param adjoint For the callback form, a function implementing the action of
#'   the conjugate transpose.
#' @param dim For the callback form, `c(nrow, ncol)`.
#' @param dtype `"double"` or `"complex"`.
#' @param properties A named logical vector of declared capabilities, for
#'   example `c(hermitian = TRUE, positive_definite = TRUE)`. A bare character
#'   vector is accepted as shorthand for all-`TRUE`. Declared properties are
#'   recorded with `source = "user_declaration"`.
#' @param cost Estimated flops per column.
#'
#' @return A `linop`.
#' @examples
#' A <- linop(diag(c(3, 2, 1)))
#' dim(A)
#' A %*% c(1, 1, 1)
#'
#' n <- 5
#' S <- linop(function(X) rbind(X[-1, , drop = FALSE], 0), dim = c(n, n),
#'            adjoint = function(X) rbind(0, X[-n, , drop = FALSE]))
#' verify(S)$overall
#' @export
linop <- function(x, ...) UseMethod("linop")

#' @rdname linop
#' @export
linop.linop <- function(x, ...) x

#' @rdname linop
#' @export
linop.matrix <- function(x, properties = NULL, ...) linop_dense(x, properties = properties)

#' @rdname linop
#' @export
linop.numeric <- function(x, properties = NULL, ...) {
  if (is.null(dim(x))) stopf("a bare vector is ambiguous; use linop(as.matrix(x)) for a column, or linop_diag()")
  linop_dense(x, properties = properties)
}

#' @rdname linop
#' @export
linop.function <- function(x, adjoint = NULL, dim = NULL, dtype = "double",
                           properties = NULL, cost = NULL, ...) {
  linop_fun(apply = x, adjoint = adjoint, dim = dim, dtype = dtype,
            properties = properties, cost = cost)
}

#' @rdname linop
#' @export
linop.default <- function(x, properties = NULL, ...) {
  if (inherits(x, "Matrix")) return(linop_sparse(x, properties = properties))
  stopf(paste0("no linop() method for class %s\n",
               "  Adapter authors: implement linop.%s() returning a linop, then run verify()."),
        paste(class(x), collapse = "/"), class(x)[1L])
}

#' Identity operator
#' @param n Dimension.
#' @param dtype `"double"` or `"complex"`.
#' @return A `linop`.
#' @export
linop_eye <- function(n, dtype = "double") linop_identity(n, dtype)

#' Diagonal operator
#' @param d The diagonal.
#' @param properties Declared capabilities, as in [linop()].
#' @return A `linop`.
#' @export
linop_scaling <- function(d, properties = NULL) linop_diag(d, properties = properties)
