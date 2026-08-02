## Section 5.1. Immutable; every operation returns a new object.

new_linop <- function(node, dim, dtype, caps = new_caps(), args = list(),
                      cost = NULL) {
  dim <- as.integer(dim)
  if (length(dim) != 2L || anyNA(dim) || any(dim < 0L)) {
    stopf("dim must be two non-negative integers")
  }
  check_dtype(dtype)
  if (isTRUE(cap_value(caps$hermitian)) || isTRUE(cap_value(caps$symmetric))) {
    if (dim[1L] != dim[2L]) stopf("a hermitian or symmetric operator must be square")
  }
  ## A double-typed operator is real by construction, unconditionally.
  if (dtype == "double" && is.na(cap_value(caps$real))) {
    caps$real <- capability(TRUE, ev_construction())
  }
  caps <- close_caps(caps)
  structure(
    list(node = node, dim = dim, dtype = dtype, caps = caps, args = args,
         cost = cost),
    class = c(paste0("linop_", node), "linop"))
}

#' @export
dim.linop <- function(x) x$dim

## nrow() and ncol() are not generic in base R; both derive from dim(). S0.1
## section 5 confirms no methods are needed.

#' @export
length.linop <- function(x) prod(x$dim)

#' Is this a linop?
#' @param x Any object.
#' @return `TRUE` or `FALSE`.
#' @export
is_linop <- function(x) inherits(x, "linop")

#' The scalar type of an operator
#' @param A A `linop`.
#' @return `"double"` or `"complex"`.
#' @export
dtype <- function(A) {
  if (!is_linop(A)) stopf("dtype() expects a linop")
  A$dtype
}

## ---------------------------------------------------------------- apply -----

## Tier 1. Authors implement this shape; solvers never call it directly.
linop_apply <- function(A, X, mode = "N") {
  check_mode(mode)
  X <- as_block(X)
  expect <- if (mode_transposes(mode)) A$dim[1L] else A$dim[2L]
  if (nrow(X) != expect) {
    stopf("non-conformable: operator is %d x %d, mode %s expects %d rows, block has %d",
          A$dim[1L], A$dim[2L], mode, expect, nrow(X))
  }
  h <- get_node(A$node)
  Y <- h$apply(A, X, mode)
  Y <- as_block(Y)
  out_rows <- if (mode_transposes(mode)) A$dim[2L] else A$dim[1L]
  if (nrow(Y) != out_rows || ncol(Y) != ncol(X)) {
    stopf("node '%s' returned a %d x %d block; expected %d x %d",
          A$node, nrow(Y), ncol(Y), out_rows, ncol(X))
  }
  target <- result_dtype(A$dtype, X)
  if (target == "complex" && !is.complex(Y)) storage.mode(Y) <- "complex"
  Y
}

## Tier 2, the internal ABI: Y <- alpha * op(A) X + beta * Y. Solvers only ever
## call this, so scratch reuse is uniform. Synthesised from tier 1 when a node
## does not provide it, which in v0.1 is every node.
linop_gemm <- function(A, X, Y = NULL, alpha = 1, beta = 0, mode = "N") {
  Z <- linop_apply(A, X, mode)
  if (!identical(alpha, 1)) Z <- alpha * Z
  if (!is.null(Y) && !identical(beta, 0)) {
    Y <- as_block(Y)
    if (!identical(dim(Y), dim(Z))) {
      stopf("gemm: Y is %d x %d but the product is %d x %d",
            nrow(Y), ncol(Y), nrow(Z), ncol(Z))
    }
    Z <- Z + beta * Y
  }
  Z
}

## ------------------------------------------------------------ estimated cost -

#' Estimated flops per column
#'
#' An ordering heuristic only. S0.2 measured cases where a wider block is
#' *slower* per column for memory-bound operators, which a flop count cannot
#' express, so this must never be read as a predictor of wall time.
#'
#' @param A A `linop`.
#' @return A numeric flop estimate.
#' @export
linop_cost <- function(A) {
  if (!is_linop(A)) stopf("linop_cost() expects a linop")
  if (!is.null(A$cost)) return(A$cost)
  get_node(A$node)$cost(A)
}
