## Section 5.1. Immutable; every operation returns a new object.

## A dimension is two non-negative whole numbers, or Inf where the operator acts
## on a sequence space. Inf is not a sentinel and not a flag: it is the extent of
## the index set, so `rev()` for a transpose, the conformability test in a product
## and the squareness test for a hermitian operator all keep working on it without
## a special case.
##
## Finite dimensions stay integer, because every existing caller compares against
## integers and a silent change of storage mode would be a wide, invisible edit. A
## dimension with an infinite entry is double, since that is the only storage that
## can hold Inf. So the type is a fact about the operator rather than a convention.
as_dim <- function(dim) {
  if (!is.numeric(dim) || length(dim) != 2L || anyNA(dim)) {
    stopf("dim must be two non-negative numbers, each a whole number or Inf")
  }
  if (any(dim < 0)) stopf("dim must be two non-negative numbers")
  fin <- is.finite(dim)
  if (any(fin & dim != trunc(dim))) {
    stopf("a finite dimension must be a whole number; got %s", fmt_dim(dim))
  }
  if (all(fin)) as.integer(dim) else as.numeric(dim)
}

## sprintf("%d", Inf) is an error, so every message that reports a shape goes
## through this rather than through %d.
fmt_dim <- function(d) {
  paste(vapply(d, function(x) if (is.finite(x)) format(x) else "Inf", character(1)),
        collapse = " x ")
}

## Whether an operator can be computed with at all. Everything numeric -- an
## apply, a materialisation, a probe, a solve, a spectrum -- needs a block to act
## on, and a sequence space has no block.
##
## The property is recursive and the operator's own shape does not settle it: an
## operator from l^2(Z) to R^4 composed with one from R^4 to l^2(Z) is 4 x 4, and
## the intermediate vector it would have to form lives in the sequence space. So
## a finite composite of infinite factors is not computable either, and a test
## that read only the root's dim would let it through to an allocation. Caught by
## test-infinite-dim.R rather than reasoned about in advance.
##
## Computed once at construction from the children's own flags, because
## linop_apply() is the hot path and walking the expression on every apply to
## re-derive a property that cannot change would be a per-iteration tree walk.
##
## A node that holds a child without ever applying it is the exception, and
## new_linop() takes `finite` for it. A finite section holds the operator it
## truncates so that the relationship is structural and explain() walks it, and
## it never forms that operator's action on anything: it reads coefficients off
## it at construction and applies its own band. Deriving would call it infinite
## and refuse the one operator the truncation exists to produce.
expr_finite <- function(dim, args) {
  if (!all(is.finite(dim))) return(FALSE)
  kids <- c(if (is_linop(args$A)) list(args$A),
            if (is.list(args$ops)) args$ops)
  for (k in kids) if (!isTRUE(k$finite)) return(FALSE)
  TRUE
}

is_finite_dim <- function(A) isTRUE(A$finite)

## The one gate. Every numeric entry point calls it, so an infinite operator is
## refused with the verb the caller actually used and with the route out named,
## rather than failing somewhere inside an allocation.
require_finite_dim <- function(A, verb) {
  if (is_finite_dim(A)) return(invisible(TRUE))
  where <- if (all(is.finite(A$dim))) {
    sprintf("this one is %s but is built from a factor that is not", fmt_dim(A$dim))
  } else {
    sprintf("this one is %s", fmt_dim(A$dim))
  }
  stopf(paste0("%s() needs an operator with finite dimensions; %s.\n",
               "  An operator on a sequence space has no matrix and no block to act on.\n",
               "  Truncate it with finite_section() and call %s() on the result."),
        verb, where, verb)
}

new_linop <- function(node, dim, dtype, caps = new_caps(), args = list(),
                      cost = NULL, finite = NULL) {
  dim <- as_dim(dim)
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
         cost = cost, finite = finite %||% expr_finite(dim, args)),
    class = c(paste0("linop_", node), "linop"))
}

#' @export
dim.linop <- function(x) x$dim

## nrow() and ncol() are not generic in base R; both derive from dim(). S0.1
## section 5 confirms no methods are needed.

## base R's length() coerces its result to integer, and prod(c(Inf, Inf)) is not
## representable there. NA_integer_ is what an unknown length is called in R, and
## it is the honest answer: the operator has no finite number of entries.
#' @export
length.linop <- function(x) {
  if (!all(is.finite(x$dim))) return(NA_integer_)
  prod(x$dim)
}

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
  require_finite_dim(A, "apply")
  X <- as_block(X)
  expect <- if (mode_transposes(mode)) A$dim[1L] else A$dim[2L]
  if (nrow(X) != expect) {
    stopf("non-conformable: operator is %s, mode %s expects %d rows, block has %d",
          fmt_dim(A$dim), mode, expect, nrow(X))
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
