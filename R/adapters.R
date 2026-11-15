## Phase 4's adapters: storage that is not a base matrix, reached through the
## linop() generic like any other class, so they cost no export. Each package is
## a Suggests, used through :: and only when an object of its class arrives, so
## the method is inert until then.
##
## An adapter computes no capability: reading the entries to check symmetry is
## exactly what out-of-core or device storage exists to avoid, so what such an
## operator has is what the caller declares, and verify() is the check.

#' @rdname linop
#' @details A `DelayedMatrix` (from the DelayedArray package, which also covers
#'   HDF5-backed matrices from HDF5Array) becomes a callback leaf applied by
#'   DelayedArray's block-processed product, so the matrix is never read whole.
#' @export
linop.DelayedMatrix <- function(x, properties = NULL, ...) {
  if (!requireNamespace("DelayedArray", quietly = TRUE)) {
    stopf("a DelayedMatrix needs the DelayedArray package")
  }
  cplx <- identical(DelayedArray::type(x), "complex")
  linop_fun(
    apply = function(X) as.matrix(x %*% X),
    adjoint = function(X) {
      if (cplx) Conj(as.matrix(crossprod(x, Conj(X))))
      else as.matrix(crossprod(x, X))
    },
    dim = dim(x), dtype = if (cplx) "complex" else "double",
    properties = properties, cost = 2 * prod(as.numeric(dim(x))))
}

#' @rdname linop
#' @details A two-dimensional `torch_tensor` becomes a callback leaf applied by
#'   torch's matrix product on the tensor's device, with each block copied to the
#'   device and back. Only `float64` and `complex128` tensors are taken: a
#'   single-precision apply carries an error near 1e-7 that the certificates'
#'   arithmetic floor, written for double precision, would not see.
#' @export
linop.torch_tensor <- function(x, properties = NULL, ...) {
  if (!requireNamespace("torch", quietly = TRUE)) stopf("a torch_tensor needs the torch package")
  if (x$dim() != 2L) stopf("a torch_tensor operator must be two-dimensional; this one has %d dimensions",
                           x$dim())
  dt <- x$dtype
  cplx <- dt == torch::torch_cdouble()
  if (!cplx && dt != torch::torch_float64()) {
    stopf(paste0("a torch_tensor operator must be float64 or complex128; this one is %s.\n",
                 "  Convert it with x$to(dtype = torch::torch_float64())."), as.character(dt))
  }
  dev <- x$device
  to_tensor <- function(X) torch::torch_tensor(X, dtype = dt, device = dev)
  from_tensor <- function(Z) as.matrix(torch::as_array(Z$cpu()))
  xh <- if (cplx) x$t()$conj() else x$t()
  linop_fun(
    apply = function(X) from_tensor(torch::torch_mm(x, to_tensor(X))),
    adjoint = function(X) from_tensor(torch::torch_mm(xh, to_tensor(X))),
    dim = as.integer(x$shape), dtype = if (cplx) "complex" else "double",
    properties = properties, cost = 2 * prod(as.numeric(x$shape)))
}

## A big.matrix is reachable a block of rows at a time and has no product of its
## own, so the adapter forms one: N reads row blocks and multiplies each, C reads
## the same blocks and accumulates their transposed products. A file-backed
## big.matrix is memory-mapped, so this never holds more than one block of it.
## The block is sized to about BIGMATRIX_BLOCK_ENTRIES entries.
BIGMATRIX_BLOCK_ENTRIES <- 4e6

#' @rdname linop
#' @details A `big.matrix` (bigmemory, including a file-backed, memory-mapped
#'   one) becomes a callback leaf that reads it a block of rows at a time, so at
#'   most one block is in memory. bigmemory stores real numbers only.
#' @export
linop.big.matrix <- function(x, properties = NULL, ...) {
  if (!requireNamespace("bigmemory", quietly = TRUE)) stopf("a big.matrix needs the bigmemory package")
  m <- nrow(x); n <- ncol(x)
  step <- max(1L, as.integer(BIGMATRIX_BLOCK_ENTRIES %/% max(n, 1L)))
  blocks <- lapply(seq(1L, m, by = step), function(s) seq.int(s, min(m, s + step - 1L)))
  linop_fun(
    apply = function(X) {
      Y <- matrix(if (is.complex(X)) 0i else 0, m, ncol(X))
      for (r in blocks) Y[r, ] <- x[r, , drop = FALSE] %*% X
      Y
    },
    adjoint = function(X) {
      Y <- matrix(if (is.complex(X)) 0i else 0, n, ncol(X))
      for (r in blocks) Y <- Y + crossprod(x[r, , drop = FALSE], X[r, , drop = FALSE])
      Y
    },
    dim = c(m, n), dtype = "double", properties = properties, cost = 2 * m * n)
}
