## Section 5.7 leaves: dense sparse fun identity diag.

## ------------------------------------------------------------------ dense ---

## Base R semantics throughout (S0.1 section 6): crossprod(M, X) is t(M) %*% X,
## the transpose and not the conjugate transpose.
dense_apply <- function(op, X, mode) {
  M <- op$args$M
  switch(mode,
    N = M %*% X,
    T = crossprod(M, X),
    C = crossprod(Conj(M), X),
    R = Conj(M) %*% X)
}

exact_caps_from_matrix <- function(M) {
  n <- nrow(M); p <- ncol(M)
  square <- n == p
  cplx <- is.complex(M)
  ## An exact check on data we already hold: computation with an identity
  ## guarantee, not a probe.
  ev <- evidence("computation", "identity", 1)
  caps <- list()
  if (square) {
    caps$symmetric <- capability(isTRUE(all.equal(M, t(M), check.attributes = FALSE)), ev)
    caps$hermitian <- if (cplx) {
      capability(isTRUE(all.equal(M, Conj(t(M)), check.attributes = FALSE)), ev)
    } else {
      caps$symmetric
    }
    offdiag <- M; diag(offdiag) <- 0
    caps$diagonal <- capability(all(offdiag == 0), ev)
    caps$triangular_upper <- capability(all(M[lower.tri(M)] == 0), ev)
    caps$triangular_lower <- capability(all(M[upper.tri(M)] == 0), ev)
  }
  caps
}

linop_dense <- function(M, properties = NULL, check = TRUE) {
  if (!is.matrix(M)) stopf("linop_dense() expects a matrix")
  if (is.integer(M) || is.logical(M)) storage.mode(M) <- "double"
  dt <- dtype_of_data(M)
  caps <- if (check) exact_caps_from_matrix(M) else list()
  caps <- utils::modifyList(caps, normalise_properties(properties))
  new_linop("dense", dim(M), dt, do.call(new_caps, caps), list(M = M),
            cost = 2 * prod(dim(M)))
}

## ----------------------------------------------------------------- sparse ---

## Real-only: Matrix 1.7.5 declares the zMatrix virtual classes but not the
## concrete ones, so a complex sparse operator cannot be constructed at all
## (S0.1 section 7). Complex structured operators go through a `fun` leaf.
sparse_apply <- function(op, X, mode) {
  M <- op$args$M
  ## The operator is real but the block need not be. Matrix has no complex
  ## sparse class, so split and recombine rather than downcast (section 5.5).
  if (is.complex(X)) {
    re <- sparse_apply(op, Re(X), mode)
    im <- sparse_apply(op, Im(X), mode)
    return(re + 1i * im)
  }
  Y <- switch(mode,
    N = M %*% X,
    T = Matrix::crossprod(M, X),
    C = Matrix::crossprod(M, X),
    R = M %*% X)
  as.matrix(Y)
}

linop_sparse <- function(M, properties = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stopf("sparse operators need the Matrix package; install it with install.packages(\"Matrix\")")
  }
  caps <- list()
  ## Symmetric storage is a contract of the class, not something probed.
  if (inherits(M, "symmetricMatrix")) {
    caps$symmetric <- capability(TRUE, ev_adapter())
    caps$hermitian <- capability(TRUE, ev_adapter())
  }
  if (inherits(M, "diagonalMatrix")) caps$diagonal <- capability(TRUE, ev_adapter())
  if (inherits(M, "triangularMatrix")) {
    up <- identical(M@uplo, "U")
    caps[[if (up) "triangular_upper" else "triangular_lower"]] <- capability(TRUE, ev_adapter())
  }
  caps <- utils::modifyList(caps, normalise_properties(properties))
  new_linop("sparse", dim(M), "double", do.call(new_caps, caps), list(M = M),
            cost = 2 * length(M@x %||% numeric(0)))
}

## -------------------------------------------------------------------- fun ---

## The author supplies the forward action and, optionally, the adjoint. The
## other two modes follow from  A^T X = conj(A^H conj(X))  and
## conj(A) X = conj(A conj(X)), so nobody writes four callbacks.
fun_apply <- function(op, X, mode) {
  f <- op$args$apply; g <- op$args$adjoint
  need_adjoint <- mode %in% c("T", "C")
  if (need_adjoint && is.null(g)) {
    stopf("this operator has no adjoint; supply adjoint = to linop() to use mode '%s'", mode)
  }
  switch(mode,
    N = f(X),
    C = g(X),
    T = Conj(g(Conj(X))),
    R = Conj(f(Conj(X))))
}

linop_fun <- function(apply, adjoint = NULL, dim, dtype = "double",
                      properties = NULL, cost = NULL) {
  if (!is.function(apply)) stopf("apply must be a function of one block argument")
  if (!is.null(adjoint) && !is.function(adjoint)) stopf("adjoint must be a function or NULL")
  if (missing(dim) || length(dim) != 2L) stopf("dim = c(nrow, ncol) is required for a callback operator")
  caps <- normalise_properties(properties)
  new_linop("fun", dim, check_dtype(dtype), do.call(new_caps, caps),
            list(apply = apply, adjoint = adjoint),
            cost = cost %||% (2 * prod(as.integer(dim))))
}

## --------------------------------------------------------------- identity ---

identity_apply <- function(op, X, mode) X

linop_identity <- function(n, dtype = "double") {
  n <- as.integer(n)
  caps <- list(
    hermitian = capability(TRUE, ev_construction()),
    symmetric = capability(TRUE, ev_construction()),
    positive_definite = capability(TRUE, ev_construction()),
    unitary = capability(TRUE, ev_construction()),
    diagonal = capability(TRUE, ev_construction()),
    full_rank = capability(TRUE, ev_construction()))
  new_linop("identity", c(n, n), check_dtype(dtype), do.call(new_caps, caps),
            list(), cost = 0)
}

## ------------------------------------------------------------------- diag ---

diag_apply <- function(op, X, mode) {
  d <- op$args$d
  if (mode %in% c("C", "R")) d <- Conj(d)
  d * X
}

linop_diag <- function(d, properties = NULL) {
  if (!is.numeric(d) && !is.complex(d)) stopf("d must be numeric or complex")
  n <- length(d)
  dt <- dtype_of_data(d)
  ev <- evidence("computation", "identity", 1)
  caps <- list(
    diagonal = capability(TRUE, ev_construction()),
    symmetric = capability(TRUE, ev_construction()),
    hermitian = capability(all(Im(d) == 0), ev),
    positive_definite = capability(all(Im(d) == 0) && all(Re(d) > 0), ev),
    full_rank = capability(all(d != 0), ev))
  caps <- utils::modifyList(caps, normalise_properties(properties))
  new_linop("diag", c(n, n), dt, do.call(new_caps, caps), list(d = d), cost = n)
}

## --------------------------------------------------------- materialisation --

dense_materialize <- function(op) op$args$M
sparse_materialize <- function(op) as.matrix(op$args$M)
identity_materialize <- function(op) diag(1, op$dim[1L])
diag_materialize <- function(op) diag(op$args$d, nrow = op$dim[1L])
fun_materialize <- function(op) {
  n <- op$dim[2L]
  E <- diag(1, n)
  if (op$dtype == "complex") storage.mode(E) <- "complex"
  linop_apply(op, E, "N")
}

register_leaf_nodes <- function() {
  linop_register_node("dense", dense_apply, dense_materialize,
                      function(op) 2 * prod(op$dim), overwrite = TRUE)
  linop_register_node("sparse", sparse_apply, sparse_materialize,
                      function(op) 2 * prod(op$dim), overwrite = TRUE)
  linop_register_node("fun", fun_apply, fun_materialize,
                      function(op) 2 * prod(op$dim), overwrite = TRUE)
  linop_register_node("identity", identity_apply, identity_materialize,
                      function(op) 0, overwrite = TRUE)
  linop_register_node("diag", diag_apply, diag_materialize,
                      function(op) op$dim[1L], overwrite = TRUE)
}
