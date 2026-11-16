## Leaves: dense sparse fun identity diag.

## ------------------------------------------------------------------ dense ---

## Base R semantics throughout: crossprod(M, X) is t(M) %*% X,
## the transpose and not the conjugate transpose.
dense_apply <- function(op, X, mode) {
  M <- op$args$M
  switch(mode,
    N = M %*% X,
    T = crossprod(M, X),
    C = cross_adjoint(M, X),
    R = conj_prod(M, X))
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
    pd <- verified_positive_definite(M)
    if (isTRUE(pd)) caps$positive_definite <- capability(TRUE, ev)
    offdiag <- M; diag(offdiag) <- 0
    caps$diagonal <- capability(all(offdiag == 0), ev)
    caps$triangular_upper <- capability(all(M[lower.tri(M)] == 0), ev)
    caps$triangular_lower <- capability(all(M[upper.tri(M)] == 0), ev)
  }
  caps
}

## Positive definiteness of a stored real symmetric matrix, proven rather than
## estimated, by Rump's approach (Verification of positive definiteness, BIT 46,
## 2006): factor M - c I in floating point and show that the factorisation's
## backward error cannot reach c. If chol() completes on the stored B, its R
## satisfies R^T R = B + E with |E| <= gamma_{n+1} |R|^T |R| (Higham, Accuracy
## and Stability of Numerical Algorithms, Theorem 10.3), and by Cauchy-Schwarz
## on the columns of R, ||E||_2 <= gamma_{n+1} sum(R^2). Forming B rounds its
## diagonal by at most u |b_ii|. So lambda_min(M) >= c - gamma_{n+1} sum(R^2) -
## u max |b_ii|, and the right side being positive is a proof, evaluated with
## the rounding of the sum itself allowed for. Exactly symmetric storage only,
## because chol() reads one triangle; no entry in the subnormal range, where
## the theorem's model of rounding does not hold.
##
## O(n^3 / 3) on a leaf whose other checks are O(n^2), so it runs only up to
## DENSE_PD_MAX_N: 0.3 ms at n = 100 and 1.2 ms at 200 with R's reference BLAS.
## It proves definiteness up to kappa near 1e12 and returns NA from about 1e14,
## where c exceeds lambda_min.
## Past that, definiteness is declared or unknown, as before. NA, never FALSE:
## a failed proof says nothing about the matrix.
DENSE_PD_MAX_N <- 200L

verified_positive_definite <- function(M) {
  n <- nrow(M)
  if (n == 0L || n > DENSE_PD_MAX_N || is.complex(M) || !all(is.finite(M)) ||
      !identical(M, t(M)) || any(diag(M) <= 0)) return(NA)
  nz <- abs(M[M != 0])
  if (any(nz < 1e-150)) return(NA)
  u <- .Machine$double.eps / 2
  gamma <- function(k) k * u / (1 - k * u)
  c0 <- 4 * gamma(n + 1) * sum(diag(M))
  B <- M - diag(c0, n)
  R <- tryCatch(suppressWarnings(chol(B)), error = function(e) NULL)
  if (is.null(R)) return(NA)
  slack <- 1 + 2 * gamma(n * n)
  if (gamma(n + 1) * sum(R * R) * slack + u * max(abs(diag(B))) < c0) TRUE else NA
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
## Complex structured operators go through a `fun` leaf.
sparse_apply <- function(op, X, mode) {
  M <- op$args$M
  ## The operator is real but the block need not be. Matrix has no complex
  ## sparse class, so split and recombine rather than downcast.
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
  ## Conj() on real storage returns the numbers it was given and allocates a
  ## duplicate to do it, so a real operator applied to a real block pays two
  ## copies of the block for the two identities above. Dropping them is the same
  ## quantity, bitwise; the complex path is unchanged.
  conj_if <- function(Z) if (is.complex(Z)) Conj(Z) else Z
  ## A callback declared real is entitled to receive real blocks, the way the
  ## sparse leaf's Matrix product is, and a complex block reaches it whenever a
  ## real operator has non-real eigenvectors. The real and imaginary parts go in
  ## as one block of twice the width, so the author sees one call, and the
  ## recombination is exact for any operator linear over the reals, including a
  ## complex one declared real by mistake.
  if (op$dtype == "double" && is.complex(X)) {
    k <- ncol(X)
    Y <- fun_apply(op, cbind(Re(X), Im(X)), mode)
    return(Y[, seq_len(k), drop = FALSE] + 1i * Y[, k + seq_len(k), drop = FALSE])
  }
  switch(mode,
    N = f(X),
    C = g(X),
    T = conj_if(g(conj_if(X))),
    R = conj_if(f(conj_if(X))))
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
  if (is.complex(d) && mode %in% c("C", "R")) d <- Conj(d)
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
