## The structured node types: the Kronecker product, block diagonal and
## stacked operators, low rank, permutation, power and exact inverse. Each is a
## node with children, so print(), explain() and the capability lattice see the
## structure, and each applies matrix-free through its children.
##
## How the closest libraries do it: LinearMaps.jl has kron, hcat, vcat, blockdiag
## and integer powers of a map; pylops has Kronecker, HStack, VStack and BlockDiag;
## SciPy's LinearOperator has a power. All apply a Kronecker product as
## vec(B X A^T) and a block operator by splitting the input, as here.
##
## The children of an n-ary node live in args$ops and the child of a unary one in
## args$A, which is what expr_finite(), node_children() and the generic walkers
## read.

## The shape of A under an apply mode.
mode_dim <- function(A, mode) if (mode_transposes(mode)) rev(A$dim) else A$dim

## Row ranges of a block of length sum(sizes), one per part.
split_rows <- function(X, sizes) {
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  lapply(seq_along(sizes), function(i) X[seq.int(starts[i], length.out = sizes[i]), , drop = FALSE])
}

coerce_operand <- function(x, verb) {
  if (is_linop(x)) return(x)
  if (is.matrix(x) || inherits(x, "Matrix")) return(linop(x))
  stopf("%s() takes linops or matrices; got %s", verb, paste(class(x), collapse = "/"))
}

ops_dtype <- function(ops) promote_all(vapply(ops, function(A) A$dtype, character(1)))
ops_cost <- function(ops) sum(vapply(ops, linop_cost, numeric(1)))

## --------------------------------------------------------- Kronecker product --

## (A (x) B) vec(X) = vec(B X A^T), with X of size ncol(B) x ncol(A) and vec()
## column major, which is base R's kronecker() and %x% ordering. Under a mode the
## same identity holds for the mode's operators: (A (x) B)^T = A^T (x) B^T, and
## likewise for the adjoint and the conjugate. A block of k columns is two block
## applies, one of B to all k reshaped columns side by side and one of A to their
## transposes.
kron_apply <- function(op, X, mode) {
  A <- op$args$ops[[1L]]; B <- op$args$ops[[2L]]
  da <- mode_dim(A, mode); db <- mode_dim(B, mode)
  k <- ncol(X)
  W <- linop_apply(B, matrix(X, nrow = db[2L]), mode)
  Wt <- matrix(aperm(array(W, c(db[1L], da[2L], k)), c(2L, 1L, 3L)), nrow = da[2L])
  Z <- linop_apply(A, Wt, mode)
  matrix(aperm(array(Z, c(da[1L], db[1L], k)), c(2L, 1L, 3L)), nrow = db[1L] * da[1L])
}

kron_materialize <- function(op) {
  kronecker(linop_materialize(op$args$ops[[1L]]), linop_materialize(op$args$ops[[2L]]))
}

#' Kronecker product of two operators
#'
#' `A (x) B`, in the ordering of base R's [base::kronecker()]: the block in row block
#' `i` and column block `j` is `A[i, j] * B`. It is applied without forming
#' either factor's entries, as `vec(B X A^T)`, at one apply of each factor per
#' block of right-hand sides. Hermitian, symmetric, definite, unitary, diagonal,
#' triangular and real carry over when both factors have them.
#'
#' @param A,B `linop`s or matrices.
#' @return A `linop` of dimension `nrow(A) nrow(B)` by `ncol(A) ncol(B)`.
#' @examples
#' A <- matrix(c(2, 1, 0, 3), 2)
#' B <- matrix(1:6, 3)
#' K <- linop_kron(A, B)
#' max(abs(as.matrix(K) - kronecker(A, B)))
#' ## the 2-D Laplacian on a 4 x 4 grid, from the 1-D one
#' L1 <- diag(2, 4)
#' L1[cbind(1:3, 2:4)] <- -1
#' L1[cbind(2:4, 1:3)] <- -1
#' L2 <- linop_kron(L1, linop_eye(4)) + linop_kron(linop_eye(4), linop(L1))
#' dim(L2)
#' @export
linop_kron <- function(A, B) {
  A <- coerce_operand(A, "linop_kron"); B <- coerce_operand(B, "linop_kron")
  require_finite_dim(A, "linop_kron"); require_finite_dim(B, "linop_kron")
  make_kron(A, B)
}

new_kron <- function(A, B) {
  ops <- list(A, B)
  new_linop("kron", A$dim * B$dim, ops_dtype(ops), do.call(new_caps, caps_kron(ops)),
            list(ops = ops),
            cost = linop_cost(A) * B$dim[1L] + linop_cost(B) * A$dim[2L])
}

## ------------------------------------------------------------ block diagonal --

blockdiag_apply <- function(op, X, mode) {
  ops <- op$args$ops
  parts <- split_rows(X, vapply(ops, function(A) mode_dim(A, mode)[2L], numeric(1)))
  do.call(rbind, lapply(seq_along(ops), function(i) linop_apply(ops[[i]], parts[[i]], mode)))
}

blockdiag_materialize <- function(op) {
  mats <- lapply(op$args$ops, linop_materialize)
  m <- sum(vapply(mats, nrow, integer(1))); n <- sum(vapply(mats, ncol, integer(1)))
  M <- matrix(if (op$dtype == "complex") 0i else 0, m, n)
  r <- 0L; c <- 0L
  for (B in mats) {
    M[r + seq_len(nrow(B)), c + seq_len(ncol(B))] <- B
    r <- r + nrow(B); c <- c + ncol(B)
  }
  M
}

#' Block diagonal operator
#'
#' The operator acting on each block of its input with its own operand, blocks
#' in the order given. The blocks may be rectangular. Every capability the
#' blocks share carries over, and full rank does when the blocks are all tall or
#' all wide.
#'
#' @param ... `linop`s or matrices, at least one.
#' @return A `linop`.
#' @examples
#' D <- linop_blockdiag(matrix(c(2, 1, 1, 2), 2), linop_scaling(c(3, 4, 5)))
#' as.matrix(D)
#' @export
linop_blockdiag <- function(...) {
  ops <- lapply(list(...), coerce_operand, verb = "linop_blockdiag")
  if (!length(ops)) stopf("linop_blockdiag() needs at least one block")
  for (A in ops) require_finite_dim(A, "linop_blockdiag")
  make_blockdiag(ops)
}

new_blockdiag <- function(ops) {
  d <- c(sum(vapply(ops, function(A) A$dim[1L], numeric(1))),
         sum(vapply(ops, function(A) A$dim[2L], numeric(1))))
  new_linop("blockdiag", d, ops_dtype(ops), do.call(new_caps, caps_blockdiag(ops)),
            list(ops = ops), cost = ops_cost(ops))
}

## ------------------------------------------------------------------- stacks --

## [A_1 ... A_k] x = sum A_i x_i; under T or C the stack turns into a row stack
## of the transposed blocks. A row stack is the transpose of that.
hstack_apply <- function(op, X, mode) stack_apply(op, X, mode, "h")
vstack_apply <- function(op, X, mode) stack_apply(op, X, mode, "v")

stack_apply <- function(op, X, mode, direction) {
  ops <- op$args$ops
  sums <- (direction == "h") != mode_transposes(mode)
  if (sums) {
    parts <- split_rows(X, vapply(ops, function(A) mode_dim(A, mode)[2L], numeric(1)))
    Y <- linop_apply(ops[[1L]], parts[[1L]], mode)
    for (i in seq_along(ops)[-1L]) Y <- Y + linop_apply(ops[[i]], parts[[i]], mode)
    Y
  } else {
    do.call(rbind, lapply(ops, function(A) linop_apply(A, X, mode)))
  }
}

hstack_materialize <- function(op) do.call(cbind, lapply(op$args$ops, linop_materialize))
vstack_materialize <- function(op) do.call(rbind, lapply(op$args$ops, linop_materialize))

new_stack <- function(ops, direction) {
  rows <- vapply(ops, function(A) A$dim[1L], numeric(1))
  cols <- vapply(ops, function(A) A$dim[2L], numeric(1))
  if (direction == "h" && any(rows != rows[1L])) {
    stopf("cbind() needs the same number of rows in every block; got %s",
          paste(rows, collapse = ", "))
  }
  if (direction == "v" && any(cols != cols[1L])) {
    stopf("rbind() needs the same number of columns in every block; got %s",
          paste(cols, collapse = ", "))
  }
  d <- if (direction == "h") c(rows[1L], sum(cols)) else c(sum(rows), cols[1L])
  new_linop(if (direction == "h") "hstack" else "vstack", d, ops_dtype(ops),
            do.call(new_caps, caps_stack(ops, direction)), list(ops = ops),
            cost = ops_cost(ops))
}

stack_operands <- function(args, verb) {
  args <- Filter(Negate(is.null), args)
  ops <- lapply(args, function(x) {
    if (is.numeric(x) && is.null(dim(x))) x <- matrix(x, nrow = if (verb == "cbind") length(x) else 1L)
    coerce_operand(x, verb)
  })
  for (A in ops) require_finite_dim(A, verb)
  ops
}

#' @export
cbind.linop <- function(..., deparse.level = 1) {
  make_stack(stack_operands(list(...), "cbind"), "h")
}

#' @export
rbind.linop <- function(..., deparse.level = 1) {
  make_stack(stack_operands(list(...), "rbind"), "v")
}

## ----------------------------------------------------------------- low rank --

## U diag(d) V^H. Its transpose is conj(V) diag(d) U^T, its adjoint
## V diag(conj d) U^H, its conjugate conj(U) diag(conj d) V^T.
lowrank_apply <- function(op, X, mode) {
  U <- op$args$U; V <- op$args$V; d <- op$args$d
  switch(mode,
    N = U %*% (d * cross_adjoint(V, X)),
    T = conj_prod(V, d * crossprod(U, X)),
    C = V %*% (Conj(d) * cross_adjoint(U, X)),
    R = conj_prod(U, Conj(d) * crossprod(V, X)))
}

lowrank_materialize <- function(op) {
  op$args$U %*% (op$args$d * t(Conj(op$args$V)))
}

#' Low-rank operator
#'
#' `U diag(d) V^H`, applied as three thin products, `(m + n) r` work per column
#' for rank `r`. With `V` left at `U` and `d` real, it is hermitian by
#' construction, and positive semidefinite when `d` is also non-negative. A rank
#' of zero is the zero operator.
#'
#' @param U An `m x r` matrix.
#' @param V An `n x r` matrix, by default `U`.
#' @param d The `r` weights, by default all one.
#' @return A `linop`.
#' @examples
#' U <- matrix(rnorm(20), 10, 2)
#' P <- linop_lowrank(U, d = c(3, 1))
#' max(abs(as.matrix(P) - U %*% diag(c(3, 1)) %*% t(U)))
#' @export
linop_lowrank <- function(U, V = U, d = rep(1, ncol(U))) {
  same <- missing(V) || identical(U, V)
  U <- as_block(U); V <- as_block(V)
  if (!(is.numeric(U) || is.complex(U)) || !(is.numeric(V) || is.complex(V))) {
    stopf("U and V must be numeric or complex matrices")
  }
  if (ncol(U) != ncol(V)) stopf("U and V must have the same number of columns; got %d and %d",
                                ncol(U), ncol(V))
  if (length(d) != ncol(U)) stopf("d must have one weight per column of U; got %d for %d",
                                  length(d), ncol(U))
  dt <- promote_all(c(dtype_of_data(U), dtype_of_data(V), dtype_of_data(d)))
  caps <- list()
  if (same && all(Im(d) == 0) && nrow(U) == nrow(V)) {
    caps$hermitian <- capability(TRUE, ev_construction())
    if (all(Re(d) >= 0)) caps$positive_semidefinite <- capability(TRUE, ev_construction())
  }
  new_linop("lowrank", c(nrow(U), nrow(V)), dt, do.call(new_caps, caps),
            list(U = U, V = V, d = d), cost = 2 * (nrow(U) + nrow(V)) * ncol(U))
}

## -------------------------------------------------------------- permutation --

## (P X) = X[p, ]: row i of the result is row p[i] of the input. P^T = P^-1 takes
## X[order(p), ].
perm_apply <- function(op, X, mode) {
  if (mode_transposes(mode)) X[op$args$inv, , drop = FALSE] else X[op$args$p, , drop = FALSE]
}

perm_materialize <- function(op) diag(1, length(op$args$p))[op$args$p, , drop = FALSE]

#' Permutation operator
#'
#' The operator taking `x` to `x[p]`. Its transpose is its inverse, so it is
#' unitary by construction, and it is symmetric exactly when `p` is an
#' involution.
#'
#' @param p A permutation of `1:n`.
#' @return A `linop`.
#' @examples
#' P <- linop_perm(c(3, 1, 2))
#' P %*% c(10, 20, 30)
#' @export
linop_perm <- function(p) {
  if (!is.numeric(p) || anyNA(p) || any(p != trunc(p))) stopf("p must be a permutation of 1:n")
  p <- as.integer(p); n <- length(p)
  if (!identical(sort(p), seq_len(n))) stopf("p must be a permutation of 1:%d", n)
  make_perm(p)
}

new_perm <- function(p) {
  caps <- list(unitary = capability(TRUE, ev_construction()),
               full_rank = capability(TRUE, ev_construction()))
  inv <- order(p)
  if (identical(inv, p)) {
    caps$symmetric <- capability(TRUE, ev_construction())
    caps$hermitian <- capability(TRUE, ev_construction())
  }
  new_linop("perm", c(length(p), length(p)), "double", do.call(new_caps, caps),
            list(p = p, inv = inv), cost = length(p))
}

## -------------------------------------------------------------------- power --

power_apply <- function(op, X, mode) {
  for (i in seq_len(op$args$k)) X <- linop_apply(op$args$A, X, mode)
  X
}

power_materialize <- function(op) {
  M <- linop_materialize(op$args$A)
  P <- M
  for (i in seq_len(op$args$k - 1L)) P <- P %*% M
  P
}

#' Integer power of an operator
#'
#' `A^k`, applied as `k` applies of `A`. `k = 0` is the identity; a negative `k`
#' is a power of the exact inverse, so it needs an operator `solve(A)` can invert.
#' A hermitian operator to an even power is positive semidefinite, and positive
#' definite when it is also established of full rank.
#'
#' @param A A square `linop` or matrix.
#' @param k A whole number.
#' @return A `linop`.
#' @examples
#' A <- linop(matrix(c(2, 1, 1, 3), 2))
#' as.matrix(linop_power(A, 3))
#' @export
linop_power <- function(A, k) {
  A <- coerce_operand(A, "linop_power")
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k != trunc(k)) {
    stopf("k must be a single whole number")
  }
  if (A$dim[1L] != A$dim[2L]) stopf("a power needs a square operator; this one is %s",
                                    fmt_dim(A$dim))
  require_finite_dim(A, "linop_power")
  make_power(A, as.integer(k))
}

new_power <- function(A, k) {
  new_linop("power", A$dim, A$dtype, do.call(new_caps, caps_power(A, k)),
            list(A = A, k = k), cost = k * linop_cost(A))
}

## ------------------------------------------------------------ exact inverse --

## An iterative solve is not a linop, because its action depends on
## the vector it is given. An inverse is one only where it is exact and fixed
## before any vector is seen: by structure (identity, nonzero diagonal,
## permutation, unitary, a scale, a view, and products, Kronecker products,
## block diagonals and powers of those), or by a factorisation of a stored
## matrix computed once. Everything else is refused and solve(A, b) named.
inverse_apply <- function(op, X, mode) {
  Minv <- op$args$Minv
  switch(mode,
    N = Minv %*% X,
    T = crossprod(Minv, X),
    C = cross_adjoint(Minv, X),
    R = conj_prod(Minv, X))
}

inverse_materialize <- function(op) op$args$Minv

INVERSE_MAX_ENTRIES <- 4e6

exact_inverse <- function(A) {
  n <- A$dim[1L]
  if (A$dim[2L] != n) {
    stopf("solve(a) needs a square operator; this one is %s", fmt_dim(A$dim))
  }
  switch(A$node,
    identity = A,
    diag = {
      if (any(A$args$d == 0)) stopf("solve(a): the diagonal has a zero entry, so the operator is singular")
      linop_diag(1 / A$args$d)
    },
    perm = new_perm(A$args$inv),
    scale = {
      if (A$args$a == 0) stopf("solve(a): the operator is scaled by zero")
      make_scale(exact_inverse(A$args$A), 1 / A$args$a)
    },
    transpose = , adjoint = , conjugate = make_view(exact_inverse(A$args$A), A$node),
    product = make_product(rev(lapply(A$args$ops, exact_inverse))),
    kron = make_kron(exact_inverse(A$args$ops[[1L]]), exact_inverse(A$args$ops[[2L]])),
    blockdiag = make_blockdiag(lapply(A$args$ops, exact_inverse)),
    power = make_power(exact_inverse(A$args$A), A$args$k),
    inverse = A$args$A,
    dense = , sparse = factored_inverse(A),
    {
      if (isTRUE(capv(A, "unitary"))) return(make_view(A, "adjoint"))
      stopf(paste0("solve(a) forms an inverse only where it is exact: a stored matrix, or a\n",
                   "  structure whose inverse is known (diagonal, permutation, unitary, and products,\n",
                   "  Kronecker products, block diagonals and powers of those). A %s node has none;\n",
                   "  an iterative solve is not a linear operator. Use solve(a, b)."), A$node)
    })
}

## A stored matrix inverted once. The result is a fixed linear map, so it is a
## linop by that test; it costs n^2 storage and is refused past a size
## where the caller should be solving rather than inverting.
factored_inverse <- function(A) {
  n <- A$dim[1L]
  if (n * n > INVERSE_MAX_ENTRIES) {
    stopf(paste0("solve(a) on a stored %d x %d matrix would form its %.3g-entry inverse.\n",
                 "  Solve against the right-hand sides instead with solve(a, b)."),
          n, n, n * n)
  }
  M <- linop_materialize(A)
  Minv <- tryCatch(solve(M), error = function(e) {
    stopf("solve(a): the stored matrix is singular to working precision (%s)",
          conditionMessage(e))
  })
  new_linop("inverse", c(n, n), A$dtype,
            do.call(new_caps, caps_inverse(A, evidence("computation", "identity", 1))),
            list(A = A, Minv = Minv), cost = 2 * n * n)
}

## --------------------------------------------------------------- registration

register_structured_nodes <- function() {
  linop_register_node("kron", kron_apply, kron_materialize, overwrite = TRUE)
  linop_register_node("blockdiag", blockdiag_apply, blockdiag_materialize, overwrite = TRUE)
  linop_register_node("hstack", hstack_apply, hstack_materialize, overwrite = TRUE)
  linop_register_node("vstack", vstack_apply, vstack_materialize, overwrite = TRUE)
  linop_register_node("lowrank", lowrank_apply, lowrank_materialize, overwrite = TRUE)
  linop_register_node("perm", perm_apply, perm_materialize, overwrite = TRUE)
  linop_register_node("power", power_apply, power_materialize, overwrite = TRUE)
  linop_register_node("inverse", inverse_apply, inverse_materialize, overwrite = TRUE)
}
