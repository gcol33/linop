## Local, applied at construction, never cost-model driven.
##
## For a real operator adjoint and transpose have identical action and Conj is
## the identity, but the node is deliberately NOT rewritten here. That
## simplification happens at apply time, so print() keeps showing what the author
## wrote and a real subgraph inside a complex expression keeps the distinction
## exactly where it starts to matter.

make_view <- function(A, kind) {
  inner <- A$node
  ## involutions
  if (kind == "transpose" && inner == "transpose") return(A$args$A)
  if (kind == "adjoint"   && inner == "adjoint")   return(A$args$A)
  if (kind == "conjugate" && inner == "conjugate") return(A$args$A)
  ## mixed pairs collapse to a conjugate
  if (kind == "adjoint"   && inner == "transpose") return(make_view(A$args$A, "conjugate"))
  if (kind == "transpose" && inner == "adjoint")   return(make_view(A$args$A, "conjugate"))
  if (kind == "conjugate" && inner == "adjoint")   return(make_view(A$args$A, "transpose"))
  if (kind == "conjugate" && inner == "transpose") return(make_view(A$args$A, "adjoint"))
  ## diagonal leaves absorb the view
  if (inner == "diag") {
    d <- A$args$d
    return(switch(kind,
      transpose = A,
      adjoint   = linop_diag(Conj(d)),
      conjugate = linop_diag(Conj(d))))
  }
  if (inner == "identity") return(A)
  ## a function-space operator's adjoint is an operator of the same kind, with
  ## the adjoint domain where the space has boundaries
  if (kind == "adjoint" && inner %in% SPACE_NODES) return(space_adjoint(A))
  new_view(A, kind)
}

make_scale <- function(A, a) {
  if (a == 1) return(A)
  if (A$node == "scale") return(make_scale(A$args$A, a * A$args$a))
  if (A$node == "diag") return(linop_diag(a * A$args$d))
  if (A$node %in% SPACE_NODES) return(space_scale(A, a))
  new_scale(A, a)
}

make_sum <- function(ops) {
  ops <- Filter(function(A) !is_zero_op(A), ops)
  if (!length(ops)) stopf("empty sum")
  ## flatten nested sums so the tree stays shallow
  flat <- list()
  for (A in ops) {
    if (A$node == "sum") flat <- c(flat, A$args$ops) else flat <- c(flat, list(A))
  }
  if (length(flat) == 1L) return(flat[[1L]])
  merged <- space_sum(flat)
  if (!is.null(merged)) return(merged)
  new_sum(flat)
}

## Function-space operators are held as coefficients from which the banded form
## is built in the right range space, so a sum of two of one kind is the sum of
## their coefficients, and its order, hence its range space, is the larger one:
## the conversion ApproxFun inserts under + happens when the band is built.
## Boundary conditions belong to the problem rather than to either term, so a
## sum of Chebyshev operators needs them identical.
SPACE_NODES <- c("fourier", "chebyshev", "basis")

add_padded <- function(x, y) {
  n <- max(length(x), length(y))
  c(x, rep(0, n - length(x))) + c(y, rep(0, n - length(y)))
}

## Centred Fourier coefficients of different degrees, aligned at frequency 0.
add_centred <- function(x, y) {
  d <- max(length(x), length(y)) %/% 2L
  pad <- function(v) { e <- d - length(v) %/% 2L; c(rep(0, e), v, rep(0, e)) }
  pad(x) + pad(y)
}

space_sum <- function(ops) {
  kinds <- vapply(ops, function(A) A$node, character(1))
  if (!all(kinds %in% SPACE_NODES) || length(unique(kinds)) != 1L) return(NULL)
  kind <- kinds[1L]
  if (kind == "basis") {
    bk <- unique(vapply(ops, function(A) A$args$kind, character(1)))
    if (length(bk) != 1L) {
      stopf("a sum of operators on %s needs one basis; got %s", "l^2(N)", paste(bk, collapse = " and "))
    }
  }
  bc <- NULL
  if (kind == "chebyshev") {
    given <- Filter(Negate(is.null), lapply(ops, function(A) A$args$boundary))
    if (length(given)) {
      bc <- given[[1L]]
      if (!all(vapply(given, function(g) identical(g, bc), logical(1)))) {
        stopf(paste0("the terms of a sum of Chebyshev operators carry different boundary conditions;
",
                     "  they belong to the problem, so give them on one term or on none"))
      }
    }
  }
  if (kind == "chebyshev") {
    ## terms per order, each a Chebyshev coefficient polynomial
    p <- max(vapply(ops, function(A) length(A$args$terms), integer(1)))
    terms <- lapply(seq_len(p), function(m) Reduce(add_padded, lapply(ops, function(A)
      if (m <= length(A$args$terms)) A$args$terms[[m]] else 0)))
    return(carry_function_records(new_chebyshev(terms, 0, bc),
                                  unlist(lapply(ops, function_records), recursive = FALSE)))
  }
  coef <- Reduce(add_padded, lapply(ops, function(A) A$args$coef))
  orthogonal <- kind == "basis" && any(vapply(ops, function(A)
    identical(A$args$multiplier_basis, "orthogonal"), logical(1)))
  mult <- switch(kind,
    fourier = Reduce(add_centred, lapply(ops, function(A) A$args$vhat)),
    basis = Reduce(add_padded, lapply(ops, space_multiplier, orthogonal = orthogonal)))
  out <- space_rebuild(ops[[1L]], coef, mult,
                        if (orthogonal) "orthogonal" else "monomial")
  carry_function_records(out, unlist(lapply(ops, function_records), recursive = FALSE))
}

space_scale <- function(A, a) {
  if (A$node == "chebyshev") {
    ## a boundary condition is a statement about the solution, not a term of the
    ## operator, so it is not scaled
    return(carry_function_records(new_chebyshev(lapply(A$args$terms, function(t) a * t), 0,
                                                 A$args$boundary), scale_function_records(function_records(A), a)))
  }
  mult <- switch(A$node, fourier = A$args$vhat, basis = A$args$v)
  carry_function_records(space_rebuild(A, a * A$args$coef, a * mult), scale_function_records(function_records(A), a))
}

space_rebuild <- function(A, coef, mult, representation = A$args$multiplier_basis %||% "monomial") {
  switch(A$node,
    fourier = linop_space("fourier", derivative = coef, multiplier = mult),
    basis = new_basis_operator(A$args$kind, coef, mult, representation))
}

make_product <- function(ops) {
  flat <- list()
  for (A in ops) {
    if (A$node == "product") flat <- c(flat, A$args$ops) else flat <- c(flat, list(A))
  }
  ## drop identities
  keep <- Filter(function(A) A$node != "identity", flat)
  if (!length(keep)) return(flat[[1L]])
  ## fuse adjacent diagonals
  fused <- list(keep[[1L]])
  for (i in seq_along(keep)[-1L]) {
    last <- fused[[length(fused)]]
    this <- keep[[i]]
    if (last$node == "diag" && this$node == "diag") {
      fused[[length(fused)]] <- linop_diag(last$args$d * this$args$d)
    } else {
      fused <- c(fused, list(this))
    }
  }
  if (length(fused) == 1L) return(fused[[1L]])
  new_product(fused)
}

is_zero_op <- function(A) isTRUE(A$node == "scale") && isTRUE(A$args$a == 0)

## The structured nodes, by the same local rules: identities and diagonals
## absorb what they can, nested nodes of one kind flatten, and an exponent of 0
## or 1 is not a node.
make_kron <- function(A, B) {
  if (A$node == "identity" && B$node == "identity") {
    return(linop_identity(A$dim[1L] * B$dim[1L], promote(A$dtype, B$dtype)))
  }
  if (A$node %in% c("diag", "identity") && B$node %in% c("diag", "identity")) {
    return(linop_diag(kronecker(diag_entries(A), diag_entries(B))))
  }
  new_kron(A, B)
}

make_blockdiag <- function(ops) {
  flat <- list()
  for (A in ops) {
    if (A$node == "blockdiag") flat <- c(flat, A$args$ops) else flat <- c(flat, list(A))
  }
  if (length(flat) == 1L) return(flat[[1L]])
  if (all(vapply(flat, function(A) A$node %in% c("diag", "identity"), logical(1)))) {
    if (all(vapply(flat, function(A) A$node == "identity", logical(1)))) {
      return(linop_identity(sum(vapply(flat, function(A) A$dim[1L], numeric(1))),
                            ops_dtype(flat)))
    }
    return(linop_diag(unlist(lapply(flat, diag_entries))))
  }
  new_blockdiag(flat)
}

make_stack <- function(ops, direction) {
  kind <- if (direction == "h") "hstack" else "vstack"
  flat <- list()
  for (A in ops) {
    if (A$node == kind) flat <- c(flat, A$args$ops) else flat <- c(flat, list(A))
  }
  if (length(flat) == 1L) return(flat[[1L]])
  new_stack(flat, direction)
}

make_perm <- function(p) {
  if (identical(p, seq_along(p))) return(linop_identity(length(p)))
  new_perm(p)
}

make_power <- function(A, k) {
  if (k < 0L) return(make_power(exact_inverse(A), -k))
  if (k == 0L) return(linop_identity(A$dim[1L], A$dtype))
  if (k == 1L || A$node == "identity") return(A)
  if (A$node == "diag") return(linop_diag(A$args$d^k))
  if (A$node == "power") return(make_power(A$args$A, A$args$k * k))
  new_power(A, k)
}

diag_entries <- function(A) {
  if (A$node == "identity") rep(if (A$dtype == "complex") 1 + 0i else 1, A$dim[1L])
  else A$args$d
}
