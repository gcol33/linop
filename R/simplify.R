## Section 5.8. Local, applied at construction, never cost-model driven.
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
  new_view(A, kind)
}

make_scale <- function(A, a) {
  if (a == 1) return(A)
  if (A$node == "scale") return(make_scale(A$args$A, a * A$args$a))
  if (A$node == "diag") return(linop_diag(a * A$args$d))
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
  new_sum(flat)
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
