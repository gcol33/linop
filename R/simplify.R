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
