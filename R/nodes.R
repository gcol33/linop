## Unary view nodes and n-ary nodes (section 5.7).

## The three view nodes share one implementation and stay three node types,
## because section 5.8 requires the printed tree to show what the author wrote.
VIEW_MODE <- c(transpose = "T", adjoint = "C", conjugate = "R")

view_apply <- function(op, X, mode) {
  linop_apply(op$args$A, X, compose_mode(mode, VIEW_MODE[[op$node]]))
}

view_materialize <- function(op) {
  M <- linop_materialize(op$args$A)
  switch(op$node, transpose = t(M), adjoint = Conj(t(M)), conjugate = Conj(M))
}

new_view <- function(A, kind) {
  d <- if (kind == "conjugate") A$dim else rev(A$dim)
  caps <- switch(kind,
    transpose = caps_transpose(A),
    adjoint   = caps_adjoint(A),
    conjugate = caps_conjugate(A))
  new_linop(kind, d, A$dtype, do.call(new_caps, caps), list(A = A),
            cost = linop_cost(A), provenance = A$provenance)
}

## ------------------------------------------------------------------ scale ---

scale_apply <- function(op, X, mode) {
  a <- op$args$a
  if (mode %in% c("C", "R")) a <- Conj(a)
  a * linop_apply(op$args$A, X, mode)
}

scale_materialize <- function(op) op$args$a * linop_materialize(op$args$A)

new_scale <- function(A, a) {
  if (length(a) != 1L) stopf("scalar multiplier must have length 1")
  dt <- promote(A$dtype, dtype_of_data(a))
  new_linop("scale", A$dim, dt, do.call(new_caps, caps_scale(A, a)),
            list(A = A, a = a), cost = linop_cost(A) + A$dim[1L],
            provenance = A$provenance)
}

## -------------------------------------------------------------------- sum ---

sum_apply <- function(op, X, mode) {
  ops <- op$args$ops
  Y <- linop_apply(ops[[1L]], X, mode)
  for (i in seq_along(ops)[-1L]) Y <- Y + linop_apply(ops[[i]], X, mode)
  Y
}

sum_materialize <- function(op) Reduce(`+`, lapply(op$args$ops, linop_materialize))

new_sum <- function(ops) {
  d <- ops[[1L]]$dim
  for (A in ops) {
    if (!identical(A$dim, d)) {
      stopf("non-conformable sum: %d x %d and %d x %d", d[1L], d[2L], A$dim[1L], A$dim[2L])
    }
  }
  dt <- promote_all(vapply(ops, function(A) A$dtype, character(1)))
  new_linop("sum", d, dt, do.call(new_caps, caps_sum(ops)), list(ops = ops),
            cost = sum(vapply(ops, linop_cost, numeric(1))),
            provenance = first_provenance(ops))
}

## ---------------------------------------------------------------- product ---

product_apply <- function(op, X, mode) {
  ops <- op$args$ops
  if (mode %in% c("N", "R")) {
    ## (A1 A2 ... Ak) X, applied right to left.
    for (i in rev(seq_along(ops))) X <- linop_apply(ops[[i]], X, mode)
  } else {
    ## (A1 ... Ak)^T = Ak^T ... A1^T, so the leftmost factor is applied last.
    for (i in seq_along(ops)) X <- linop_apply(ops[[i]], X, mode)
  }
  X
}

product_materialize <- function(op) {
  Reduce(`%*%`, lapply(op$args$ops, linop_materialize))
}

new_product <- function(ops) {
  for (i in seq_along(ops)[-1L]) {
    if (ops[[i - 1L]]$dim[2L] != ops[[i]]$dim[1L]) {
      stopf("non-conformable product: %d x %d times %d x %d",
            ops[[i - 1L]]$dim[1L], ops[[i - 1L]]$dim[2L],
            ops[[i]]$dim[1L], ops[[i]]$dim[2L])
    }
  }
  d <- c(ops[[1L]]$dim[1L], ops[[length(ops)]]$dim[2L])
  dt <- promote_all(vapply(ops, function(A) A$dtype, character(1)))
  new_linop("product", d, dt, do.call(new_caps, caps_product(ops)), list(ops = ops),
            cost = sum(vapply(ops, linop_cost, numeric(1))),
            provenance = first_provenance(ops))
}

first_provenance <- function(ops) {
  for (A in ops) if (!is.null(A$provenance)) return(A$provenance)
  NULL
}

register_composite_nodes <- function() {
  for (kind in names(VIEW_MODE)) {
    linop_register_node(kind, view_apply, view_materialize,
                        function(op) linop_cost(op$args$A), overwrite = TRUE)
  }
  linop_register_node("scale", scale_apply, scale_materialize,
                      function(op) linop_cost(op$args$A) + op$dim[1L], overwrite = TRUE)
  linop_register_node("sum", sum_apply, sum_materialize,
                      function(op) sum(vapply(op$args$ops, linop_cost, numeric(1))),
                      overwrite = TRUE)
  linop_register_node("product", product_apply, product_materialize,
                      function(op) sum(vapply(op$args$ops, linop_cost, numeric(1))),
                      overwrite = TRUE)
}
