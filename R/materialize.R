## Section 5.9. Three categories: user-operator materialisation is forbidden
## without an explicit request; tall solver workspace is allowed; projected
## problems are allowed. What follows is the explicit-request path.

linop_materialize <- function(A) {
  ## Recursive entry point, so it gates too: as.matrix() on a FINITE composite
  ## of infinite factors passes the gate at the verb and would allocate here.
  require_finite_dim(A, "materialise")
  h <- get_node(A$node)
  if (!is.null(h$materialize)) return(h$materialize(A))
  fun_materialize(A)
}

#' Materialise an operator as a dense matrix
#'
#' An explicit request. Guarded by size, because the whole point of the package
#' is that this normally never happens.
#'
#' @param x A `linop`.
#' @param max_entries Guard on `nrow * ncol`. Raise it deliberately.
#' @param ... Unused.
#' @return A base matrix.
#' @export
as.matrix.linop <- function(x, max_entries = 1e7, ...) {
  require_finite_dim(x, "as.matrix")
  n <- prod(as.numeric(x$dim))
  if (n > max_entries) {
    stopf(paste0("refusing to materialise a %d x %d operator (%.3g entries, about %.1f GB).\n",
                 "  Raise the guard with as.matrix(A, max_entries = ...) if you mean it."),
          x$dim[1L], x$dim[2L], n, n * 8 / 2^30)
  }
  linop_materialize(x)
}

#' @export
as.array.linop <- function(x, ...) as.matrix(x, ...)

#' Materialise an operator as a sparse matrix
#' @param A A `linop`.
#' @param max_entries Guard on `nrow * ncol`.
#' @return A `Matrix` object.
#' @export
as_sparse <- function(A, max_entries = 1e7) {
  require_finite_dim(A, "as_sparse")
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stopf("as_sparse() needs the Matrix package")
  }
  if (A$node == "sparse") return(A$args$M)
  Matrix::Matrix(as.matrix(A, max_entries = max_entries), sparse = TRUE)
}

## ------------------------------------------------------------- collapse -----

#' Collapse cheap subexpressions into dense leaves
#'
#' Cost optimisation as an inspectable transformation rather than hidden
#' behaviour, so a wrong cost model costs a suboptimal plan and never a surprise
#' allocation.
#'
#' @param A A `linop`.
#' @param expected_applies How many times the operator will be applied.
#' @param max_entries Largest leaf this is allowed to build.
#' @return A `linop` carrying a `report` attribute.
#' @export
collapse <- function(A, expected_applies = 100, max_entries = 1e6) {
  require_finite_dim(A, "collapse")
  if (!is_linop(A)) stopf("collapse() expects a linop")
  report <- list(collapsed = character(), retained = character(),
                 memory_added = 0, break_even = NA_real_)

  walk <- function(op) {
    if (op$node %in% c("dense", "sparse", "fun", "identity", "diag")) return(op)
    n <- prod(as.numeric(op$dim))
    dense_cost <- 2 * n
    if (n <= max_entries && linop_cost(op) > dense_cost) {
      saved_per_apply <- linop_cost(op) - dense_cost
      bytes <- n * if (op$dtype == "complex") 16 else 8
      if (saved_per_apply * expected_applies > 0) {
        report$collapsed <<- c(report$collapsed,
          sprintf("%s of %s", op$node, fmt_dim(op$dim)))
        report$memory_added <<- report$memory_added + bytes
        report$break_even <<- ceiling(bytes / 8 / max(saved_per_apply, 1))
        return(linop_dense(linop_materialize(op), check = FALSE))
      }
    }
    report$retained <<- c(report$retained,
      sprintf("%s %s expression", fmt_dim(op$dim), op$node))
    switch(op$node,
      sum = make_sum(lapply(op$args$ops, walk)),
      product = make_product(lapply(op$args$ops, walk)),
      scale = make_scale(walk(op$args$A), op$args$a),
      transpose = , adjoint = , conjugate = make_view(walk(op$args$A), op$node),
      op)
  }

  out <- walk(A)
  attr(out, "report") <- report
  out
}

#' @export
print.linop_collapse_report <- function(x, ...) {
  cat(sprintf("collapsed:     %s\n", if (length(x$collapsed)) paste(x$collapsed, collapse = ", ") else "nothing"))
  cat(sprintf("retained:      %s\n", if (length(x$retained)) paste(x$retained, collapse = ", ") else "nothing"))
  cat(sprintf("memory added:  %s\n", format_bytes(x$memory_added)))
  cat(sprintf("break-even:    %s applications\n", x$break_even))
  invisible(x)
}

format_bytes <- function(b) {
  if (b < 1024) return(sprintf("%d B", as.integer(b)))
  if (b < 1024^2) return(sprintf("%.1f KB", b / 1024))
  if (b < 1024^3) return(sprintf("%.1f MB", b / 1024^2))
  sprintf("%.2f GB", b / 1024^3)
}

#' Report the structure and cost of an expression
#' @param A A `linop`.
#' @return Invisibly, a data frame of nodes.
#' @export
explain <- function(A) {
  if (!is_linop(A)) stopf("explain() expects a linop")
  rows <- list()
  walk <- function(op, depth) {
    rows[[length(rows) + 1L]] <<- data.frame(
      depth = depth, node = op$node,
      nrow = op$dim[1L], ncol = op$dim[2L],
      dtype = op$dtype, cost = linop_cost(op), stringsAsFactors = FALSE)
    kids <- node_children(op)
    for (k in kids) walk(k, depth + 1L)
  }
  walk(A, 0L)
  df <- do.call(rbind, rows)
  cat(sprintf("operator %s, dtype %s, estimated %.3g flops per column\n",
              fmt_dim(A$dim), A$dtype, linop_cost(A)))
  cat(sprintf("%d nodes, depth %d\n\n", nrow(df), max(df$depth)))
  print(df, row.names = FALSE)
  invisible(df)
}

node_children <- function(op) {
  switch(op$node,
    sum = , product = op$args$ops,
    transpose = , adjoint = , conjugate = , scale = , section = list(op$args$A),
    list())
}
