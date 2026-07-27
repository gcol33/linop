#' @export
print.linop <- function(x, ...) {
  cat(sprintf("<linop> %d x %d, %s\n", x$dim[1L], x$dim[2L], x$dtype))
  known <- known_caps(x)
  if (length(known)) cat("  ", paste(known, collapse = ", "), "\n", sep = "")
  cat("\n")
  print_tree(x, prefix = "", last = TRUE, root = TRUE)
  invisible(x)
}

known_caps <- function(A) {
  out <- character()
  for (nm in CAPABILITY_NAMES) {
    v <- capv(A, nm)
    if (isTRUE(v)) out <- c(out, nm)
    else if (isFALSE(v)) out <- c(out, paste0("not ", nm))
  }
  out
}

node_label <- function(op) {
  base <- switch(op$node,
    dense = sprintf("dense[%d x %d]", op$dim[1L], op$dim[2L]),
    sparse = sprintf("sparse[%d x %d]", op$dim[1L], op$dim[2L]),
    fun = sprintf("fun[%d x %d]%s", op$dim[1L], op$dim[2L],
                  if (is.null(op$args$adjoint)) " (no adjoint)" else ""),
    identity = sprintf("I[%d]", op$dim[1L]),
    diag = sprintf("diag[%d]", op$dim[1L]),
    scale = sprintf("scale by %s", format(op$args$a)),
    op$node)
  base
}

print_tree <- function(op, prefix, last, root = FALSE) {
  if (root) {
    cat(node_label(op), "\n", sep = "")
  } else {
    cat(prefix, if (last) "`- " else "|- ", node_label(op), "\n", sep = "")
    prefix <- paste0(prefix, if (last) "   " else "|  ")
  }
  kids <- node_children(op)
  for (i in seq_along(kids)) {
    print_tree(kids[[i]], prefix, i == length(kids))
  }
}

#' @export
summary.linop <- function(object, ...) {
  caps <- lapply(CAPABILITY_NAMES, function(nm) {
    c1 <- cap(object, nm)
    data.frame(capability = nm,
               value = c1$value,
               evidence = format_evidence(c1$evidence),
               stringsAsFactors = FALSE)
  })
  out <- list(dim = object$dim, dtype = object$dtype, node = object$node,
              cost = linop_cost(object), caps = do.call(rbind, caps),
              provenance = object$provenance)
  class(out) <- "summary.linop"
  out
}

#' @export
print.summary.linop <- function(x, ...) {
  cat(sprintf("<linop> %d x %d, %s, root node '%s'\n", x$dim[1L], x$dim[2L], x$dtype, x$node))
  cat(sprintf("estimated cost: %.3g flops per column\n\n", x$cost))
  df <- x$caps
  df$evidence[is.na(df$value)] <- ""
  print(df, row.names = FALSE, right = FALSE)
  if (!is.null(x$provenance)) {
    cat(sprintf("\nprovenance: provider '%s'\n", x$provenance$provider))
  }
  invisible(x)
}
