#' @export
print.linop <- function(x, ...) {
  cat(sprintf("<linop> %s, %s\n", fmt_dim(x$dim), x$dtype))
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
    dense = sprintf("dense[%s]", fmt_dim(op$dim)),
    sparse = sprintf("sparse[%s]", fmt_dim(op$dim)),
    fun = sprintf("fun[%s]%s", fmt_dim(op$dim),
                  if (is.null(op$args$adjoint)) " (no adjoint)" else ""),
    identity = sprintf("I[%d]", op$dim[1L]),
    diag = sprintf("diag[%d]", op$dim[1L]),
    scale = sprintf("scale by %s", format(op$args$a)),
    jacobi = sprintf("jacobi on l^2(Z), sigma_ess %s",
                     paste(sprintf("[%.6g, %.6g]", jacobi_bands(op)[, 1L],
                                   jacobi_bands(op)[, 2L]), collapse = " u ")),
    section = sprintf("section |j| <= %d", op$args$n),
    kron = "kron",
    fourier = sprintf("fourier operator on l^2(Z), order %d, multiplier of degree %d",
                      op$args$order, op$args$d),
    fourier_section = sprintf("section |k| <= %d", op$args$n),
    basis = sprintf("%s operator on l^2(N), order %d, polynomial of degree %d",
                    op$args$kind, op$args$order, op$args$degree),
    basis_section = sprintf("section 0 <= n <= %d", op$args$n),
    chebyshev = sprintf("chebyshev operator of order %d, %s -> %s, domain %s%s",
                        op$args$order, op$args$domain, op$args$range, cheb_domain_text(op),
                        if (isTRUE(op$args$self_adjoint)) ", self-adjoint in L2" else ""),
    blockdiag = sprintf("blockdiag of %d", length(op$args$ops)),
    hstack = sprintf("cbind of %d", length(op$args$ops)),
    vstack = sprintf("rbind of %d", length(op$args$ops)),
    lowrank = sprintf("lowrank[%s], rank %d", fmt_dim(op$dim), ncol(op$args$U)),
    perm = sprintf("perm[%d]", op$dim[1L]),
    power = sprintf("power %d", op$args$k),
    inverse = sprintf("inverse of a stored %s matrix", fmt_dim(op$dim)),
    op$node)
  records <- op$args$function_coefficients
  if (length(records)) {
    base <- paste0(base, "; ", paste(vapply(records, function(x)
      sprintf("%s: length %d, tail/vscale %.3e", x$name, length(x$coefficients), x$tail),
      character(1)), collapse = "; "))
  }
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
              cost = linop_cost(object), caps = do.call(rbind, caps))
  class(out) <- "summary.linop"
  out
}

#' @export
print.summary.linop <- function(x, ...) {
  cat(sprintf("<linop> %s, %s, root node '%s'\n", fmt_dim(x$dim), x$dtype, x$node))
  cat(sprintf("estimated cost: %.3g flops per column\n\n", x$cost))
  df <- x$caps
  df$evidence[is.na(df$value)] <- ""
  print(df, row.names = FALSE, right = FALSE)
  invisible(x)
}
