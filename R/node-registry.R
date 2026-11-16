## Section 5.7. The registry is the mechanism the ecosystem proposition depends
## on, so internal node types go through exactly the same door external ones do.

.linop_nodes <- new.env(parent = emptyenv())

## Apply modes, following the BLAS convention. Four rather than one, because
## section 5.4 keeps transpose, adjoint and conjugate distinct.
##   N  A X          T  A^T X          C  A^H X          R  conj(A) X
APPLY_MODES <- c("N", "T", "C", "R")

check_mode <- function(mode) {
  if (!is_scalar_string(mode) || !mode %in% APPLY_MODES) {
    stopf("mode must be one of: %s", paste(APPLY_MODES, collapse = ", "))
  }
  mode
}

## Composition of modes: applying `outer` to an operator already viewed under
## `inner`. The group is Klein four: N identity, T and R commute, C = T then R.
MODE_COMPOSE <- matrix(
  c("N", "T", "C", "R",
    "T", "N", "R", "C",
    "C", "R", "N", "T",
    "R", "C", "T", "N"),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("N", "T", "C", "R"), c("N", "T", "C", "R")))

compose_mode <- function(outer, inner) MODE_COMPOSE[outer, inner]

## Does this mode swap the operator's dimensions?
mode_transposes <- function(mode) mode %in% c("T", "C")

#' Register a node type
#'
#' @param node Node name.
#' @param apply `function(op, X, mode)` returning the action of the operator
#'   under `mode`, one of `"N"`, `"T"`, `"C"`, `"R"`.
#' @param materialize `function(op)` returning a dense matrix, or `NULL`.
#' @param cost `function(op)` returning estimated flops per column.
#' @param certify `function(ctx)` returning the certificate [eigs()] should
#'   attach to a result on this node type, or `NULL` for the eigenpair
#'   certificate every other node gets. A finite section supplies one because its
#'   result is a statement about the operator it truncates rather than about
#'   itself, which is a different set of rows.
#' @param spectrum `function(A, k, control, args)` answering [eigs()] on an
#'   operator of this type that has no matrix, or `NULL` to be refused. It is
#'   what says a node on a sequence space can be made computable at all; the
#'   node type owns how, since choosing a truncation is a statement about its own
#'   mathematics rather than about the eigensolver.
#' @param section `function(A, n)` returning the truncation of an operator of
#'   this type to the indices `-n` to `n`, or `NULL` when [finite_section()]
#'   should refuse it. The truncation is a node of its own, holding this one.
#' @param solve `function(A, b, args)` answering [solve()] on an operator of
#'   this type that has no matrix, or `NULL` to be refused.
#' @param methods Named list of `function(A, k, control, args)`, each answering
#'   [eigs()] for the `method =` value it is named by, on an operator of this
#'   type whether or not it has a matrix. A method that only a node type's own
#'   mathematics defines lives here rather than in the eigensolver.
#' @param overwrite Replace an existing registration.
#' @param converged Optional function of `(A, rn, theta, tol, floor_const)`
#'   testing the residual against the node's certificate target.
#' @return Invisibly, the node name.
#' @noRd
linop_register_node <- function(node, apply, materialize = NULL, cost = NULL,
                                certify = NULL, spectrum = NULL,
                                section = NULL, solve = NULL, methods = list(),
                                overwrite = FALSE, converged = NULL) {
  if (!is_scalar_string(node)) stopf("node must be a single string")
  if (!is.function(apply)) stopf("apply must be a function(op, X, mode)")
  if (!overwrite && !is.null(.linop_nodes[[node]])) {
    stopf("node '%s' is already registered; pass overwrite = TRUE to replace it", node)
  }
  .linop_nodes[[node]] <- list(
    node = node, apply = apply,
    materialize = materialize, certify = certify, spectrum = spectrum,
    section = section, solve = solve, methods = methods, converged = converged,
    cost = cost %||% function(op) prod(op$dim))
  invisible(node)
}

node_method_owners <- function(method) {
  nodes <- sort(ls(.linop_nodes))
  nodes[vapply(nodes, function(x) !is.null(.linop_nodes[[x]]$methods[[method]]), logical(1))]
}

get_node <- function(node) {
  h <- .linop_nodes[[node]]
  if (is.null(h)) stopf("no handler registered for node type '%s'", node)
  h
}

#' Node types currently registered
#' @return A character vector of node names.
#' @noRd
linop_nodes <- function() sort(ls(.linop_nodes))
