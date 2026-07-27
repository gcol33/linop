## Section 5.11. Core carries provenance and never inspects it. Hard-coding a
## schema naming HilbertOperator, Discretisation, lift and project would leak
## Layer 3 concepts into core and invites retained object graphs, closures
## capturing large environments, serialization problems and unstable hashes.

#' Attach an opaque provenance envelope
#'
#' @param A A `linop`.
#' @param provider Package name of the provider.
#' @param payload Anything the provider needs. Core never looks inside.
#' @return A `linop`.
#' @export
set_provenance <- function(A, provider, payload) {
  if (!is_linop(A)) stopf("set_provenance() expects a linop")
  if (!is_scalar_string(provider)) stopf("provider must be a single string")
  A$provenance <- list(provider = provider, payload = payload)
  A
}

#' Remove the provenance envelope
#' @param A A `linop`.
#' @return A `linop`.
#' @export
strip_provenance <- function(A) {
  if (!is_linop(A)) stopf("strip_provenance() expects a linop")
  A$provenance <- NULL
  A
}

#' The provenance envelope, or NULL
#' @param A A `linop`.
#' @return The envelope list, or `NULL`.
#' @export
provenance <- function(A) A$provenance

#' Lift a finite vector into the provider's space
#' @param p A provenance envelope.
#' @param x A vector or matrix in the finite space.
#' @param ... Passed to methods.
#' @return Provider-defined.
#' @export
provenance_lift <- function(p, x, ...) UseMethod("provenance_lift")

#' @export
provenance_lift.default <- function(p, x, ...) {
  stopf("no provenance_lift() method registered by provider '%s'", p$provider %||% "<none>")
}

#' Refine a discretisation
#' @param p A provenance envelope.
#' @param n_new The new discretisation size.
#' @param ... Passed to methods.
#' @return Provider-defined.
#' @export
provenance_refine <- function(p, n_new, ...) UseMethod("provenance_refine")

#' @export
provenance_refine.default <- function(p, n_new, ...) {
  stopf("no provenance_refine() method registered by provider '%s'", p$provider %||% "<none>")
}

#' Residual against the original operator
#' @param p A provenance envelope.
#' @param result A finite solver result.
#' @param ... Passed to methods.
#' @return Provider-defined.
#' @export
provenance_original_residual <- function(p, result, ...) UseMethod("provenance_original_residual")

#' @export
provenance_original_residual.default <- function(p, result, ...) {
  stopf("no provenance_original_residual() method registered by provider '%s'", p$provider %||% "<none>")
}

#' One-line description of a provenance envelope
#' @param p A provenance envelope.
#' @param ... Passed to methods.
#' @return A string.
#' @export
provenance_summary <- function(p, ...) UseMethod("provenance_summary")

#' @export
provenance_summary.default <- function(p, ...) {
  sprintf("provenance from '%s' (no summary method registered)", p$provider %||% "<none>")
}
