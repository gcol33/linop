## Section 5.3, capability set. `hermitian` and `symmetric` are separate because
## A = A^H and A = A^T differ for complex operators, and a merged flag sends a
## solver down the wrong path (see S0.1 section 6).

CAPABILITY_NAMES <- c(
  "hermitian", "symmetric", "positive_definite", "positive_semidefinite",
  "unitary", "diagonal", "triangular_upper", "triangular_lower",
  "full_rank", "real"
)

new_caps <- function(...) {
  caps <- list(...)
  bad <- setdiff(names(caps), CAPABILITY_NAMES)
  if (length(bad)) stopf("unknown capability: %s", paste(bad, collapse = ", "))
  out <- stats::setNames(vector("list", length(CAPABILITY_NAMES)), CAPABILITY_NAMES)
  for (nm in CAPABILITY_NAMES) {
    out[[nm]] <- if (!is.null(caps[[nm]])) caps[[nm]] else cap_unknown()
  }
  out
}

#' Read a capability off an operator
#'
#' @param A A `linop`.
#' @param name A capability name.
#' @return The `linop_capability`, or an unknown capability if absent.
#' @export
cap <- function(A, name) {
  if (!name %in% CAPABILITY_NAMES) stopf("unknown capability: %s", name)
  A$caps[[name]] %||% cap_unknown()
}

## Convenience for the propagation rules: the three-valued value only.
capv <- function(A, name) cap_value(cap(A, name))

cape <- function(A, name) cap(A, name)$evidence

## Normalise the user-facing `properties=` argument. A named logical is the
## documented form; a bare character vector is shorthand for all-TRUE.
normalise_properties <- function(properties) {
  if (is.null(properties) || !length(properties)) return(list())
  if (is.character(properties)) {
    properties <- stats::setNames(rep(TRUE, length(properties)), properties)
  }
  if (!is.logical(properties) || is.null(names(properties))) {
    stopf("properties must be a named logical vector, or a character vector shorthand")
  }
  bad <- setdiff(names(properties), CAPABILITY_NAMES)
  if (length(bad)) {
    stopf("unknown properties: %s\n  known: %s",
          paste(bad, collapse = ", "), paste(CAPABILITY_NAMES, collapse = ", "))
  }
  out <- list()
  for (nm in names(properties)) {
    v <- properties[[nm]]
    out[[nm]] <- if (is.na(v)) cap_unknown() else capability(v, ev_declared())
  }
  out
}

## A solver stating what it needs of the operator. Unknown is not false: a
## capability nobody established is reported as such, with the declaration that
## would settle it, rather than as a refusal to believe the operator.
require_capability <- function(A, name, method) {
  v <- capv(A, name)
  if (isTRUE(v)) return(invisible(TRUE))
  stopf(paste0("method '%s' requires the operator to be %s; it declares %s = %s.\n",
               "  %s"),
        method, name, name, format(v),
        if (is.na(v)) {
          sprintf(paste0("Unknown is not false. Declare it with ",
                         "linop(..., properties = c(%s = TRUE)) if you can justify it."), name)
        } else {
          "The operator declares the opposite."
        })
}

## Implications that hold by definition, applied once at construction so a user
## who declares positive_definite does not also have to declare
## positive_semidefinite.
close_caps <- function(caps) {
  pd <- cap_value(caps$positive_definite)
  if (isTRUE(pd) && is.na(cap_value(caps$positive_semidefinite))) {
    caps$positive_semidefinite <- capability(TRUE, ev_construction(
      depends_on = list(caps$positive_definite$evidence)))
  }
  dg <- cap_value(caps$diagonal)
  if (isTRUE(dg)) {
    for (nm in c("triangular_upper", "triangular_lower")) {
      if (is.na(cap_value(caps[[nm]]))) {
        caps[[nm]] <- capability(TRUE, ev_construction(
          depends_on = list(caps$diagonal$evidence)))
      }
    }
  }
  ## A real operator that is hermitian is symmetric, and conversely.
  if (isTRUE(cap_value(caps$real))) {
    if (isTRUE(cap_value(caps$hermitian)) && is.na(cap_value(caps$symmetric))) {
      caps$symmetric <- capability(TRUE, ev_construction(
        depends_on = list(caps$hermitian$evidence, caps$real$evidence)))
    }
    if (isTRUE(cap_value(caps$symmetric)) && is.na(cap_value(caps$hermitian))) {
      caps$hermitian <- capability(TRUE, ev_construction(
        depends_on = list(caps$symmetric$evidence, caps$real$evidence)))
    }
  }
  caps
}
