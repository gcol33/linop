## Section 5.3. Evidence is three independent fields, not a strength ranking,
## because the classes are not comparable: a probabilistic bound can carry a
## rigorous confidence statement while a deterministic estimate carries no
## guarantee, and an adapter contract can be stronger than a hundred probes.

EVIDENCE_SOURCES <- c("construction", "adapter_contract", "user_declaration",
                      "probe", "computation", "theorem")

EVIDENCE_GUARANTEES <- c("identity", "deterministic_bound", "probabilistic_bound",
                         "estimate", "heuristic")

#' Describe the argument behind a capability value
#'
#' @param source One of `construction`, `adapter_contract`, `user_declaration`,
#'   `probe`, `computation`, `theorem`.
#' @param guarantee One of `identity`, `deterministic_bound`,
#'   `probabilistic_bound`, `estimate`, `heuristic`.
#' @param confidence `1`, a stated probability, or `NA`.
#' @param depends_on List of evidence objects this one rests on. Empty when the
#'   rule is unconditional.
#' @return An object of class `linop_evidence`.
#' @export
evidence <- function(source, guarantee, confidence = NA_real_, depends_on = list()) {
  if (!is_scalar_string(source) || !source %in% EVIDENCE_SOURCES) {
    stopf("source must be one of: %s", paste(EVIDENCE_SOURCES, collapse = ", "))
  }
  if (!is_scalar_string(guarantee) || !guarantee %in% EVIDENCE_GUARANTEES) {
    stopf("guarantee must be one of: %s", paste(EVIDENCE_GUARANTEES, collapse = ", "))
  }
  if (!is.list(depends_on)) stopf("depends_on must be a list of evidence objects")
  for (d in depends_on) {
    if (!inherits(d, "linop_evidence")) stopf("every depends_on element must be evidence()")
  }
  structure(
    list(source = source, guarantee = guarantee,
         confidence = as.numeric(confidence)[1L], depends_on = depends_on),
    class = "linop_evidence")
}

#' A capability value together with its evidence
#'
#' @param value `TRUE`, `FALSE`, or `NA`. `NA` means unknown and is never read
#'   as `FALSE`.
#' @param evidence An object from [evidence()], or `NULL` when `value` is `NA`.
#' @return An object of class `linop_capability`.
#' @export
capability <- function(value, evidence = NULL) {
  if (!is.logical(value) || length(value) != 1L) {
    stopf("capability value must be a single logical (TRUE, FALSE or NA)")
  }
  if (is.na(value)) {
    evidence <- NULL
  } else if (is.null(evidence)) {
    stopf("a capability with a known value needs evidence")
  } else if (!inherits(evidence, "linop_evidence")) {
    stopf("evidence must come from evidence()")
  }
  structure(list(value = value, evidence = evidence), class = "linop_capability")
}

cap_unknown <- function() capability(NA)

cap_value <- function(cap) if (is.null(cap)) NA else cap$value

## Unconditional construction: the rule holds whatever the inputs claim.
ev_construction <- function(depends_on = list()) {
  evidence("construction", "identity", 1, depends_on)
}

ev_declared <- function() evidence("user_declaration", "identity", 1)
ev_adapter <- function() evidence("adapter_contract", "identity", 1)
ev_probe <- function() evidence("probe", "heuristic", NA_real_)

#' State what evidence a solver will accept
#'
#' @param sources Acceptable evidence sources.
#' @param guarantees Acceptable guarantee kinds.
#' @param min_confidence Minimum confidence, or `NA` to not require one.
#' @return An object of class `linop_requirement`.
#' @export
requirement <- function(sources = EVIDENCE_SOURCES,
                        guarantees = EVIDENCE_GUARANTEES,
                        min_confidence = NA_real_) {
  bad <- setdiff(sources, EVIDENCE_SOURCES)
  if (length(bad)) stopf("unknown evidence source: %s", paste(bad, collapse = ", "))
  bad <- setdiff(guarantees, EVIDENCE_GUARANTEES)
  if (length(bad)) stopf("unknown guarantee: %s", paste(bad, collapse = ", "))
  structure(list(sources = sources, guarantees = guarantees,
                 min_confidence = as.numeric(min_confidence)[1L]),
            class = "linop_requirement")
}

#' Does this evidence satisfy a requirement?
#'
#' Recurses into `depends_on`, so a composite built by an unconditional rule from
#' user-declared inputs fails any requirement the inputs would have failed
#' directly. Without the recursion, `crossprod(A) + crossprod(B)` would report
#' bare `construction` and launder the declarations underneath it.
#'
#' @param actual Evidence from [evidence()], or `NULL`.
#' @param required A [requirement()].
#' @return `TRUE` or `FALSE`.
#' @export
evidence_satisfies <- function(actual, required) {
  if (!inherits(required, "linop_requirement")) stopf("required must come from requirement()")
  if (is.null(actual)) return(FALSE)
  if (!inherits(actual, "linop_evidence")) stopf("actual must be evidence() or NULL")

  if (!actual$source %in% required$sources) return(FALSE)
  if (!actual$guarantee %in% required$guarantees) return(FALSE)
  if (!is.na(required$min_confidence)) {
    if (is.na(actual$confidence) || actual$confidence < required$min_confidence) return(FALSE)
  }
  for (d in actual$depends_on) {
    if (!evidence_satisfies(d, required)) return(FALSE)
  }
  TRUE
}

#' @export
print.linop_evidence <- function(x, ...) {
  cat(format_evidence(x), "\n", sep = "")
  invisible(x)
}

format_evidence <- function(x, depth = 0L) {
  if (is.null(x)) return("<none>")
  conf <- if (is.na(x$confidence)) "NA" else format(x$confidence)
  s <- sprintf("%s / %s / conf %s", x$source, x$guarantee, conf)
  if (length(x$depends_on)) {
    kids <- vapply(x$depends_on, format_evidence, character(1), depth = depth + 1L)
    s <- paste0(s, " <- [", paste(kids, collapse = "; "), "]")
  }
  s
}

#' @export
print.linop_capability <- function(x, ...) {
  cat(sprintf("%s  (%s)\n", x$value, format_evidence(x$evidence)))
  invisible(x)
}

## The weakest link across a set of evidence, used when a rule is conditional on
## several inputs: the composite depends on all of them.
combine_evidence <- function(...) {
  deps <- Filter(Negate(is.null), list(...))
  ev_construction(depends_on = deps)
}
