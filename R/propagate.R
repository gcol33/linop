## Section 5.6, corrected by S0.1 section 6. Every propagated capability records
## source = "construction" with depends_on holding the input capabilities the
## rule relied on, empty when the rule is unconditional.

## Structural identity, used to recognise A^H A and A A^H. Conservative: it may
## say FALSE for two operators that happen to be equal, never TRUE for two that
## differ.
same_operator <- function(A, B) {
  if (!is_linop(A) || !is_linop(B)) return(FALSE)
  if (!identical(A$node, B$node)) return(FALSE)
  if (!identical(A$dim, B$dim) || !identical(A$dtype, B$dtype)) return(FALSE)
  switch(A$node,
    dense    = identical(A$args$M, B$args$M),
    sparse   = identical(A$args$M, B$args$M),
    diag     = identical(A$args$d, B$args$d),
    identity = TRUE,
    fun      = identical(A$args$apply, B$args$apply) &&
               identical(A$args$adjoint, B$args$adjoint),
    transpose = , adjoint = , conjugate = same_operator(A$args$A, B$args$A),
    scale    = identical(A$args$a, B$args$a) && same_operator(A$args$A, B$args$A),
    sum      = ,
    product  = length(A$args$ops) == length(B$args$ops) &&
               all(vapply(seq_along(A$args$ops),
                          function(i) same_operator(A$args$ops[[i]], B$args$ops[[i]]),
                          logical(1))),
    ## Externally registered node types: fall back to comparing the stored
    ## arguments, so a third-party leaf still gets A^H A recognised.
    identical(A$args, B$args))
}

## Is B the adjoint (conjugate transpose) of A, provably and structurally?
is_adjoint_of <- function(B, A) {
  if (identical(B$node, "adjoint") && same_operator(B$args$A, A)) return(TRUE)
  ## For a real operator, transpose and adjoint have the same action.
  if (identical(B$node, "transpose") && isTRUE(capv(A, "real")) &&
      same_operator(B$args$A, A)) return(TRUE)
  FALSE
}

is_transpose_of <- function(B, A) {
  if (identical(B$node, "transpose") && same_operator(B$args$A, A)) return(TRUE)
  if (identical(B$node, "adjoint") && isTRUE(capv(A, "real")) &&
      same_operator(B$args$A, A)) return(TRUE)
  FALSE
}

## ------------------------------------------------------------------- sum ----

caps_sum <- function(ops) {
  caps <- list()
  for (nm in c("hermitian", "symmetric", "positive_definite", "positive_semidefinite", "real")) {
    vals <- lapply(ops, function(A) cap(A, nm))
    v <- tv_all(lapply(vals, cap_value))
    if (isTRUE(v)) {
      caps[[nm]] <- capability(TRUE, ev_construction(
        depends_on = Filter(Negate(is.null), lapply(vals, function(c) c$evidence))))
    } else {
      caps[[nm]] <- cap_unknown()
    }
  }
  ## Diagonal is closed under addition; triangularity likewise.
  for (nm in c("diagonal", "triangular_upper", "triangular_lower")) {
    vals <- lapply(ops, function(A) cap(A, nm))
    if (isTRUE(tv_all(lapply(vals, cap_value)))) {
      caps[[nm]] <- capability(TRUE, ev_construction(
        depends_on = Filter(Negate(is.null), lapply(vals, function(c) c$evidence))))
    }
  }
  caps
}

## ---------------------------------------------------------------- scale -----

caps_scale <- function(A, a) {
  caps <- list()
  a_real <- Im(a) == 0
  a_pos <- a_real && Re(a) > 0

  h <- cap(A, "hermitian")
  if (isTRUE(cap_value(h)) && a_real) {
    caps$hermitian <- capability(TRUE, ev_construction(depends_on = list(h$evidence)))
  }
  s <- cap(A, "symmetric")
  if (isTRUE(cap_value(s))) {
    caps$symmetric <- capability(TRUE, ev_construction(depends_on = list(s$evidence)))
  }
  for (nm in c("positive_definite", "positive_semidefinite")) {
    p <- cap(A, nm)
    if (isTRUE(cap_value(p)) && a_pos) {
      caps[[nm]] <- capability(TRUE, ev_construction(depends_on = list(p$evidence)))
    }
  }
  for (nm in c("diagonal", "triangular_upper", "triangular_lower", "full_rank")) {
    p <- cap(A, nm)
    if (isTRUE(cap_value(p)) && (nm != "full_rank" || a != 0)) {
      caps[[nm]] <- capability(TRUE, ev_construction(depends_on = list(p$evidence)))
    }
  }
  r <- cap(A, "real")
  if (isTRUE(cap_value(r)) && a_real) {
    caps$real <- capability(TRUE, ev_construction(depends_on = list(r$evidence)))
  }
  caps
}

## ------------------------------------------------------- transpose family ---

## adjoint(A): hermitian preserved, symmetric only if real, PD preserved.
caps_adjoint <- function(A) {
  caps <- list()
  keep <- function(nm) {
    p <- cap(A, nm)
    if (isTRUE(cap_value(p))) capability(TRUE, ev_construction(depends_on = list(p$evidence)))
    else NULL
  }
  for (nm in c("hermitian", "positive_definite", "positive_semidefinite",
               "unitary", "diagonal", "full_rank", "real")) {
    v <- keep(nm); if (!is.null(v)) caps[[nm]] <- v
  }
  s <- cap(A, "symmetric"); r <- cap(A, "real")
  if (isTRUE(cap_value(s)) && isTRUE(cap_value(r))) {
    caps$symmetric <- capability(TRUE, ev_construction(
      depends_on = list(s$evidence, r$evidence)))
  }
  ## Triangularity flips.
  tu <- cap(A, "triangular_upper"); tl <- cap(A, "triangular_lower")
  if (isTRUE(cap_value(tu))) caps$triangular_lower <- capability(TRUE, ev_construction(list(tu$evidence)))
  if (isTRUE(cap_value(tl))) caps$triangular_upper <- capability(TRUE, ev_construction(list(tl$evidence)))
  caps
}

## t(A): symmetric preserved, hermitian only if real.
caps_transpose <- function(A) {
  caps <- list()
  for (nm in c("symmetric", "unitary", "diagonal", "full_rank", "real")) {
    p <- cap(A, nm)
    if (isTRUE(cap_value(p))) caps[[nm]] <- capability(TRUE, ev_construction(list(p$evidence)))
  }
  h <- cap(A, "hermitian"); r <- cap(A, "real")
  if (isTRUE(cap_value(h)) && isTRUE(cap_value(r))) {
    caps$hermitian <- capability(TRUE, ev_construction(list(h$evidence, r$evidence)))
    for (nm in c("positive_definite", "positive_semidefinite")) {
      p <- cap(A, nm)
      if (isTRUE(cap_value(p))) caps[[nm]] <- capability(TRUE, ev_construction(list(p$evidence, r$evidence)))
    }
  }
  tu <- cap(A, "triangular_upper"); tl <- cap(A, "triangular_lower")
  if (isTRUE(cap_value(tu))) caps$triangular_lower <- capability(TRUE, ev_construction(list(tu$evidence)))
  if (isTRUE(cap_value(tl))) caps$triangular_upper <- capability(TRUE, ev_construction(list(tl$evidence)))
  caps
}

## Conj(A): everything preserved.
caps_conjugate <- function(A) {
  caps <- list()
  for (nm in c("hermitian", "symmetric", "positive_definite", "positive_semidefinite",
               "unitary", "diagonal", "triangular_upper", "triangular_lower",
               "full_rank", "real")) {
    p <- cap(A, nm)
    if (isTRUE(cap_value(p))) caps[[nm]] <- capability(TRUE, ev_construction(list(p$evidence)))
  }
  caps
}

## --------------------------------------------------------------- product ----

## A %*% B is NA in general. The two exceptions are structural and hold
## unconditionally: A^H A is hermitian and PSD whatever is known about A, and
## likewise A A^H.
caps_product <- function(ops) {
  caps <- list()
  if (length(ops) == 2L) {
    L <- ops[[1L]]; R <- ops[[2L]]
    gram_h <- is_adjoint_of(L, R)     # A^H A
    gram_t <- is_adjoint_of(R, L)     # A A^H
    if (gram_h || gram_t) {
      caps$hermitian <- capability(TRUE, ev_construction())
      caps$positive_semidefinite <- capability(TRUE, ev_construction())
      inner <- if (gram_h) R else L
      fr <- cap(inner, "full_rank")
      if (isTRUE(cap_value(fr))) {
        caps$positive_definite <- capability(TRUE, ev_construction(list(fr$evidence)))
      }
      if (isTRUE(capv(inner, "real"))) {
        caps$symmetric <- capability(TRUE, ev_construction(
          list(cap(inner, "real")$evidence)))
      }
    }
    ## t(A) %*% A is symmetric always; hermitian only when A is real, in which
    ## case the branch above has already fired.
    if (is_transpose_of(L, R) || is_transpose_of(R, L)) {
      caps$symmetric <- capability(TRUE, ev_construction())
    }
  }
  ## Structural closures that survive any product.
  for (nm in c("diagonal", "triangular_upper", "triangular_lower", "unitary", "real")) {
    vals <- lapply(ops, function(A) cap(A, nm))
    if (isTRUE(tv_all(lapply(vals, cap_value)))) {
      caps[[nm]] <- capability(TRUE, ev_construction(
        depends_on = Filter(Negate(is.null), lapply(vals, function(c) c$evidence))))
    }
  }
  caps
}
