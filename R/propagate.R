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
    perm     = identical(A$args$p, B$args$p),
    lowrank  = identical(A$args$U, B$args$U) && identical(A$args$V, B$args$V) &&
               identical(A$args$d, B$args$d),
    power    = identical(A$args$k, B$args$k) && same_operator(A$args$A, B$args$A),
    inverse  = same_operator(A$args$A, B$args$A),
    sum      = , product = , kron = , blockdiag = , hstack = ,
    vstack   = length(A$args$ops) == length(B$args$ops) &&
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

## A capability every operand holds, with all their evidence underneath it, or
## NULL. The closure rule of a sum, a product, a Kronecker product or a block
## diagonal is this test and nothing more.
cap_all <- function(ops, nm) {
  vals <- lapply(ops, function(A) cap(A, nm))
  if (!isTRUE(tv_all(lapply(vals, cap_value)))) return(NULL)
  capability(TRUE, ev_construction(
    depends_on = Filter(Negate(is.null), lapply(vals, function(c) c$evidence))))
}

## ------------------------------------------------------------------- sum ----

caps_sum <- function(ops) {
  caps <- list()
  for (nm in c("hermitian", "symmetric", "positive_definite", "positive_semidefinite", "real")) {
    caps[[nm]] <- cap_all(ops, nm) %||% cap_unknown()
  }
  ## Diagonal is closed under addition; triangularity likewise.
  for (nm in c("diagonal", "triangular_upper", "triangular_lower")) {
    caps[[nm]] <- cap_all(ops, nm)
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
      ## X^H X is definite exactly when X has full column rank, X X^H when it has
      ## full row rank. full_rank is rank min(m, n), which is the right one only
      ## when the dimension the Gram matrix lives on is the smaller: a wide X of
      ## full rank has a singular X^H X.
      side <- if (gram_h) inner$dim[2L] else inner$dim[1L]
      if (isTRUE(cap_value(fr)) && side <= min(inner$dim)) {
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
    caps[[nm]] <- cap_all(ops, nm)
  }
  caps
}

## ------------------------------------------------------ Kronecker product ---

## The eigenvalues of A (x) B are the products lambda_i mu_j and its singular
## values the products sigma_i tau_j, (A (x) B)^H = A^H (x) B^H, and the block
## structure keeps diagonal and triangular shapes. So every capability below
## holds when it holds for both factors. Full rank needs both factors square or
## both tall or both wide: a 1 x 2 and a 2 x 1 of full rank give a 2 x 2 of
## rank 1.
caps_kron <- function(ops) {
  caps <- list()
  for (nm in c("hermitian", "symmetric", "positive_definite", "positive_semidefinite",
               "unitary", "diagonal", "triangular_upper", "triangular_lower", "real")) {
    caps[[nm]] <- cap_all(ops, nm)
  }
  if (same_orientation(ops)) caps$full_rank <- cap_all(ops, "full_rank")
  caps
}

same_orientation <- function(ops) {
  tall <- vapply(ops, function(A) A$dim[1L] >= A$dim[2L], logical(1))
  wide <- vapply(ops, function(A) A$dim[1L] <= A$dim[2L], logical(1))
  all(tall) || all(wide)
}

## ---------------------------------------------------------- block diagonal ---

## A block diagonal operator acts on each block separately, so its spectrum,
## singular values, adjoint and shape are those of the blocks. Full rank under
## the same orientation condition as the Kronecker product.
caps_blockdiag <- function(ops) caps_kron(ops)

## ------------------------------------------------------------------ stacks ---

## [A_1 ... A_k] has the sum of the ranges as its range, so one block of full
## row rank (full rank and no taller than wide) makes a wide stack of full rank.
## Stacking rows is the transpose of that. Nothing about symmetry follows.
caps_stack <- function(ops, direction) {
  caps <- list(real = cap_all(ops, "real"))
  rows <- vapply(ops, function(A) A$dim[1L], numeric(1))
  cols <- vapply(ops, function(A) A$dim[2L], numeric(1))
  m <- if (direction == "h") rows[1L] else sum(rows)
  n <- if (direction == "h") sum(cols) else cols[1L]
  if ((direction == "h" && m <= n) || (direction == "v" && m >= n)) {
    for (A in ops) {
      fits <- if (direction == "h") A$dim[1L] <= A$dim[2L] else A$dim[1L] >= A$dim[2L]
      fr <- cap(A, "full_rank")
      if (fits && isTRUE(cap_value(fr))) {
        caps$full_rank <- capability(TRUE, ev_construction(list(fr$evidence)))
        break
      }
    }
  }
  caps
}

## ------------------------------------------------------------------- power ---

## A^k for k >= 1: hermitian, symmetric, unitary, diagonal, triangular, real and
## full rank survive any power, and a positive definite operator stays so. A
## hermitian operator to an even power is positive semidefinite, and definite
## when it is also of full rank, since its eigenvalues are lambda^k.
caps_power <- function(A, k) {
  caps <- list()
  for (nm in c("hermitian", "symmetric", "positive_definite", "positive_semidefinite",
               "unitary", "diagonal", "triangular_upper", "triangular_lower",
               "full_rank", "real")) {
    caps[[nm]] <- cap_all(list(A), nm)
  }
  h <- cap(A, "hermitian")
  if (k %% 2L == 0L && isTRUE(cap_value(h))) {
    caps$positive_semidefinite <- capability(TRUE, ev_construction(list(h$evidence)))
    fr <- cap(A, "full_rank")
    if (isTRUE(cap_value(fr))) {
      caps$positive_definite <- capability(TRUE, ev_construction(
        list(h$evidence, fr$evidence)))
    }
  }
  caps
}

## ----------------------------------------------------------------- inverse ---

## A^-1 inherits every shape and symmetry A has: the inverse of a triangular,
## diagonal, unitary, hermitian or positive definite matrix is one. It is of
## full rank by existing, and that evidence is the factorisation that showed it.
caps_inverse <- function(A, invertible_evidence) {
  caps <- list()
  for (nm in c("hermitian", "symmetric", "positive_definite", "positive_semidefinite",
               "unitary", "diagonal", "triangular_upper", "triangular_lower", "real")) {
    caps[[nm]] <- cap_all(list(A), nm)
  }
  caps$full_rank <- capability(TRUE, invertible_evidence)
  caps
}
