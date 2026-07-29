`%||%` <- function(x, y) if (is.null(x)) y else x

## Section 5.2: X is always a matrix; a vector is a block of width 1.
## `dim<-` rather than as.matrix(): S0.2 measured 300 ns against 2300 ns, and this
## sits on the hot path of every apply.
as_block <- function(X) {
  if (is.null(dim(X))) dim(X) <- c(length(X), 1L)
  X
}

## Solvers hand back what a matrix would have returned: a single-column result
## from a vector input is a vector again.
undo_block <- function(Y, was_vector) {
  if (was_vector && is.matrix(Y) && ncol(Y) == 1L) dim(Y) <- NULL
  Y
}

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

## Column-wise reductions over a block, written in real arithmetic so they work
## on complex and real storage through one path. Re(<x, y>) = Re(x).Re(y) +
## Im(x).Im(y), which is the only part of the inner product the Krylov scalars
## use, since the quantities they form are real for a hermitian operator.
## x * x rather than x^2, which is the same bits at about twice the speed,
## measured 1.76x to 1.98x over n = 1e4 to 1e6 by
## dev_notes/spikes/allocation-sweep.R.
col_norms <- function(X) {
  if (is.complex(X)) {
    re <- Re(X); im <- Im(X)
    sqrt(colSums(re * re + im * im))
  } else {
    sqrt(colSums(X * X))
  }
}

col_dot <- function(X, Y) {
  if (is.complex(X) || is.complex(Y)) colSums(Re(X) * Re(Y) + Im(X) * Im(Y))
  else colSums(X * Y)
}

## Im(<x, y>), the part col_dot() discards. Zero for real storage, and nonzero
## for a quantity a hermitian operator requires to vanish, which is what makes it
## a contradiction test rather than a diagnostic.
col_dot_im <- function(X, Y) {
  if (is.complex(X) || is.complex(Y)) colSums(Re(X) * Im(Y) - Im(X) * Re(Y))
  else numeric(ncol(as_block(Y)))
}

## The whole inner product, <x, y> = conj(x)^T y, of which col_dot() is the real
## part and col_dot_im() the imaginary one. A hermitian operator makes the
## quantities its recurrences form real, which is why CG and MINRES need only the
## real part and are written in real arithmetic to get it. Arnoldi has no such
## guarantee: the Gram-Schmidt coefficient <v_i, A v_j> of a non-hermitian
## complex operator is complex, and discarding its imaginary part would
## orthogonalise against the wrong direction.
col_cdot <- function(X, Y) {
  if (is.complex(X) || is.complex(Y)) colSums(Conj(X) * Y) else colSums(X * Y)
}

## X^H Y and conj(X) Y as matrix products. Conj() on double storage returns the
## numbers it was given and allocates a full duplicate to do it, so the real
## branch here is not an approximation of the complex one: it is the same
## quantity, bitwise, without the copy.
##
## The copy is worth avoiding because these sit on apply and orthogonalisation
## paths where X is the largest array in the run. Orthogonalising against a
## stored basis duplicates the whole basis once per pass, 610 Mb at n = 1e6 and
## ncv = 80, measured at 18 to 26 per cent of eigs() runtime in self time by
## dev_notes/spikes/eigs-orth-crossover.R.
cross_adjoint <- function(X, Y) {
  if (is.complex(X)) crossprod(Conj(X), Y) else crossprod(X, Y)
}

conj_prod <- function(X, Y) {
  if (is.complex(X)) Conj(X) %*% Y else X %*% Y
}

## Y[, j] * v[j], without transposing twice.
##
## rep(v, each = nrow(X)) materialises a full n x k array to hold k distinct
## values, so a block one column wide pays an allocation the size of itself to
## express a scalar multiply. The guarded branch is the same arithmetic, bitwise,
## and a single right-hand side is the shape most solves have: this sits at 23 to
## 47 per cent of self time across the six solver loops, and the guard measures
## 3.6x to 8.9x on the helper (dev_notes/spikes/allocation-sweep.R).
scale_cols <- function(X, v) {
  if (length(v) == 1L) X * v else X * rep(v, each = nrow(X))
}

## A block of zeros in the storage mode of an existing one, so a recurrence that
## starts from zero does not silently demote a complex iterate to real.
zero_block <- function(X) {
  Z <- matrix(0, nrow(X), ncol(X))
  if (is.complex(X)) storage.mode(Z) <- "complex"
  Z
}

## A report is reproducible only if it fixes a seed, and a solver has no business
## moving the user's stream while doing it.
with_preserved_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  }
  set.seed(seed)
  expr
}

is_scalar_string <- function(x) is.character(x) && length(x) == 1L && !is.na(x)

## Three-valued logic. NA means unknown and never reads as FALSE (section 5.3).
tv_and <- function(a, b) {
  if (isFALSE(a) || isFALSE(b)) return(FALSE)
  if (is.na(a) || is.na(b)) return(NA)
  TRUE
}

tv_all <- function(xs) {
  if (!length(xs)) return(TRUE)
  Reduce(tv_and, xs)
}

## TRUE only when known TRUE; used where a rule needs a proven precondition.
known_true <- function(x) isTRUE(x)
