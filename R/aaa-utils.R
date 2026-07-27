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
col_norms <- function(X) {
  if (is.complex(X)) sqrt(colSums(Re(X)^2 + Im(X)^2)) else sqrt(colSums(X^2))
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

## Y[, j] * v[j], without transposing twice.
scale_cols <- function(X, v) X * rep(v, each = nrow(X))

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
