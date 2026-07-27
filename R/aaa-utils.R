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
