## Section 5.5. No path silently downcasts.

DTYPES <- c("double", "complex")

check_dtype <- function(dtype) {
  if (!is_scalar_string(dtype) || !dtype %in% DTYPES) {
    stopf("dtype must be one of: %s", paste(DTYPES, collapse = ", "))
  }
  dtype
}

#' Promote two scalar types
#'
#' @param a,b Scalar type labels.
#' @return `"complex"` if either is complex, else `"double"`.
#' @keywords internal
promote <- function(a, b) {
  check_dtype(a); check_dtype(b)
  if (a == "complex" || b == "complex") "complex" else "double"
}

promote_all <- function(xs) Reduce(promote, xs)

dtype_of_data <- function(x) if (is.complex(x)) "complex" else "double"

## A real operator applied to a complex block promotes the result rather than
## dropping the imaginary part.
cast_block <- function(X, dtype) {
  if (dtype == "complex" && !is.complex(X)) {
    storage.mode(X) <- "complex"
  } else if (dtype == "double" && is.complex(X)) {
    stopf("refusing to downcast a complex block to double")
  }
  X
}

## The dtype an apply must produce given the operator dtype and the input block.
result_dtype <- function(op_dtype, X) promote(op_dtype, dtype_of_data(X))
