## The same system solved three ways as n grows: through the operator's action,
## through a dense matrix, and through a sparse one.
##
## The operator is tridiag(-1, 3, -1), the 1-D Laplacian shifted by one, and the
## shift is the whole reason it is not the plain Laplacian. The Laplacian's
## condition number grows like n^2, so CG's iteration count grows with the ladder
## and a table built on it would be measuring conditioning while claiming to
## measure size. Shifted, the spectrum sits in (1, 5) for every n, the iteration
## count is flat down the ladder -- which the results record, so the claim is
## checked rather than asserted -- and what is left moving is size.
##
## The sparse column is the one to read carefully. A tridiagonal matrix has an
## O(n) direct solve, and it wins here at every size. That is the honest shape of
## the comparison: a matrix-free method is not a faster route to a banded system,
## it is the only route to an operator that has no matrix. The dense column shows
## where that starts to matter.
##
## object_mb is object.size() of whatever the route holds: the closure and its
## counter for the operator, the matrix for the other two. It is the quantity the
## ladder is about, and it is exact.

MF_LADDER <- c(500L, 1000L, 2000L, 4000L, 8000L, 16000L, 50000L, 100000L)

## Above this the dense route is not attempted. At n = 16000 the matrix alone is
## 1.9 GB before the factorisation copies it; the row records that, and the size
## it would have been, rather than quietly ending the column.
MF_DENSE_MAX <- 8000L
MF_TOL <- 1e-8

mf_apply <- function(X) laplacian_1d_apply(X) + X

mf_dense <- function(n) {
  M <- diag(3, n)
  idx <- seq_len(n - 1L)
  M[cbind(idx, idx + 1L)] <- -1
  M[cbind(idx + 1L, idx)] <- -1
  M
}

mf_sparse <- function(n) {
  Matrix::bandSparse(n, n, k = c(-1, 0, 1),
                     diagonals = list(rep(-1, n - 1L), rep(3, n),
                                      rep(-1, n - 1L)),
                     symmetric = FALSE)
}

bench_matrixfree <- function() {
  out <- bench_rows()
  have_matrix <- requireNamespace("Matrix", quietly = TRUE)
  if (!have_matrix) cat("  Matrix is not installed; the sparse column is absent\n")

  ## The sparse direct solve is what the forward errors below are measured
  ## against, so it is checked against the closed form once, at the smallest n,
  ## before anything rests on it.
  if (have_matrix) {
    n0 <- MF_LADDER[1L]
    set.seed(21)
    b0 <- stats::rnorm(n0)
    ref <- as.vector(Matrix::solve(mf_sparse(n0), b0))
    cat(sprintf("  sparse LU against the closed form at n = %d: %.2e\n",
                n0, rel_error(ref, shifted_laplacian_solve(n0, -1, b0))))
  }

  for (n in MF_LADDER) {
    set.seed(21)
    b <- stats::rnorm(n)
    cn <- new_counter()
    A <- counted_linop(mf_apply, mf_apply, c(n, n), cn,
                       properties = c(hermitian = TRUE, positive_definite = TRUE))

    truth <- if (have_matrix) as.vector(Matrix::solve(mf_sparse(n), b)) else NULL

    free <- bench_run(solve(A, b, method = "cg", tol = MF_TOL, maxit = 10L * n,
                            details = TRUE),
                      counter = cn)
    out$add(n = n, route = "matrix free, cg", reps = BENCH_REPS,
            seconds = free$seconds, iterations = free$value$iterations,
            applies = free$applies,
            object_mb = as.numeric(utils::object.size(A)) / 2^20,
            backward_error = max(free$value$certificate$values$backward_error),
            forward_error = if (is.null(truth)) NA_real_
                            else rel_error(free$value$x, truth),
            note = "")
    cat(sprintf("    n = %6d   matrix free %7.3f s  %3d iterations  %5d applies\n",
                n, free$seconds, free$value$iterations, free$applies))

    if (have_matrix) {
      sp <- bench_run(Matrix::solve(mf_sparse(n), b))
      out$add(n = n, route = "sparse, form and solve", reps = BENCH_REPS,
              seconds = sp$seconds, iterations = NA_integer_,
              applies = NA_integer_,
              object_mb = as.numeric(utils::object.size(mf_sparse(n))) / 2^20,
              backward_error = NA_real_,
              forward_error = 0, note = "the reference the forward errors use")
      cat(sprintf("    n = %6d   sparse      %7.3f s\n", n, sp$seconds))
    }

    if (n <= MF_DENSE_MAX) {
      ## One timed run once the factorisation is the dominant cost; the count is
      ## recorded in the reps column rather than left to be assumed.
      dreps <- if (n <= 2000L) BENCH_REPS else 1L
      form <- bench_run(mf_dense(n), reps = dreps)
      M <- form$value
      dense <- bench_run(solve(M, b), reps = dreps)
      out$add(n = n, route = "dense, form the matrix", reps = dreps,
              seconds = form$seconds, iterations = NA_integer_,
              applies = NA_integer_,
              object_mb = as.numeric(utils::object.size(M)) / 2^20,
              backward_error = NA_real_,
              forward_error = NA_real_, note = "")
      out$add(n = n, route = "dense, LU solve", reps = dreps,
              seconds = dense$seconds, iterations = NA_integer_,
              applies = NA_integer_, object_mb = NA_real_,
              backward_error = NA_real_,
              forward_error = if (is.null(truth)) NA_real_
                              else rel_error(dense$value, truth),
              note = "the factorisation copies the matrix")
      cat(sprintf("    n = %6d   dense       %7.3f s to form, %7.3f s to solve\n",
                  n, form$seconds, dense$seconds))
      rm(M)
    } else {
      out$add(n = n, route = "dense, form the matrix", reps = NA_integer_,
              seconds = NA_real_, iterations = NA_integer_, applies = NA_integer_,
              object_mb = 8 * (as.numeric(n)^2) / 2^20,
              backward_error = NA_real_, forward_error = NA_real_,
              note = sprintf("not run: the matrix is %.1f GB before the factorisation copies it",
                             8 * (as.numeric(n)^2) / 2^30))
    }
    invisible(gc(verbose = FALSE))
  }

  out$collect()
}
