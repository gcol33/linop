## Several right-hand sides in lockstep against the same right-hand sides one
## after another.
##
## The claim CG established and the other six inherited is that k columns cost
## one block apply per step rather than k applies, because each column's
## recurrence is independent. What that does not buy is fewer column-applies: the
## block runs every column until the last one converges, so a column that was
## finished early keeps being carried. Both numbers are here, and they move in
## opposite directions.
##
## Two fixtures, because the two routes differ for two unrelated reasons. The
## matrix-free operator is an R closure, so the block saves interpreter calls;
## the stored one reaches BLAS, so the block turns k matrix-vector products into
## one matrix-matrix product. Neither fixture measures the other's effect.
##
## The answers are compared as well as the timings. A route that is faster and
## returns something else is not a route.

BLOCK_WIDTHS <- c(1L, 2L, 4L, 8L, 16L)
BLOCK_TOL <- 1e-8

block_cases <- function() {
  list(
    list(label = "laplacian_1d(200), matrix free",
         n = 200L,
         build = function(cn) {
           counted_linop(laplacian_1d_apply, laplacian_1d_apply, c(200L, 200L), cn,
                         properties = c(hermitian = TRUE, positive_definite = TRUE))
         }),
    list(label = "spd_prescribed(300, kappa 1e4), stored",
         n = 300L,
         build = function(cn) {
           M <- spd_prescribed(300L, 10^seq(-4, 0, length.out = 300L), seed = 5L)
           counted_dense(M, cn, properties = c(hermitian = TRUE,
                                               positive_definite = TRUE))
         }))
}

bench_blocks <- function() {
  out <- bench_rows()

  for (case in block_cases()) {
    cat(sprintf("  %s\n", case$label))
    cn <- new_counter()
    A <- case$build(cn)
    n <- case$n
    maxit <- 10L * n

    for (k in BLOCK_WIDTHS) {
      set.seed(11)
      B <- matrix(stats::rnorm(n * k), n, k)

      lock <- bench_run(solve(A, B, method = "cg", tol = BLOCK_TOL, maxit = maxit,
                              details = TRUE),
                        counter = cn)
      each <- bench_run(lapply(seq_len(k), function(j)
                          solve(A, B[, j, drop = FALSE], method = "cg",
                                tol = BLOCK_TOL, maxit = maxit, details = TRUE)),
                        counter = cn)

      X_lock <- lock$value$x
      X_each <- do.call(cbind, lapply(each$value, function(f) f$x))
      agree <- max(vapply(seq_len(k), function(j)
        rel_error(X_each[, j], X_lock[, j]), numeric(1)))

      out$add(fixture = case$label, n = n, k = k, route = "lockstep",
              iterations = lock$value$iterations, applies = lock$applies,
              column_applies = lock$cols, seconds = lock$seconds,
              agreement_vs_lockstep = NA_real_)
      out$add(fixture = case$label, n = n, k = k, route = "one column at a time",
              iterations = sum(vapply(each$value, function(f) f$iterations, integer(1))),
              applies = each$applies, column_applies = each$cols,
              seconds = each$seconds, agreement_vs_lockstep = agree)

      cat(sprintf("    k = %2d   lockstep %6d applies %7.3f s   separate %6d applies %7.3f s   agree to %.1e\n",
                  k, lock$applies, lock$seconds, each$applies, each$seconds, agree))
    }
  }

  out$collect()
}
