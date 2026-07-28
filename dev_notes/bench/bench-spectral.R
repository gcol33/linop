## What the two spectral verbs cost, against operators whose spectra are known in
## closed form, with a dense factorisation of the same operator on the row below.
##
## ncv is varied rather than fixed because it is the knob the reference label is
## about: a round holds its whole basis and orthogonalises against all of it, so
## a larger subspace costs O(ncv) vectors of storage and O(ncv^2 n) of
## arithmetic. Every setting is run at the same budget, so what the ladder
## compares is the subspace and not maxit.
##
## On laplacian_1d(400) the subspace is the only knob that matters, and the two
## ends of that spectrum are equally hard: the relative gap is 4.6e-5 at both. At
## ncv = 40 the run stalls at 240 iterations with nothing converged, and raising
## maxit from 1000 to 5000, or loosening tol from 1e-10 to 1e-8, changes neither
## the iteration count nor the value error. At ncv = 80 the same request reaches
## 6 of 6 at 1.1e-16 in 477 iterations. The stall detector stopping a quarter of
## the way into the budget is therefore reporting a fact rather than giving up
## early: there was nothing left in that subspace to find.

##
## The shift-invert row prices the whole transformation, inner solves included:
## the counter sits on the operator, and every MINRES step inside
## (A - sigma I)^-1 goes through it. That is the number to read when deciding
## whether a shift is worth taking matrix-free, and it is not small.

SPECTRAL_NCV <- c(20L, 40L, 80L)
SPECTRAL_TOL <- 1e-10

spectral_dense <- function(n) {
  M <- diag(2, n)
  idx <- seq_len(n - 1L)
  M[cbind(idx, idx + 1L)] <- -1
  M[cbind(idx + 1L, idx)] <- -1
  M
}

bench_spectral <- function() {
  out <- bench_rows()

  ## fit is the spectral result where one exists and NULL for the dense routes,
  ## whose value is a bare vector of eigen- or singular values.
  add_row <- function(problem, verb, k, which, ncv, route, run, values, truth,
                      fit = NULL, note = "") {
    err <- max(abs(values - truth) / abs(truth))
    out$add(problem = problem, verb = verb, k = k, which = which,
            ncv = if (is.na(ncv)) NA_integer_ else as.integer(ncv),
            route = route,
            iterations = if (is.null(fit)) NA_integer_ else fit$iterations,
            restarts = if (is.null(fit)) NA_integer_ else fit$restarts,
            nconv = if (is.null(fit)) NA_integer_ else fit$nconv,
            applies = run$applies, seconds = run$seconds,
            backward_error = if (is.null(fit)) NA_real_
                             else max(fit$certificate$values$backward_error),
            value_error = err, note = note)
    cat(sprintf("    %-34s %6s applies %8.3f s  %s converged  values to %.2e\n",
                route, format(run$applies), run$seconds,
                if (is.null(fit)) "  -" else sprintf("%d/%d", fit$nconv, fit$k),
                err))
  }

  ## ----------------------------------------- eigs, both ends, three subspaces
  n <- 400L
  lam <- laplacian_1d_eigenvalues(n)
  M <- spectral_dense(n)
  k <- 6L
  maxit <- 1000L

  for (which in c("largest", "smallest")) {
    truth <- if (which == "largest") sort(lam, decreasing = TRUE)[seq_len(k)]
             else sort(lam)[seq_len(k)]
    cat(sprintf("  eigs(laplacian_1d(%d), k = %d, which = '%s')\n", n, k, which))
    for (ncv in SPECTRAL_NCV) {
      cnt <- new_counter()
      A <- counted_linop(laplacian_1d_apply, laplacian_1d_apply, c(n, n), cnt,
                         properties = c(hermitian = TRUE, positive_definite = TRUE))
      run <- bench_run(eigs(A, k = k, which = which, ncv = ncv, maxit = maxit,
                            tol = SPECTRAL_TOL),
                       counter = cnt)
      add_row(sprintf("laplacian_1d(%d)", n), "eigs", k, which, ncv,
              sprintf("lanczos, ncv = %d", ncv), run, run$value$values, truth,
              fit = run$value)
    }
    dense <- bench_run(eigen(M, symmetric = TRUE, only.values = TRUE)$values)
    dv <- if (which == "largest") sort(dense$value, decreasing = TRUE)[seq_len(k)]
          else sort(dense$value)[seq_len(k)]
    add_row(sprintf("laplacian_1d(%d)", n), "eigen()", k, which, NA_integer_,
            "dense eigen(), whole spectrum", dense, dv, truth,
            note = "returns all n values")
  }

  ## ------------------------------------------------------------ shift-invert
  ## Smaller, and with the shift away from the clustered end: sigma at the bottom
  ## of this spectrum puts kappa(A - sigma I) past 1e5 and the inner solve stops
  ## being the cheap part of anything.
  ns <- 100L
  sigma <- 1.0
  lam_s <- laplacian_1d_eigenvalues(ns)
  ks <- 4L
  truth_s <- lam_s[order(abs(lam_s - sigma))][seq_len(ks)]
  cat(sprintf("  eigs(laplacian_1d(%d), k = %d, sigma = %g)\n", ns, ks, sigma))
  cnt <- new_counter()
  As <- counted_linop(laplacian_1d_apply, laplacian_1d_apply, c(ns, ns), cnt,
                      properties = c(hermitian = TRUE, positive_definite = TRUE))
  run <- bench_run(eigs(As, k = ks, sigma = sigma, ncv = 20L, maxit = 200L,
                        tol = SPECTRAL_TOL),
                   counter = cnt)
  add_row(sprintf("laplacian_1d(%d)", ns), "eigs", ks, sprintf("sigma = %g", sigma),
          20L, "lanczos on (A - sigma I)^-1", run,
          run$value$values[order(abs(run$value$values - sigma))][seq_len(ks)],
          truth_s, fit = run$value,
          note = "applies include every inner minres step")

  ## ------------------------------------------- svds, matrix free and rank deficient
  nd <- 400L
  sv <- sort(diff_1d_singular_values(nd), decreasing = TRUE)[seq_len(k)]
  D <- diag(1, nd)
  D <- D[-1L, , drop = FALSE] - D[-nd, , drop = FALSE]
  cat(sprintf("  svds(diff_1d(%d), k = %d, which = 'largest')\n", nd, k))
  for (ncv in SPECTRAL_NCV) {
    cnt <- new_counter()
    Ad <- counted_linop(function(X) X[-1L, , drop = FALSE] - X[-nd, , drop = FALSE],
                        function(X) rbind(0, X) - rbind(X, 0),
                        c(nd - 1L, nd), cnt)
    run <- bench_run(svds(Ad, k = k, which = "largest", ncv = ncv, maxit = maxit,
                          tol = SPECTRAL_TOL),
                     counter = cnt)
    add_row(sprintf("diff_1d(%d)", nd), "svds", k, "largest", ncv,
            sprintf("golub-kahan, ncv = %d", ncv), run, run$value$d, sv,
            fit = run$value)
  }
  dense <- bench_run(svd(D, nu = 0, nv = 0)$d)
  add_row(sprintf("diff_1d(%d)", nd), "svd()", k, "largest", NA_integer_,
          "dense svd(), whole spectrum", dense, dense$value[seq_len(k)], sv,
          note = "returns all singular values")

  ## Stored, rectangular, with the singular values dialled in.
  m <- 600L
  nr <- 200L
  f <- lsq_prescribed(m, nr, 10^seq(0, -3, length.out = nr), seed = 7L)
  cat(sprintf("  svds(lsq_prescribed(%d x %d), k = %d)\n", m, nr, k))
  for (ncv in SPECTRAL_NCV) {
    cnt <- new_counter()
    Al <- counted_dense(f$A, cnt)
    run <- bench_run(svds(Al, k = k, which = "largest", ncv = ncv, maxit = maxit,
                          tol = SPECTRAL_TOL),
                     counter = cnt)
    add_row(sprintf("lsq_prescribed(%d x %d)", m, nr), "svds", k, "largest", ncv,
            sprintf("golub-kahan, ncv = %d", ncv), run, run$value$d,
            f$sigma[seq_len(k)], fit = run$value)
  }
  dense <- bench_run(svd(f$A, nu = 0, nv = 0)$d)
  add_row(sprintf("lsq_prescribed(%d x %d)", m, nr), "svd()", k, "largest",
          NA_integer_, "dense svd(), whole spectrum", dense,
          dense$value[seq_len(k)], f$sigma[seq_len(k)],
          note = "returns all singular values")

  out$collect()
}
