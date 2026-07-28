## Which knob is binding when eigs() stalls.
##
## On laplacian_1d(400) asking for 6 pairs at ncv = 40, the run stops after 240
## of at most 1000 iterations with nothing converged, and the certificate says
## "the subspace stopped improving". That reads like a detector giving up early.
## This decides it: hold everything else and raise the budget, then loosen the
## tolerance, then raise the subspace.
##
## Both ends of this spectrum are equally hard, and the gap line below says why:
## laplacian_1d is symmetric about the middle of its spectrum, so the relative
## gap at the top equals the one at the bottom.

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

n <- 400L
A <- laplacian_1d(n)
lam <- laplacian_1d_eigenvalues(n)
lo <- sort(lam)
hi <- sort(lam, decreasing = TRUE)

cat(sprintf("laplacian_1d(%d)\n", n))
cat(sprintf("  four smallest %s\n",
            paste(formatC(lo[1:4], format = "e", digits = 3), collapse = "  ")))
cat(sprintf("  four largest  %s\n",
            paste(formatC(hi[1:4], format = "e", digits = 6), collapse = "  ")))
cat(sprintf("  relative gap, top    (l_n - l_n-1) / (l_n - l_1) = %.3e\n",
            (hi[1] - hi[2]) / (max(lam) - min(lam))))
cat(sprintf("  relative gap, bottom (l_2 - l_1)   / (l_n - l_1) = %.3e\n\n",
            (lo[2] - lo[1]) / (max(lam) - min(lam))))

probe <- function(which, ncv, maxit, tol) {
  f <- eigs(A, k = 6L, which = which, ncv = ncv, maxit = maxit, tol = tol)
  truth <- if (which == "largest") hi[1:6] else lo[1:6]
  cat(sprintf("%-9s ncv %3d  maxit %5d  tol %.0e -> %4d iterations, %2d restarts, %d/%d converged, values to %.2e\n",
              which, ncv, maxit, tol, f$iterations, f$restarts, f$nconv, f$k,
              max(abs(f$values - truth) / abs(truth))))
}

## Five times the budget, then a looser tolerance, then a larger subspace.
probe("largest", 20L, 1000L, 1e-10)
probe("largest", 40L, 1000L, 1e-10)
probe("largest", 40L, 5000L, 1e-10)
probe("largest", 40L, 5000L, 1e-08)
probe("largest", 80L, 5000L, 1e-10)
cat("\n")
probe("smallest", 40L, 1000L, 1e-10)
probe("smallest", 40L, 5000L, 1e-10)
probe("smallest", 80L, 5000L, 1e-10)
