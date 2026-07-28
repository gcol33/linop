## The measurements behind dev_notes/lsqr-and-the-least-squares-certificate.md.
##
## Two properties of LSQR that neither CG, MINRES nor GMRES showed, both about
## how little of the recurrence survives to be reported:
##
##   1. two arithmetically equivalent implementations part company by O(1)
##      relative in mid-flight, and agree again once both converge;
##   2. finite termination at n steps does not survive floating point.
##
## Run from the package root with the Windows R.
devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")
lsqr <- linop:::lsqr_solve

## Paige and Saunders 1982, single column, no stopping test, no preconditioner.
reference_lsqr <- function(M, b, steps) {
  n <- ncol(M)
  Mh <- Conj(t(M))
  x <- rep(0, n)
  beta <- sqrt(sum(Mod(b)^2)); u <- b / beta
  v <- as.vector(Mh %*% u); alpha <- sqrt(sum(Mod(v)^2)); v <- v / alpha
  w <- v
  phibar <- beta; rhobar <- alpha
  for (j in seq_len(steps)) {
    ut <- as.vector(M %*% v) - alpha * u
    beta <- sqrt(sum(Mod(ut)^2))
    if (beta == 0) break
    u <- ut / beta
    vt <- as.vector(Mh %*% u) - beta * v
    alpha <- sqrt(sum(Mod(vt)^2))
    rho <- sqrt(rhobar^2 + beta^2)
    cs <- rhobar / rho; sn <- beta / rho
    theta <- sn * alpha; rhobar <- -cs * alpha
    phi <- cs * phibar; phibar <- sn * phibar
    x <- x + (phi / rho) * w
    if (alpha == 0) break
    v <- vt / alpha
    w <- v - (theta / rho) * w
  }
  matrix(x, n, 1L)
}

## 45 x 15, singular values from 0.5 to 20, so kappa = 40 and the operator is
## well conditioned. Nothing below is about ill conditioning.
fixture <- function(seed) {
  f <- lsq_prescribed(45, 15, exp(seq(log(0.5), log(20), length.out = 15)), seed = seed)
  set.seed(600L + seed)
  list(f = f, A = linop(f$A), b = matrix(stats::rnorm(45), 45, 1))
}

cat("== 1. divergence from the published recurrence, relative, by step ==\n")
steps <- c(2, 4, 6, 8, 10, 12, 14, 16)
tab <- matrix(NA_real_, length(steps), 5,
              dimnames = list(paste("step", steps), paste("seed", 1:5)))
for (si in 1:5) {
  fx <- fixture(si)
  for (ki in seq_along(steps)) {
    mine <- lsqr(fx$A, fx$b, tol = 0, maxit = steps[ki])$x
    ref <- reference_lsqr(fx$f$A, fx$b, steps[ki])
    tab[ki, si] <- max(Mod(mine - ref)) / max(Mod(ref))
  }
}
print(signif(tab, 3))

cat("\n== 2. the published recurrence against closed-form truth, by step ==\n")
steps2 <- c(14, 15, 16, 18, 20, 25, 30)
tab2 <- matrix(NA_real_, length(steps2), 5,
               dimnames = list(paste("step", steps2), paste("seed", 1:5)))
mine_it <- integer(5)
for (si in 1:5) {
  fx <- fixture(si)
  truth <- lsq_prescribed_solve(fx$f, fx$b)
  for (ki in seq_along(steps2)) {
    ref <- reference_lsqr(fx$f$A, fx$b, steps2[ki])
    tab2[ki, si] <- max(Mod(ref - truth)) / max(Mod(truth))
  }
  mine_it[si] <- lsqr(fx$A, fx$b, tol = 1e-12, maxit = 500L)$iterations
}
print(signif(tab2, 3))
cat(sprintf("this implementation reaches 1e-12 in %s iterations\n",
            paste(mine_it, collapse = ", ")))

cat("\n== 3. the arithmetic floor on the least-squares line ==\n")
f <- lsq_prescribed(40, 12, seq(1, 4, length.out = 12), seed = 28)
A <- linop(f$A)
set.seed(29)
b <- matrix(stats::rnorm(40), 40, 1)
achieved <- lsqr(A, b, tol = 0, maxit = 200L)$certificate$values$backward_error
cat(sprintf("achieved backward error %.4e, c eps = %.4e\n",
            achieved, 4 * .Machine$double.eps))
wf <- lsqr(A, b, tol = achieved / 2, maxit = 200L)
wo <- lsqr(A, b, tol = achieved / 2, maxit = 200L, floor_const = 0)
cat(sprintf("identical iterates: %s\n", identical(wf$x, wo$x)))
for (nm in c("backward error", "convergence")) {
  cat(sprintf("%-16s with floor %-10s without %s\n", nm,
              linop:::cert_status(wf$certificate, nm),
              linop:::cert_status(wo$certificate, nm)))
}
