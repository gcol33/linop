## Reproduces every number in dev_notes/cg-and-the-arithmetic-floor.md.
##   Rscript dev_notes/spikes/cg-floor.R

pkg <- "C:/Users/Gilles Colling/Documents/dev/linop"
suppressMessages(library(devtools))
devtools::load_all(pkg, quiet = TRUE)
source(file.path(pkg, "tests/testthat/helper-linop.R"))

cg <- linop:::cg_solve
norm2 <- linop:::norm2
cert_status <- linop:::cert_status

cat("R:", R.version.string, "\n\n")

## --------------------------------------------------------------- recovery ---
cat("=== recovery on the 1-D Dirichlet Laplacian, tol = 1e-12 ===\n")
for (n in c(20, 60, 200)) {
  A <- laplacian_1d(n)
  set.seed(1)
  x_true <- matrix(stats::rnorm(n), n, 1)
  fit <- cg(A, A %*% x_true, tol = 1e-12)
  lambda <- laplacian_1d_eigenvalues(n)
  cat(sprintf("n = %3d | kappa %9.3g | %4d it | rel resid %.2e | max |x-x*| %.2e\n",
              n, max(lambda) / min(lambda), fit$iterations,
              max(fit$certificate$values$residual), max(abs(fit$x - x_true))))
}

## ------------------------------------------------------ the Krylov statement -
cat("\n=== distinct eigenvalues bound the iteration count ===\n")
for (d in c(3, 5, 10)) {
  n <- 30
  A <- as_spd_linop(spd_prescribed(n, rep(seq_len(d), each = n / d), seed = 6))
  set.seed(2)
  fit <- cg(A, matrix(stats::rnorm(n), n, 1), tol = 1e-10)
  cat(sprintf("%2d distinct eigenvalues | %2d iterations\n", d, fit$iterations))
}

## ------------------------------------------------------- the arithmetic floor
cat("\n=== the arithmetic floor, tol = 0 ===\n")
n <- 30
A <- as_spd_linop(spd_prescribed(n, seq(1, 4, length.out = n), seed = 5))
set.seed(2)
b <- matrix(stats::rnorm(n), n, 1)
with_floor <- cg(A, b, tol = 0, maxit = 200L)
no_floor <- cg(A, b, tol = 0, maxit = 200L, floor_const = 0)
cat(sprintf("relative residual achieved  %.3e\n",
            max(with_floor$certificate$values$residual)))
cat(sprintf("relative floor at c = 4     %.3e\n",
            max(with_floor$certificate$values$floor)))
cat(sprintf("c = 4 -> residual %s, overall %s\n",
            cert_status(with_floor$certificate, "residual"),
            with_floor$certificate$overall))
cat(sprintf("c = 0 -> residual %s, overall %s\n",
            cert_status(no_floor$certificate, "residual"),
            no_floor$certificate$overall))
cat(sprintf("iterates bitwise identical: %s\n", identical(with_floor$x, no_floor$x)))

## ------------------------------------------------------------- norm routes ---
cat("\n=== norm routes ===\n")
operators <- list(
  "linop_eye(9)"          = linop_eye(9),
  "scaling(2,-5,3)"       = linop_scaling(c(2, -5, 3)),
  "3 * scaling(2,-5,3)"   = 3 * linop_scaling(c(2, -5, 3)),
  "dense 20 x 20"         = linop(rmat(20, 20, seed = 3)),
  "laplacian n = 200"     = laplacian_1d(200))
for (nm in names(operators)) {
  est <- norm2(operators[[nm]])
  cat(sprintf("%-22s %-10s %10.6f  exact: %s\n", nm, est$method, est$value,
              evidence_satisfies(est$evidence, requirement(guarantees = "identity"))))
}
cat(sprintf("%-22s %-10s %10.6f\n", "laplacian n = 200", "truth",
            max(laplacian_1d_eigenvalues(200))))

## an estimate under a structural rule stays an estimate
scaled <- norm2(4 * laplacian_1d(200))
cat(sprintf("4 * laplacian: %s, exact at the top: %s\n",
            format_evidence(scaled$evidence),
            evidence_satisfies(scaled$evidence, requirement(guarantees = "identity"))))

## ------------------------------------------------------ preconditioning ------
cat("\n=== jacobi on a badly scaled system ===\n")
n <- 60
s <- 10^seq(-2, 2, length.out = n)
M <- diag(s) %*% kms_matrix(n, 0.5) %*% diag(s)
M <- (M + t(M)) / 2
A <- as_spd_linop(M)
d <- diag(M)
P <- preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                    positive_definite = TRUE, dim = c(n, n))
sv <- svd(M, nu = 0L, nv = 0L)$d
set.seed(21)
b <- matrix(stats::rnorm(n), n, 1)
plain <- cg(A, b, tol = 1e-12)
prec <- cg(A, b, tol = 1e-12, preconditioner = P)
cat(sprintf("kappa %.3g\n", sv[1L] / sv[n]))
cat(sprintf("unpreconditioned %4d it, resid %.3e, certificate %s\n",
            plain$iterations, max(plain$certificate$values$residual),
            plain$certificate$overall))
cat(sprintf("jacobi           %4d it, resid %.3e, certificate %s\n",
            prec$iterations, max(prec$certificate$values$residual),
            prec$certificate$overall))

## ------------------------------------------------------------- lockstep ------
cat("\n=== lockstep against per-column, 8 right-hand sides ===\n")
n <- 40
A <- as_spd_linop(spd_prescribed(n, exp(seq(log(0.5), log(50), length.out = n)), seed = 9))
set.seed(11)
B <- matrix(stats::rnorm(n * 8), n, 8)
block <- cg(A, B, tol = 1e-11)
singles <- lapply(seq_len(8), function(j) cg(A, B[, j, drop = FALSE], tol = 1e-11))
same <- vapply(seq_len(8), function(j) identical(block$x[, j], singles[[j]]$x[, 1L]),
               logical(1))
cat(sprintf("bitwise identical columns: %d of 8\n", sum(same)))
cat(sprintf("block applies %d, per-column applies would be %d\n",
            block$iterations, sum(vapply(singles, function(f) f$iterations, integer(1)))))

## ------------------------------------ drift, and the two lines that disagree -
cat("\n=== the recurrence believing itself, kappa 1e6 ===\n")
reference_believing <- function(M, b, tol, maxit) {
  x <- rep(0, length(b)); r <- b; p <- r
  rho <- sum(r * r); nb <- sqrt(sum(b^2)); it <- 0L
  while (sqrt(rho) > tol * nb && it < maxit) {
    it <- it + 1L
    q <- as.numeric(M %*% p)
    alpha <- rho / sum(p * q)
    x <- x + alpha * p
    r <- r - alpha * q
    rho_new <- sum(r * r)
    p <- r + (rho_new / rho) * p
    rho <- rho_new
  }
  list(x = x, iterations = it, believed = sqrt(rho) / nb)
}
n <- 50
M <- spd_prescribed(n, exp(seq(0, log(1e6), length.out = n)), seed = 1)
A <- as_spd_linop(M)
set.seed(7)
b <- matrix(stats::rnorm(n), n, 1)
mine <- cg(A, b, tol = 1e-12)
ref <- reference_believing(M, as.numeric(b), tol = 1e-12, maxit = 500L)
ref_true <- sqrt(sum((as.numeric(b) - M %*% ref$x)^2)) / sqrt(sum(b^2))
cat(sprintf("textbook, self-trusting | %3d it | reports %.3e | truly %.3e\n",
            ref$iterations, ref$believed, ref_true))
cat(sprintf("this implementation     | %3d it (%d restart) | reports %.3e | truly %.3e\n",
            mine$iterations, mine$restarts,
            max(mine$certificate$values$residual),
            sqrt(sum((b - M %*% mine$x)^2)) / sqrt(sum(b^2))))
checks <- mine$certificate$checks
cat("\n", checks$detail[checks$check == "arithmetic floor"], "\n", sep = "")
cat(sprintf("residual %.4e, backward error %.4e, converged %s\n\n",
            max(mine$certificate$values$residual),
            max(mine$certificate$values$backward_error),
            mine$converged))
print(mine$certificate)

## --------------------------------------------------- reference agreement -----
cat("\n=== agreement with a textbook recurrence ===\n")
reference_cg <- function(M, b, tol, maxit) {
  x <- rep(0, length(b)); r <- b; p <- r
  rho <- sum(r * r); nb <- sqrt(sum(b^2)); it <- 0L
  while (sqrt(rho) > tol * nb && it < maxit) {
    it <- it + 1L
    q <- as.numeric(M %*% p)
    alpha <- rho / sum(p * q)
    x <- x + alpha * p
    r <- r - alpha * q
    rho_new <- sum(r * r)
    p <- r + (rho_new / rho) * p
    rho <- rho_new
  }
  list(x = x, iterations = it)
}
for (seed in 1:5) {
  M <- spd_prescribed(50, exp(seq(log(0.2), log(20), length.out = 50)), seed = seed)
  set.seed(500L + seed)
  b <- stats::rnorm(50)
  mine <- cg(as_spd_linop(M), b, tol = 1e-11)
  ref <- reference_cg(M, b, tol = 1e-11, maxit = 500L)
  cat(sprintf("seed %d | mine %3d it, reference %3d it | max gap %.2e\n",
              seed, mine$iterations, ref$iterations, max(abs(mine$x - ref$x))))
}
