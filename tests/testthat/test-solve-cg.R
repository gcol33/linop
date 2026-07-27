## Gate 2 for the first of the seven methods: recovery against closed-form truth
## run to convergence, certificate coverage over 20 seeds, the rate its condition
## number predicts, and the refusals the section 4.3 table requires.

cg <- linop:::cg_solve
cert_status <- linop:::cert_status

## ------------------------------------------------------------------ recovery

test_that("cg recovers the solution of the 1-D Dirichlet Laplacian", {
  n <- 60
  A <- laplacian_1d(n)
  set.seed(1)
  x_true <- matrix(stats::rnorm(n), n, 1)
  b <- A %*% x_true

  fit <- cg(A, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$x, x_true, tolerance = 1e-9)
  expect_equal(fit$method, "cg")
  expect_equal(fit$restarts, 0L)
})

test_that("cg reproduces the closed-form KMS inverse, all columns at once", {
  n <- 50
  rho <- 0.7
  A <- as_spd_linop(kms_matrix(n, rho))
  fit <- cg(A, diag(1, n), tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$x, kms_inverse(n, rho), tolerance = 1e-7)
})

test_that("cg solves a complex hermitian positive definite system", {
  n <- 25
  A <- as_spd_linop(hpd_prescribed(n, seq(1, 20, length.out = n), seed = 4))
  set.seed(3)
  x_true <- matrix(complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)), n, 1)
  b <- A %*% x_true

  fit <- cg(A, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_true(is.complex(fit$x))
  expect_equal(fit$x, x_true, tolerance = 1e-8)
})

test_that("a vector goes in and a vector comes out", {
  A <- laplacian_1d(20)
  set.seed(4)
  x_true <- stats::rnorm(20)
  fit <- cg(A, as.numeric(A %*% x_true), tol = 1e-12)
  expect_null(dim(fit$x))
  expect_equal(fit$x, x_true, tolerance = 1e-9)
})

test_that("a zero right-hand side and a warm start both cost nothing", {
  A <- laplacian_1d(20)

  z <- cg(A, rep(0, 20))
  expect_equal(z$iterations, 0L)
  expect_true(all(z$x == 0))
  expect_true(z$converged)

  set.seed(4)
  x_true <- stats::rnorm(20)
  w <- cg(A, as.numeric(A %*% x_true), x0 = x_true, tol = 1e-8)
  expect_equal(w$iterations, 0L)
})

## ------------------------------------------------------------ Krylov identity

test_that("cg terminates in as many steps as the operator has distinct eigenvalues", {
  ## The Krylov statement, not a tolerance statement: with d distinct
  ## eigenvalues the CG polynomial of degree d annihilates the error exactly, so
  ## a run needing many more than d steps is wrong rather than slow.
  n <- 30
  A <- as_spd_linop(spd_prescribed(n, rep(c(1, 2, 3, 4, 5), each = 6), seed = 6))
  set.seed(2)
  b <- matrix(stats::rnorm(n), n, 1)
  fit <- cg(A, b, tol = 1e-10)
  expect_lte(fit$iterations, 6L)
  expect_true(fit$converged)
})

test_that("cg converges at least as fast as its condition number predicts", {
  n <- 60
  A <- laplacian_1d(n)
  lambda <- laplacian_1d_eigenvalues(n)
  kappa <- max(lambda) / min(lambda)
  rate <- (sqrt(kappa) - 1) / (sqrt(kappa) + 1)

  set.seed(7)
  b <- matrix(stats::rnorm(n), n, 1)
  fit <- cg(A, b, tol = 1e-12, history = TRUE)

  ## ||e_m||_A <= 2 rate^m ||e_0||_A, and ||r|| <= sqrt(lambda_max) ||e||_A with
  ## ||e_0||_A <= ||r_0|| / sqrt(lambda_min), so the residual obeys the same
  ## geometric rate with one factor of sqrt(kappa).
  ratio <- as.numeric(fit$history) / sqrt(sum(b^2))
  m <- seq_along(ratio)
  expect_true(all(ratio <= sqrt(kappa) * 2 * rate^m * (1 + 1e-8)))
  expect_lte(fit$iterations, 2L * n)
})

test_that("columns run in lockstep are the iterates of per-column cg", {
  n <- 40
  A <- as_spd_linop(spd_prescribed(n, exp(seq(log(0.5), log(50), length.out = n)), seed = 9))
  set.seed(11)
  B <- matrix(stats::rnorm(n * 4), n, 4)

  block <- cg(A, B, tol = 1e-11)
  for (j in seq_len(4)) {
    one <- cg(A, B[, j, drop = FALSE], tol = 1e-11)
    expect_identical(block$x[, j], one$x[, 1L],
                     info = sprintf("column %d diverged from its own solve", j))
  }
})

## ------------------------------------------------------ the arithmetic floor

test_that("a fully converged solve does not certify as fail, and the floor is what stops it", {
  ## S0.6: the bound keeps decaying past machine epsilon while the true error
  ## plateaus, so without a roundoff term the results that certify as fail are
  ## exactly the most converged ones. tol = 0 asks for an exact solve, which no
  ## floating-point iteration delivers.
  n <- 30
  A <- as_spd_linop(spd_prescribed(n, seq(1, 4, length.out = n), seed = 5))
  set.seed(2)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- cg(A, b, tol = 0, maxit = 200L)
  expect_equal(cert_status(fit$certificate, "residual"), "qualified")
  expect_false(fit$certificate$overall == "fail")

  bare <- cg(A, b, tol = 0, maxit = 200L, floor_const = 0)
  expect_equal(cert_status(bare$certificate, "residual"), "fail")
  expect_equal(bare$certificate$overall, "fail")

  ## the two ran the same iteration; only the certificate differs
  expect_identical(fit$x, bare$x)
})

test_that("a line that meets its tolerance only through the floor says estimate", {
  n <- 30
  A <- as_spd_linop(spd_prescribed(n, seq(1, 4, length.out = n), seed = 5))
  set.seed(2)
  b <- matrix(stats::rnorm(n), n, 1)

  clean <- cg(A, b, tol = 1e-8)$certificate
  row <- clean$checks[clean$checks$check == "residual", ]
  expect_equal(row$status, "pass")
  expect_equal(row$guarantee, "identity")

  floored <- cg(A, b, tol = 0, maxit = 200L)$certificate
  row <- floored$checks[floored$checks$check == "residual", ]
  expect_equal(row$guarantee, "estimate")
  expect_true(is.na(row$confidence))
})

test_that("the certificate reports the true residual, not the recurrence residual", {
  n <- 40
  M <- spd_prescribed(n, exp(seq(log(1e-2), log(1e2), length.out = n)), seed = 12)
  A <- as_spd_linop(M)
  set.seed(13)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- cg(A, b, tol = 1e-10)
  independent <- sqrt(sum((b - M %*% fit$x)^2)) / sqrt(sum(b^2))
  expect_equal(fit$certificate$values$residual, independent, tolerance = 1e-12)
})

test_that("forward error is not_checked, and the certificate says why", {
  A <- laplacian_1d(20)
  set.seed(8)
  cert <- cg(A, matrix(stats::rnorm(20), 20, 1), tol = 1e-10)$certificate
  row <- cert$checks[cert$checks$check == "forward error", ]
  expect_equal(row$status, "not_checked")
  expect_match(row$detail, "A\\^-1")
  expect_true("forward error" %in% cert$without_deterministic_bound)
  ## a certificate resting on an unchecked line is not a pass
  expect_equal(cert$overall, "qualified")
})

## -------------------------------------------------- certificate coverage, 20

test_that("over 20 seeds the certificate never understates what it reports", {
  n <- 40
  lambda <- exp(seq(log(1e-2), log(1e2), length.out = n))
  seeds <- seq_len(20)

  recovered <- logical(length(seeds))
  residual_agrees <- logical(length(seeds))
  bounds_true_omega <- logical(length(seeds))
  estimate_is_conservative <- logical(length(seeds))

  for (i in seq_along(seeds)) {
    s <- seeds[i]
    M <- spd_prescribed(n, lambda, seed = s)
    A <- as_spd_linop(M)
    set.seed(1000L + s)
    x_true <- matrix(stats::rnorm(n), n, 1)
    b <- M %*% x_true

    fit <- cg(A, b, tol = 1e-11)
    r <- b - M %*% fit$x
    sv <- svd(M, nu = 0L, nv = 0L)$d

    recovered[i] <- max(Mod(fit$x - x_true)) <= 1e-6 * max(Mod(x_true))
    residual_agrees[i] <- abs(fit$certificate$values$residual -
                              sqrt(sum(r^2)) / sqrt(sum(b^2))) <= 1e-12

    ## Rigal-Gaches, with the exact norm
    true_omega <- sqrt(sum(r^2)) /
      (sv[1L] * sqrt(sum(fit$x^2)) + sqrt(sum(b^2)))
    bounds_true_omega[i] <-
      fit$certificate$values$backward_error >= true_omega * (1 - 1e-9)

    ## Every route to ||A|| returns a lower bound, so an estimated norm shrinks
    ## the denominator and the reported backward error can only grow. Forcing the
    ## iteration route checks that the error is on the safe side.
    est <- cg(A, b, tol = 1e-11, norm_control = list(exact_max = 0))
    estimate_is_conservative[i] <-
      est$certificate$values$backward_error >=
      fit$certificate$values$backward_error * (1 - 1e-9)
  }

  failing <- function(ok) paste(seeds[!ok], collapse = ", ")
  expect_true(all(recovered), info = paste("seeds:", failing(recovered)))
  expect_true(all(residual_agrees), info = paste("seeds:", failing(residual_agrees)))
  expect_true(all(bounds_true_omega), info = paste("seeds:", failing(bounds_true_omega)))
  expect_true(all(estimate_is_conservative),
              info = paste("seeds:", failing(estimate_is_conservative)))
})

test_that("an exact norm and an estimated one are labelled differently", {
  n <- 40
  A <- as_spd_linop(spd_prescribed(n, seq(1, 9, length.out = n), seed = 3))
  set.seed(6)
  b <- matrix(stats::rnorm(n), n, 1)

  exact <- cg(A, b, tol = 1e-10)$certificate
  row <- exact$checks[exact$checks$check == "arithmetic floor", ]
  expect_equal(row$guarantee, "identity")
  expect_false("arithmetic floor" %in% exact$without_deterministic_bound)

  estimated <- cg(A, b, tol = 1e-10, norm_control = list(exact_max = 0))$certificate
  row <- estimated$checks[estimated$checks$check == "arithmetic floor", ]
  expect_equal(row$guarantee, "estimate")
  expect_true("arithmetic floor" %in% estimated$without_deterministic_bound)
})

## ------------------------------------------------------------- preconditioner

test_that("a jacobi preconditioner turns a system cg cannot finish into one it can", {
  n <- 60
  s <- 10^seq(-2, 2, length.out = n)
  M <- diag(s) %*% kms_matrix(n, 0.5) %*% diag(s)
  M <- (M + t(M)) / 2
  A <- as_spd_linop(M)
  sv <- svd(M, nu = 0L, nv = 0L)$d
  kappa <- sv[1L] / sv[n]
  expect_gt(kappa, 1e8)

  set.seed(21)
  b <- matrix(stats::rnorm(n), n, 1)
  tol <- 1e-12

  ## Unpreconditioned CG loses orthogonality long before it reaches this
  ## tolerance and stalls two digits short. The result is a failing certificate,
  ## returned rather than thrown: the iteration reports what it achieved.
  plain <- cg(A, b, tol = tol)
  expect_false(plain$converged)
  expect_equal(plain$iterations, 10L * n)
  expect_equal(cert_status(plain$certificate, "residual"), "fail")
  expect_equal(cert_status(plain$certificate, "convergence"), "fail")
  expect_equal(plain$certificate$overall, "fail")
  expect_match(plain$certificate$checks$detail[
    plain$certificate$checks$check == "convergence"], "budget exhausted")

  d <- diag(M)
  P <- preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                      positive_definite = TRUE, dim = c(n, n))
  prec <- cg(A, b, tol = tol, preconditioner = P)
  expect_true(prec$converged)
  expect_lt(prec$iterations, 100L)
  expect_lte(max(prec$certificate$values$residual), tol)

  ## Correct to what the residual permits and no further. At this kappa a
  ## relative residual of 1e-12 still leaves several digits of the solution free,
  ## which is the case the forward-error line is not_checked for.
  x_true <- solve(M, as.numeric(b))
  forward <- sqrt(sum((as.numeric(prec$x) - x_true)^2)) / sqrt(sum(x_true^2))
  expect_lt(forward, 100 * kappa * tol)
})

reference_cg <- function(M, b, tol, maxit) {
  x <- rep(0, length(b))
  r <- b
  p <- r
  rho <- sum(r * r)
  nb <- sqrt(sum(b^2))
  it <- 0L
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

test_that("the outer loop catches a recurrence that has drifted into believing itself", {
  ## At kappa 1e6 the recurrence residual declares the tolerance met while the
  ## true residual is fifty times larger. A solver that trusts its own recurrence
  ## reports a convergence it did not reach; recomputing b - A x is what stops
  ## that, and the restart is what the recomputation triggers.
  n <- 50
  M <- spd_prescribed(n, exp(seq(0, log(1e6), length.out = n)), seed = 1)
  A <- as_spd_linop(M)
  set.seed(7)
  b <- matrix(stats::rnorm(n), n, 1)
  tol <- 1e-12

  fit <- cg(A, b, tol = tol)
  expect_gte(fit$restarts, 1L)

  ref <- reference_cg(M, as.numeric(b), tol = tol, maxit = 500L)
  ref_true <- sqrt(sum((as.numeric(b) - M %*% ref$x)^2)) / sqrt(sum(b^2))
  expect_lt(ref$believed, tol)          # what a self-trusting recurrence reports
  expect_gt(ref_true, 10 * tol)         # what it actually achieved

  ## the reported number is the recomputed one, wherever that lands
  independent <- sqrt(sum((b - M %*% fit$x)^2)) / sqrt(sum(b^2))
  expect_equal(max(fit$certificate$values$residual), independent, tolerance = 1e-12)
  expect_gt(independent, tol)
})

test_that("the residual and the backward error can disagree, and both are right", {
  ## Same system. ||A|| ||x|| dwarfs ||b||, so the residual sits above the
  ## requested tolerance and below the arithmetic floor, while the backward error
  ## is at machine epsilon: x is the exact solution of a system indistinguishable
  ## from the one asked about. Reporting one number would have to pick a story.
  n <- 50
  A <- as_spd_linop(spd_prescribed(n, exp(seq(0, log(1e6), length.out = n)), seed = 1))
  set.seed(7)
  fit <- cg(A, matrix(stats::rnorm(n), n, 1), tol = 1e-12)
  v <- fit$certificate$values

  expect_gt(max(v$residual), 1e-12)
  expect_lt(max(v$residual), max(v$floor))
  expect_equal(cert_status(fit$certificate, "residual"), "qualified")

  expect_lt(max(v$backward_error), 1e-14)
  expect_equal(cert_status(fit$certificate, "backward error"), "pass")
  expect_true(fit$converged)
})

test_that("cg agrees with a textbook reference implementation, iterate for iterate", {
  ## The lockstep block, the shrinking active set and the true-residual restart
  ## are all supposed to leave the iterates untouched. A plain textbook
  ## recurrence on a dense matrix is the check that they do, at a condition
  ## number where no restart is triggered and the two paths must coincide.
  for (seed in 1:5) {
    M <- spd_prescribed(50, exp(seq(log(0.2), log(20), length.out = 50)), seed = seed)
    set.seed(500L + seed)
    b <- stats::rnorm(50)
    mine <- cg(as_spd_linop(M), b, tol = 1e-11)
    ref <- reference_cg(M, b, tol = 1e-11, maxit = 500L)
    expect_equal(mine$restarts, 0L, info = sprintf("seed %d", seed))
    expect_equal(mine$iterations, ref$iterations,
                 info = sprintf("seed %d", seed))
    ## bitwise, not merely close: the block path and the textbook path perform
    ## the same operations in the same order, and a gap here would mean one of
    ## the two structural departures had changed the recurrence after all
    expect_identical(mine$x, ref$x, info = sprintf("seed %d", seed))
  }
})

test_that("cg accepts a preconditioner on any of the three sides", {
  ## Section 4.3 leaves the CG row unrestricted. The mathematical coincidence
  ## behind it is Saad section 9.2, cited in the implementation; what this
  ## asserts is the contract, that none of the three is refused.
  n <- 20
  A <- laplacian_1d(n)
  set.seed(22)
  b <- matrix(stats::rnorm(n), n, 1)
  for (sd in c("left", "right", "split")) {
    P <- preconditioner(function(R) R / 2, fixed = TRUE, hermitian = TRUE,
                        positive_definite = TRUE, side = sd, dim = c(n, n))
    expect_silent(cg(A, b, tol = 1e-10, preconditioner = P))
  }
})

test_that("cg refuses a preconditioner whose declared properties it requires and lacks", {
  n <- 20
  A <- laplacian_1d(n)
  b <- rep(1, n)

  flexible <- preconditioner(function(R) R, fixed = FALSE, hermitian = TRUE,
                             positive_definite = TRUE)
  expect_error(cg(A, b, preconditioner = flexible), "fgmres")

  undeclared <- preconditioner(function(R) R)
  expect_error(cg(A, b, preconditioner = undeclared), "hermitian")

  half <- preconditioner(function(R) R, hermitian = TRUE)
  expect_error(cg(A, b, preconditioner = half), "positive_definite")
})

test_that("a preconditioner that returns the wrong shape is caught", {
  n <- 10
  A <- laplacian_1d(n)
  bad <- preconditioner(function(R) R[-1L, , drop = FALSE], fixed = TRUE,
                        hermitian = TRUE, positive_definite = TRUE)
  expect_error(cg(A, rep(1, n), preconditioner = bad), "returned a 9 x 1 block")
})

## ------------------------------------------------------------------- refusals

test_that("cg names the capability it needed", {
  n <- 10
  ## symmetric and hermitian are established exactly by the dense leaf;
  ## positive definiteness is not checked by anything, so it stays unknown
  A <- linop(spd_prescribed(n, seq(1, 2, length.out = n), seed = 1))
  expect_true(is.na(capv(A, "positive_definite")))
  expect_error(cg(A, rep(1, n)), "positive_definite")
  expect_error(cg(A, rep(1, n)), "Unknown is not false")
})

test_that("cg refuses a rectangular operator and says what a rectangular system is", {
  R <- linop(rmat(5, 4, seed = 1))
  expect_error(cg(R, rep(1, 4)), "square")
  expect_error(cg(R, rep(1, 4)), "least-squares")
})

test_that("cg contradicts a false positive_definite declaration by name", {
  ## p = e_1 on the first step, so p^H A p = -1 immediately: the refusal is
  ## structural rather than a matter of which seed was drawn.
  A <- linop_scaling(c(-1, 1, 1, 1), properties = c(positive_definite = TRUE))
  expect_error(cg(A, c(1, 0, 0, 0)), "declares positive_definite = TRUE")
  expect_error(cg(A, c(1, 0, 0, 0)), "contradicts it")
})

test_that("cg contradicts a false declaration on the preconditioner by name", {
  A <- laplacian_1d(4)
  P <- preconditioner(function(R) -R, fixed = TRUE, hermitian = TRUE,
                      positive_definite = TRUE)
  expect_error(cg(A, c(1, 0, 0, 0), preconditioner = P), "preconditioner declares")
})

test_that("a non-conformable right-hand side is refused", {
  A <- laplacian_1d(6)
  expect_error(cg(A, rep(1, 5)), "non-conformable")
  expect_error(cg(A, rep(1, 6), x0 = rep(0, 5)), "x0 is")
})

## ---------------------------------------------------------- evidence for auto

test_that("the requirement method = auto would apply admits proof and refuses declaration", {
  n <- 8
  req <- linop:::CG_PD_REQUIREMENT

  ## a bare declaration is the caller's claim, not an argument
  declared <- linop(function(X) X, adjoint = function(X) X, dim = c(n, n),
                    properties = c(hermitian = TRUE, positive_definite = TRUE))
  expect_false(evidence_satisfies(linop:::cape(declared, "positive_definite"), req))

  ## the signs of d prove it, on data the operator already holds
  proved <- linop_scaling(seq_len(n))
  expect_true(evidence_satisfies(linop:::cape(proved, "positive_definite"), req))

  ## and a probe never does, whatever its source
  expect_false(evidence_satisfies(linop:::ev_probe(), req))

  ## naming the method accepts the declaration; the requirement is for auto
  expect_silent(cg(declared, rep(1, n), tol = 1e-10))
})
