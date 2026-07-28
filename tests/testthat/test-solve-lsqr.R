## Gate 2 for the fifth of the seven methods, and the first whose problem is a
## minimisation rather than an equation: recovery against closed-form truth on
## compatible, incompatible and rank-deficient right-hand sides, agreement with
## the published recurrence, certificate coverage over 20 seeds, and the refusals
## the section 4.3 table requires now that the least-squares rows are narrowed.

lsqr <- linop:::lsqr_solve
cert_status <- linop:::cert_status

## ------------------------------------------------------------ the fixtures --
## A closed form is code, and a wrong one passes plausible tests, so the two
## rectangular fixtures are checked against their own definitions before any
## solver is run on them.

test_that("lsq_prescribed really has the singular value decomposition it claims", {
  for (dt in c("double", "complex")) {
    n <- 12
    sigma <- seq(1, 8, length.out = n)
    f <- lsq_prescribed(30, n, sigma, seed = 2, dtype = dt)

    ## A V = U diag(sigma) and A^H U = V diag(sigma), pairing included. A check
    ## against a sorted spectrum would pass with every vector mismatched.
    expect_lt(max(Mod(f$A %*% f$V - f$U %*% diag(sigma))), 1e-12)
    expect_lt(max(Mod(crossprod(Conj(f$A), f$U) - f$V %*% diag(sigma))), 1e-12)
    expect_lt(max(Mod(crossprod(Conj(f$U), f$U) - diag(1, n))), 1e-12)
    expect_lt(max(Mod(crossprod(Conj(f$V), f$V) - diag(1, n))), 1e-12)
  }
})

test_that("the closed-form least-squares solution and residual are what they say", {
  n <- 10
  f <- lsq_prescribed(40, n, seq(1, 5, length.out = n), seed = 7)
  set.seed(3)
  b <- matrix(stats::rnorm(40), 40, 1)
  x <- lsq_prescribed_solve(f, b)
  r <- lsq_prescribed_residual(f, b)

  ## the normal equations, and the residual as the part of b outside the range
  expect_lt(max(Mod(crossprod(f$A, b - f$A %*% x))), 1e-10)
  expect_lt(max(Mod(r - (b - f$A %*% x))), 1e-10)
  ## and it is a genuine least-squares problem rather than a compatible one
  expect_gt(sqrt(sum(r^2)), 1)
})

test_that("diff_1d is the first-difference operator with the adjoint it declares", {
  n <- 25
  D <- diff_1d(n)
  M <- as.matrix(D)
  expect_equal(dim(D), c(n - 1L, n))
  expect_equal(M, diff(diag(1, n)))
  ## the declared adjoint against the transpose of the materialised operator
  set.seed(5)
  Y <- matrix(stats::rnorm((n - 1L) * 3), n - 1L, 3)
  expect_equal(as.matrix(linop:::linop_apply(D, Y, "C")), t(M) %*% Y)
  ## its singular values, against the closed form
  expect_equal(sort(svd(M, nu = 0L, nv = 0L)$d),
               sort(diff_1d_singular_values(n))[-1L], tolerance = 1e-12)
})

test_that("diff_1d_min_norm_solve solves the system and is the one of least norm", {
  n <- 30
  D <- diff_1d(n)
  set.seed(6)
  b <- matrix(stats::rnorm(n - 1L), n - 1L, 1)
  x <- diff_1d_min_norm_solve(n, b)

  expect_lt(max(Mod(as.matrix(D %*% x) - b)), 1e-12)
  ## the nullspace is the constants, so least norm means orthogonal to them
  expect_lt(abs(sum(x)), 1e-12)
  ## every other solution is longer, by construction
  for (shift in c(-2, -0.5, 0.5, 2)) {
    expect_gt(sqrt(sum((x + shift)^2)), sqrt(sum(x^2)))
  }
})

## ------------------------------------------------------------------ recovery

test_that("lsqr recovers the least-squares solution of an incompatible system", {
  m <- 60; n <- 20
  f <- lsq_prescribed(m, n, seq(1, 10, length.out = n), seed = 3)
  A <- linop(f$A)
  set.seed(11)
  b <- matrix(stats::rnorm(m), m, 1)

  fit <- lsqr(A, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$method, "lsqr")
  expect_equal(fit$x, lsq_prescribed_solve(f, b), tolerance = 1e-9)
  ## and the residual it leaves is the distance from b to the range of A, which
  ## no method can reduce
  expect_equal(sqrt(sum((b - f$A %*% fit$x)^2)),
               sqrt(sum(lsq_prescribed_residual(f, b)^2)), tolerance = 1e-9)
})

test_that("lsqr recovers the solution of a compatible rectangular system", {
  m <- 50; n <- 18
  f <- lsq_prescribed(m, n, seq(1, 6, length.out = n), seed = 8)
  A <- linop(f$A)
  set.seed(12)
  x_true <- matrix(stats::rnorm(n), n, 1)

  fit <- lsqr(A, f$A %*% x_true, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$x, x_true, tolerance = 1e-9)
  expect_equal(cert_status(fit$certificate, "residual"), "pass")
})

test_that("lsqr solves a complex rectangular least-squares problem", {
  m <- 40; n <- 15
  f <- lsq_prescribed(m, n, seq(1, 7, length.out = n), seed = 4, dtype = "complex")
  A <- linop(f$A)
  set.seed(13)
  b <- matrix(complex(real = stats::rnorm(m), imaginary = stats::rnorm(m)), m, 1)

  fit <- lsqr(A, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_true(is.complex(fit$x))
  expect_equal(fit$x, lsq_prescribed_solve(f, b), tolerance = 1e-8)
})

test_that("lsqr returns the minimum-norm solution of a rank-deficient system", {
  ## D is onto, so every right-hand side is compatible and the solution set is a
  ## line. Started from zero the iterates stay in the range of D^H, which is the
  ## orthogonal complement of the nullspace, so the limit is the shortest one.
  n <- 40
  D <- diff_1d(n)
  set.seed(14)
  b <- matrix(stats::rnorm(n - 1L), n - 1L, 1)

  fit <- lsqr(D, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$x, diff_1d_min_norm_solve(n, b), tolerance = 1e-9)
  expect_lt(abs(sum(fit$x)), 1e-8)
})

test_that("lsqr solves a square indefinite system, where cg cannot", {
  n <- 50
  sigma <- 4 * sin(25 * pi / (2 * (n + 1)))^2
  A <- shifted_laplacian_1d(n, sigma + 0.35)
  set.seed(15)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- lsqr(A, b, tol = 1e-11, maxit = 4000L)
  expect_true(fit$converged)
  expect_equal(fit$x, shifted_laplacian_solve(n, sigma + 0.35, b), tolerance = 1e-6)
})

test_that("a vector goes in and a vector comes out", {
  f <- lsq_prescribed(30, 8, seq(1, 4, length.out = 8), seed = 16)
  A <- linop(f$A)
  set.seed(17)
  b <- stats::rnorm(30)
  fit <- lsqr(A, b, tol = 1e-12)
  expect_null(dim(fit$x))
  expect_equal(fit$x, as.numeric(lsq_prescribed_solve(f, b)), tolerance = 1e-9)
})

test_that("a zero right-hand side and a warm start both cost nothing", {
  f <- lsq_prescribed(30, 8, seq(1, 4, length.out = 8), seed = 18)
  A <- linop(f$A)

  z <- lsqr(A, rep(0, 30))
  expect_equal(z$iterations, 0L)
  expect_true(all(z$x == 0))
  expect_true(z$converged)

  set.seed(19)
  x_true <- matrix(stats::rnorm(8), 8, 1)
  w <- lsqr(A, f$A %*% x_true, x0 = x_true, tol = 1e-8)
  expect_equal(w$iterations, 0L)
})

## --------------------------------------------------- the published recurrence

## Paige and Saunders 1982, single column, no stopping test and no preconditioner:
## the bidiagonalisation and the rotations exactly as published, so agreement is a
## statement about the recurrence rather than about the bookkeeping around it.
reference_lsqr <- function(M, b, steps) {
  n <- ncol(M)
  Mh <- Conj(t(M))
  x <- rep(0, n)
  if (is.complex(M) || is.complex(b)) storage.mode(x) <- "complex"

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
    cs <- rhobar / rho
    sn <- beta / rho
    theta <- sn * alpha
    rhobar <- -cs * alpha
    phi <- cs * phibar
    phibar <- sn * phibar

    x <- x + (phi / rho) * w
    ## alpha = 0 is the exhausted Krylov space: the iterate above is already the
    ## least-squares solution and there is no direction left to normalise.
    if (alpha == 0) break
    v <- vt / alpha
    w <- v - (theta / rho) * w
  }
  list(x = matrix(x, n, 1L), iterations = steps)
}

test_that("the iteration is the published one, step for step", {
  ## tol = 0 removes every stopping test, so both run the same fixed number of
  ## steps and any gap is the recurrence itself.
  ##
  ## Four steps, and the tolerance below cannot be carried much further out. Two
  ## arithmetically equivalent LSQR implementations do not stay close: the
  ## relative gap here grows about three orders of magnitude every two steps and
  ## reaches 1e-2 by step 12 of a 15-column problem, because the bidiagonal
  ## vectors lose orthogonality and what fills the lost direction is whichever
  ## rounding noise arrived first. It closes again as both converge on the same
  ## least-squares solution, which the recovery tests above check. The numbers
  ## are in dev_notes/lsqr-and-the-least-squares-certificate.md, and they are the
  ## reason the certificate recomputes rather than reporting the recurrence.
  for (seed in 1:5) {
    f <- lsq_prescribed(45, 15, exp(seq(log(0.5), log(20), length.out = 15)),
                        seed = seed)
    set.seed(600L + seed)
    b <- matrix(stats::rnorm(45), 45, 1)
    mine <- lsqr(linop(f$A), b, tol = 0, maxit = 4L)
    ref <- reference_lsqr(f$A, b, 4L)$x
    expect_equal(mine$iterations, 4L, info = sprintf("seed %d", seed))
    expect_lt(max(Mod(mine$x - ref)) / max(Mod(ref)), 1e-13,
              label = sprintf("seed %d", seed))
  }
})

test_that("a solve whose recurrence has stopped being reproducible still lands on the answer", {
  ## The other half of the same fact. Run to convergence rather than to a fixed
  ## step count, and the two implementations agree to the accuracy asked for,
  ## because what they agree on is the least-squares solution rather than the
  ## path to it.
  for (seed in 1:5) {
    f <- lsq_prescribed(45, 15, exp(seq(log(0.5), log(20), length.out = 15)),
                        seed = seed)
    set.seed(600L + seed)
    b <- matrix(stats::rnorm(45), 45, 1)
    truth <- lsq_prescribed_solve(f, b)
    mine <- lsqr(linop(f$A), b, tol = 1e-12, maxit = 500L)
    ## 30 steps on a 15-column operator, and the margin is not generosity. In
    ## exact arithmetic the Krylov space is complete at step 15 and the iterate
    ## there is the answer; in floating point the bidiagonal vectors have lost
    ## orthogonality by then and the published recurrence is still wrong in the
    ## second digit, reaching 1e-9 only around step 20 and 1e-15 around step 25.
    ## Finite termination is the property this method loses first, which is why
    ## the iteration budget defaults to a multiple of n rather than to n.
    ref <- reference_lsqr(f$A, b, 30L)$x
    expect_lt(max(Mod(mine$x - truth)) / max(Mod(truth)), 1e-9,
              label = sprintf("seed %d, this implementation", seed))
    expect_lt(max(Mod(ref - truth)) / max(Mod(truth)), 1e-9,
              label = sprintf("seed %d, the published one", seed))
  }
})

test_that("columns run in lockstep are the iterates of per-column lsqr", {
  f <- lsq_prescribed(50, 16, exp(seq(log(0.3), log(30), length.out = 16)), seed = 9)
  A <- linop(f$A)
  set.seed(21)
  B <- matrix(stats::rnorm(50 * 4), 50, 4)

  block <- lsqr(A, B, tol = 1e-11)
  for (j in seq_len(4)) {
    one <- lsqr(A, B[, j, drop = FALSE], tol = 1e-11)
    expect_identical(block$x[, j], one$x[, 1L],
                     info = sprintf("column %d diverged from its own solve", j))
  }
})

test_that("a block whose columns converge at different rates still runs in lockstep", {
  ## One column inside the range of A and one outside it stop on different tests
  ## at different steps, which is the case the shrinking active set exists for.
  f <- lsq_prescribed(40, 12, seq(1, 9, length.out = 12), seed = 22)
  A <- linop(f$A)
  set.seed(23)
  B <- cbind(f$A %*% matrix(stats::rnorm(12), 12, 1),
             matrix(stats::rnorm(40), 40, 1))

  block <- lsqr(A, B, tol = 1e-11)
  for (j in 1:2) {
    one <- lsqr(A, B[, j, drop = FALSE], tol = 1e-11)
    expect_identical(block$x[, j], one$x[, 1L], info = sprintf("column %d", j))
  }
})

## --------------------------------------------------------------- certificate

test_that("an incompatible system does not fail the residual line", {
  ## The whole reason the certificate has two readings. A converged
  ## least-squares solve leaves a residual that no x can reduce, and testing it
  ## against tol would report the correct answer as a failure.
  m <- 50; n <- 12
  f <- lsq_prescribed(m, n, seq(1, 5, length.out = n), seed = 24)
  set.seed(25)
  b <- matrix(stats::rnorm(m), m, 1)

  cert <- lsqr(linop(f$A), b, tol = 1e-12)$certificate
  expect_equal(cert_status(cert, "residual"), "not_checked")
  expect_equal(cert_status(cert, "backward error"), "pass")
  expect_equal(cert_status(cert, "convergence"), "pass")
  expect_false(cert$overall == "fail")
  ## and it says why, rather than leaving a reader to infer it
  detail <- cert$checks$detail[cert$checks$check == "residual"]
  expect_match(detail, "not in the range of A")
})

test_that("the reported backward error is one an exhibited perturbation achieves", {
  ## Stewart's dA = -(r r^H A) / ||r||^2 makes x the exact least-squares solution
  ## of A + dA, so the number on the certificate is a distance to a problem that
  ## exists rather than a bound on one. Checking that perturbation directly is
  ## what separates the two.
  m <- 45; n <- 14
  f <- lsq_prescribed(m, n, seq(1, 6, length.out = n), seed = 26)
  M <- f$A
  set.seed(27)
  b <- matrix(stats::rnorm(m), m, 1)

  fit <- lsqr(linop(M), b, tol = 1e-10)
  x <- fit$x
  r <- b - M %*% x
  dA <- -(r %*% (crossprod(Conj(r), M))) / sum(Mod(r)^2)
  Ap <- M + dA

  ## x satisfies the normal equations of the perturbed problem exactly
  expect_lt(max(Mod(crossprod(Conj(Ap), b - Ap %*% x))),
            1e-8 * max(Mod(crossprod(Conj(M), b))))
  ## and the size of that perturbation is the number the certificate reported
  rel_dA <- svd(dA, nu = 0L, nv = 0L)$d[1L] / svd(M, nu = 0L, nv = 0L)$d[1L]
  expect_lt(rel_dA, fit$certificate$values$backward_error * (1 + 1e-6))
})

test_that("over 20 seeds the certificate never understates what it reports", {
  n <- 16
  sigma <- exp(seq(log(1e-1), log(1e2), length.out = n))
  seeds <- seq_len(20)

  recovered <- logical(length(seeds))
  optimality_agrees <- logical(length(seeds))
  norm_is_a_lower_bound <- logical(length(seeds))
  bounds_true_backward_error <- logical(length(seeds))
  estimate_is_conservative <- logical(length(seeds))

  for (i in seq_along(seeds)) {
    s <- seeds[i]
    f <- lsq_prescribed(50, n, sigma, seed = s)
    M <- f$A
    A <- linop(M)
    set.seed(2000L + s)
    b <- matrix(stats::rnorm(50), 50, 1)

    fit <- lsqr(A, b, tol = 1e-11, maxit = 2000L)
    r <- b - M %*% fit$x
    atr <- crossprod(Conj(M), r)
    sv <- svd(M, nu = 0L, nv = 0L)$d

    recovered[i] <- max(Mod(fit$x - lsq_prescribed_solve(f, b))) <=
      1e-6 * max(Mod(lsq_prescribed_solve(f, b)))

    ## the certificate recomputes both quantities itself, so they have to agree
    ## with the same quantities formed densely
    optimality_agrees[i] <-
      abs(fit$certificate$values$optimality -
          sqrt(sum(Mod(atr)^2)) / (fit$certificate$values$norm * sqrt(sum(Mod(r)^2)))) <= 1e-12

    ## every route to ||A|| returns a lower bound
    norm_is_a_lower_bound[i] <- fit$certificate$values$norm <= sv[1L] * (1 + 1e-9)

    ## with the exact norm the backward error is smaller, so the reported one
    ## overstates rather than understates
    true_bw <- sqrt(sum(Mod(atr)^2)) / (sv[1L] * sqrt(sum(Mod(r)^2)))
    bounds_true_backward_error[i] <-
      fit$certificate$values$backward_error >= true_bw * (1 - 1e-9)

    ## and forcing the iteration route to ||A|| can only push it further that way
    est <- lsqr(A, b, tol = 1e-11, maxit = 2000L, norm_control = list(exact_max = 0))
    estimate_is_conservative[i] <-
      est$certificate$values$backward_error >=
      fit$certificate$values$backward_error * (1 - 1e-9)
  }

  expect_true(all(recovered))
  expect_true(all(optimality_agrees))
  expect_true(all(norm_is_a_lower_bound))
  expect_true(all(bounds_true_backward_error))
  expect_true(all(estimate_is_conservative))
})

test_that("the arithmetic floor is load-bearing on the least-squares line too", {
  ## S0.6 for the residual, and the same statement for the optimality ratio. The
  ## tolerance is set from what the solve actually delivered rather than guessed,
  ## because the window the floor covers is one c eps wide and a fixed tolerance
  ## either sits inside it or does not depending on the fixture. Halving the
  ## achieved value puts the request strictly below what the arithmetic can
  ## reach, which is the regime the floor exists for, and both runs are then the
  ## same iterates certified twice.
  f <- lsq_prescribed(40, 12, seq(1, 4, length.out = 12), seed = 28)
  A <- linop(f$A)
  set.seed(29)
  b <- matrix(stats::rnorm(40), 40, 1)

  achieved <- lsqr(A, b, tol = 0, maxit = 200L)$certificate$values$backward_error
  expect_lt(achieved, 2 * linop:::SOLVE_FLOOR_CONST * .Machine$double.eps)

  with_floor <- lsqr(A, b, tol = achieved / 2, maxit = 200L)
  without <- lsqr(A, b, tol = achieved / 2, maxit = 200L, floor_const = 0)

  expect_identical(with_floor$x, without$x)
  expect_equal(cert_status(with_floor$certificate, "backward error"), "qualified")
  expect_equal(cert_status(with_floor$certificate, "convergence"), "pass")
  expect_equal(cert_status(without$certificate, "backward error"), "fail")
  expect_equal(cert_status(without$certificate, "convergence"), "fail")
})

test_that("a spent budget comes back as a certificate rather than as an error", {
  f <- lsq_prescribed(60, 40, exp(seq(log(1e-4), log(1), length.out = 40)), seed = 30)
  A <- linop(f$A)
  set.seed(31)
  b <- matrix(stats::rnorm(60), 60, 1)

  fit <- lsqr(A, b, tol = 1e-14, maxit = 3L)
  expect_false(fit$converged)
  expect_equal(cert_status(fit$certificate, "convergence"), "fail")
  expect_match(fit$certificate$checks$detail[fit$certificate$checks$check == "convergence"],
               "budget exhausted")
})

## -------------------------------------------------------------- the operator

test_that("lsqr needs an adjoint, which is where it parts company with gmres", {
  ## Every step applies A^H, so an operator supplying only a forward action is
  ## solvable by gmres and by nothing else in the package.
  A <- linop(function(X) rbind(X, 0), dim = c(11L, 10L))
  expect_error(lsqr(A, rep(1, 11)), "adjoint")
})

test_that("a non-conformable right-hand side is refused", {
  f <- lsq_prescribed(20, 6, seq(1, 3, length.out = 6), seed = 32)
  A <- linop(f$A)
  expect_error(lsqr(A, rep(1, 19)), "non-conformable")
  expect_error(lsqr(A, rep(1, 20), x0 = rep(0, 5)), "x0 is")
  expect_error(lsqr(A, rep(1, 20), maxit = 0L), "maxit")
  expect_error(lsqr(matrix(1, 3, 2), rep(1, 3)), "expects a linop")
})

test_that("the square methods refuse the shape lsqr exists for", {
  f <- lsq_prescribed(20, 6, seq(1, 3, length.out = 6), seed = 33)
  A <- linop(f$A)
  b <- rep(1, 20)
  expect_error(linop:::cg_solve(A, b), "square")
  expect_error(linop:::minres_solve(A, b), "square")
  expect_error(linop:::gmres_solve(A, b), "square")
  expect_silent(lsqr(A, b, maxit = 5L))
})

## ---------------------------------------------------------- preconditioning

test_that("left and split preconditioners are refused, and the message says why", {
  f <- lsq_prescribed(20, 6, seq(1, 3, length.out = 6), seed = 34)
  A <- linop(f$A)
  for (sd in c("left", "split")) {
    P <- preconditioner(function(R) R, side = sd, hermitian = TRUE)
    expect_error(lsqr(A, rep(1, 20), preconditioner = P),
                 "acts on the domain of A")
  }
})

test_that("a right preconditioner without M^-H is refused by name", {
  ## The iteration runs on A M^-1, whose adjoint is M^-H A^H, so this is a
  ## requirement of the method rather than a preference. A hermitian M supplies
  ## it for nothing and is accepted.
  f <- lsq_prescribed(20, 6, seq(1, 3, length.out = 6), seed = 35)
  A <- linop(f$A)
  b <- rep(1, 20)

  bare <- preconditioner(function(R) R, side = "right")
  expect_error(lsqr(A, b, preconditioner = bare), "M\\^-H")

  herm <- preconditioner(function(R) R, side = "right", hermitian = TRUE)
  expect_silent(lsqr(A, b, preconditioner = herm, maxit = 5L))

  supplied <- preconditioner(function(R) R, apply_inverse_adjoint = function(R) R,
                             side = "right")
  expect_silent(lsqr(A, b, preconditioner = supplied, maxit = 5L))
})

test_that("column equilibration converges where the unscaled operator does not", {
  ## A right preconditioner on a least-squares problem is a change of variable in
  ## the domain, so the test is that it changes the iterates and that the change
  ## survives more than one seed at the same budget and the same tolerance. The
  ## fixture is an operator whose columns differ by 1e6 in norm, where D^-1 is
  ## exactly the scaling that undoes it.
  n <- 20
  d <- 10^seq(-3, 3, length.out = n)
  f <- lsq_prescribed(60, n, seq(1, 4, length.out = n), seed = 36)
  M <- f$A %*% diag(d)
  A <- linop(M)
  P <- preconditioner(function(R) R / d, side = "right", hermitian = TRUE)

  helped <- logical(6)
  differs <- logical(6)
  for (s in 1:6) {
    set.seed(700L + s)
    b <- matrix(stats::rnorm(60), 60, 1)
    plain <- lsqr(A, b, tol = 1e-10, maxit = 40L)
    prec <- lsqr(A, b, tol = 1e-10, preconditioner = P, maxit = 40L)
    helped[s] <- prec$converged && !plain$converged
    differs[s] <- !isTRUE(all.equal(prec$x, plain$x))
  }
  expect_true(all(differs))
  expect_true(all(helped))
})

test_that("a preconditioned solve lands on the same answer as an unpreconditioned one", {
  ## The preconditioner is a change of variable and must not move the minimiser.
  n <- 14
  d <- 10^seq(-2, 2, length.out = n)
  f <- lsq_prescribed(50, n, seq(1, 5, length.out = n), seed = 37)
  M <- f$A %*% diag(d)
  A <- linop(M)
  P <- preconditioner(function(R) R / d, side = "right", hermitian = TRUE)
  set.seed(38)
  b <- matrix(stats::rnorm(50), 50, 1)

  prec <- lsqr(A, b, tol = 1e-12, preconditioner = P, maxit = 2000L)
  truth <- diag(1 / d) %*% lsq_prescribed_solve(f, b)
  expect_true(prec$converged)
  expect_equal(prec$x, truth, tolerance = 1e-7)
})

## ------------------------------------------------------------ the condition limit

test_that("the condition limit stops a bidiagonal that has stopped meaning anything", {
  ## The same failure the gmres suite records, reached through the other
  ## projected problem: past the point where the bidiagonal is numerically
  ## singular the recurrence estimate keeps falling while the true solution runs
  ## away.
  n <- 30
  M <- outer(seq(0.1, 2, length.out = 50), seq_len(n) - 1L, "^")
  A <- linop(M)

  for (s in 1:4) {
    set.seed(800L + s)
    b <- matrix(stats::rnorm(50), 50, 1)
    run <- function(cl) lsqr(A, b, tol = 1e-300, maxit = 30L, conlim = cl)
    limited <- run(linop:::KRYLOV_CONDITION_LIMIT)
    unlimited <- run(Inf)
    expect_lte(sqrt(sum(limited$x^2)), sqrt(sum(unlimited$x^2)),
               label = sprintf("seed %d", s))
  }
})

test_that("the condition limit never fires on a well-posed problem", {
  ## A stopping rule that ends an iteration has to be silent where it is not
  ## needed, so where the bidiagonal stays healthy the two settings agree bit for
  ## bit rather than merely closely.
  fixtures <- list(
    tall = lsq_prescribed(60, 20, seq(1, 10, length.out = 20), seed = 40)$A,
    clustered = lsq_prescribed(50, 15, rep(c(1, 2, 3), each = 5), seed = 41)$A,
    spread = lsq_prescribed(70, 25, exp(seq(log(0.05), log(20), length.out = 25)),
                            seed = 42)$A)
  for (nm in names(fixtures)) {
    M <- fixtures[[nm]]
    A <- linop(M)
    set.seed(43)
    b <- matrix(stats::rnorm(nrow(M)), nrow(M), 1)
    a <- lsqr(A, b, tol = 1e-11, maxit = 1200L)
    z <- lsqr(A, b, tol = 1e-11, maxit = 1200L, conlim = Inf)
    expect_identical(a$x, z$x, info = nm)
    expect_equal(a$iterations, z$iterations, info = nm)
  }
})

## ------------------------------------------------------------------- history

test_that("history records one row per iteration and one column per right-hand side", {
  f <- lsq_prescribed(40, 10, seq(1, 5, length.out = 10), seed = 44)
  A <- linop(f$A)
  set.seed(45)
  B <- matrix(stats::rnorm(40 * 3), 40, 3)
  fit <- lsqr(A, B, tol = 1e-10, history = TRUE)
  expect_equal(ncol(fit$history), 3L)
  expect_equal(nrow(fit$history), fit$iterations)
  expect_null(lsqr(A, B, tol = 1e-10)$history)
})
