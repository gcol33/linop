## Gate 2 for the sixth of the seven methods: recovery against closed-form truth
## on compatible, incompatible and rank-deficient right-hand sides, agreement
## with the published recurrence including on ill-conditioned problems,
## certificate coverage over 20 seeds, the refusals the section 4.3 table
## requires, and the one property that distinguishes LSMR from LSQR rather than
## being inherited from it.
##
## The rectangular fixtures are checked against their own definitions in
## test-solve-lsqr.R, which runs first, so they are used here rather than
## re-verified.

lsmr <- linop:::lsmr_solve
lsqr <- linop:::lsqr_solve
cert_status <- linop:::cert_status

## ------------------------------------------------------------------ recovery

test_that("lsmr recovers the least-squares solution of an incompatible system", {
  m <- 60; n <- 20
  f <- lsq_prescribed(m, n, seq(1, 10, length.out = n), seed = 3)
  A <- linop(f$A)
  set.seed(11)
  b <- matrix(stats::rnorm(m), m, 1)

  fit <- lsmr(A, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$method, "lsmr")
  expect_equal(fit$x, lsq_prescribed_solve(f, b), tolerance = 1e-9)
  ## and the residual it leaves is the distance from b to the range of A, which
  ## no method can reduce
  expect_equal(sqrt(sum((b - f$A %*% fit$x)^2)),
               sqrt(sum(lsq_prescribed_residual(f, b)^2)), tolerance = 1e-9)
})

test_that("lsmr recovers the solution of a compatible rectangular system", {
  m <- 50; n <- 18
  f <- lsq_prescribed(m, n, seq(1, 6, length.out = n), seed = 8)
  A <- linop(f$A)
  set.seed(12)
  x_true <- matrix(stats::rnorm(n), n, 1)

  fit <- lsmr(A, f$A %*% x_true, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$x, x_true, tolerance = 1e-9)
  expect_equal(cert_status(fit$certificate, "residual"), "pass")
})

test_that("lsmr solves a complex rectangular least-squares problem", {
  m <- 40; n <- 15
  f <- lsq_prescribed(m, n, seq(1, 7, length.out = n), seed = 4, dtype = "complex")
  A <- linop(f$A)
  set.seed(13)
  b <- matrix(complex(real = stats::rnorm(m), imaginary = stats::rnorm(m)), m, 1)

  fit <- lsmr(A, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_true(is.complex(fit$x))
  expect_equal(fit$x, lsq_prescribed_solve(f, b), tolerance = 1e-8)
})

test_that("lsmr returns the minimum-norm solution of a rank-deficient system", {
  ## D is onto, so every right-hand side is compatible and the solution set is a
  ## line. Started from zero the iterates stay in the range of D^H, which is the
  ## orthogonal complement of the nullspace, so the limit is the shortest one.
  n <- 40
  D <- diff_1d(n)
  set.seed(14)
  b <- matrix(stats::rnorm(n - 1L), n - 1L, 1)

  fit <- lsmr(D, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$x, diff_1d_min_norm_solve(n, b), tolerance = 1e-9)
  expect_lt(abs(sum(fit$x)), 1e-8)
})

test_that("lsmr solves a square indefinite system, where cg cannot", {
  n <- 50
  sigma <- 4 * sin(25 * pi / (2 * (n + 1)))^2
  A <- shifted_laplacian_1d(n, sigma + 0.35)
  set.seed(15)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- lsmr(A, b, tol = 1e-11, maxit = 4000L)
  expect_true(fit$converged)
  expect_equal(fit$x, shifted_laplacian_solve(n, sigma + 0.35, b), tolerance = 1e-6)
})

test_that("a vector goes in and a vector comes out", {
  f <- lsq_prescribed(30, 8, seq(1, 4, length.out = 8), seed = 16)
  A <- linop(f$A)
  set.seed(17)
  b <- stats::rnorm(30)
  fit <- lsmr(A, b, tol = 1e-12)
  expect_null(dim(fit$x))
  expect_equal(fit$x, as.numeric(lsq_prescribed_solve(f, b)), tolerance = 1e-9)
})

test_that("a zero right-hand side and a warm start both cost nothing", {
  f <- lsq_prescribed(30, 8, seq(1, 4, length.out = 8), seed = 18)
  A <- linop(f$A)

  z <- lsmr(A, rep(0, 30))
  expect_equal(z$iterations, 0L)
  expect_true(all(z$x == 0))
  expect_true(z$converged)

  set.seed(19)
  x_true <- matrix(stats::rnorm(8), 8, 1)
  w <- lsmr(A, f$A %*% x_true, x0 = x_true, tol = 1e-8)
  expect_equal(w$iterations, 0L)
})

test_that("a Krylov space that closes at the first step is solved there", {
  ## b = A v_1 for a single right singular vector, so A^H u_2 already lies in the
  ## span and alpha_2 is zero. Near-breakdown as the exact case rather than as a
  ## perturbation of one.
  f <- lsq_prescribed(12, 4, c(1, 2, 3, 4), seed = 1)
  v1 <- f$V[, 1L, drop = FALSE]
  fit <- lsmr(linop(f$A), f$A %*% v1, tol = 1e-12)
  expect_equal(fit$iterations, 1L)
  expect_true(fit$converged)
  expect_equal(fit$x, v1, tolerance = 1e-12)
})

## --------------------------------------------------- the published recurrence

## Fong and Saunders 2011, single column, no damping, no stopping test and no
## preconditioner: both factorisations and the residual estimate exactly as
## published, so agreement is a statement about the recurrence rather than about
## the bookkeeping around it.
reference_lsmr <- function(M, b, steps) {
  n <- ncol(M)
  Mh <- Conj(t(M))
  x <- rep(0, n)
  if (is.complex(M) || is.complex(b)) storage.mode(x) <- "complex"

  beta <- sqrt(sum(Mod(b)^2)); u <- b / beta
  v <- as.vector(Mh %*% u); alpha <- sqrt(sum(Mod(v)^2)); v <- v / alpha

  zetabar <- alpha * beta
  alphabar <- alpha
  rho <- 1; rhobar <- 1; cbar <- 1; sbar <- 0; zeta <- 0
  h <- v; hbar <- rep(0, n)
  if (is.complex(x)) { storage.mode(h) <- "complex"; storage.mode(hbar) <- "complex" }

  betadd <- beta; betad <- 0; rhodold <- 1; tautildeold <- 0; thetatilde <- 0
  normr <- numeric(steps); normar <- numeric(steps)

  for (k in seq_len(steps)) {
    ut <- as.vector(M %*% v) - alpha * u
    beta <- sqrt(sum(Mod(ut)^2))
    u <- if (beta > 0) ut / beta else ut
    vt <- as.vector(Mh %*% u) - beta * v
    alpha_new <- sqrt(sum(Mod(vt)^2))

    rhoold <- rho
    rho <- sqrt(alphabar^2 + beta^2)
    cs <- alphabar / rho
    sn <- beta / rho
    thetanew <- sn * alpha_new
    alphabar <- cs * alpha_new

    rhobarold <- rhobar
    zetaold <- zeta
    thetabar <- sbar * rho
    rhotemp <- cbar * rho
    rhobar <- sqrt(rhotemp^2 + thetanew^2)
    cbar <- rhotemp / rhobar
    sbar <- thetanew / rhobar
    zeta <- cbar * zetabar
    zetabar <- -sbar * zetabar

    hbar <- h - (thetabar * rho / (rhoold * rhobarold)) * hbar
    x <- x + (zeta / (rho * rhobar)) * hbar
    v <- if (alpha_new > 0) vt / alpha_new else vt
    h <- v - (thetanew / rho) * h
    alpha <- alpha_new

    ## the third rotation sequence, section 3.2, for ||r_k||
    betahat <- cs * betadd
    betadd <- -sn * betadd
    thetatildeold <- thetatilde
    rhotildeold <- sqrt(rhodold^2 + thetabar^2)
    ctildeold <- rhodold / rhotildeold
    stildeold <- thetabar / rhotildeold
    thetatilde <- stildeold * rhobar
    rhodold <- ctildeold * rhobar
    betad <- -stildeold * betad + ctildeold * betahat
    tautildeold <- (zetaold - thetatildeold * tautildeold) / rhotildeold
    taud <- (zeta - thetatilde * tautildeold) / rhodold

    normr[k] <- sqrt((betad - taud)^2 + betadd^2)
    normar[k] <- abs(zetabar)
  }
  list(x = matrix(x, n, 1L), normr = normr, normar = normar)
}

test_that("the iteration is the published one, step for step", {
  ## tol = 0 removes every stopping test, so both run the same fixed number of
  ## steps and any gap is the recurrence itself.
  ##
  ## Four steps, for the reason the LSQR suite records: two arithmetically
  ## equivalent implementations of a Golub-Kahan method do not stay close. The
  ## measured gap here grows from 1e-15 at step 4 to 3e-11 at step 8 and past
  ## 1e-2 by step 12 of a 15-column problem, then closes again as both converge
  ## on the same least-squares solution. SciPy's own LSMR diverges from this one
  ## on the same schedule; the numbers are in
  ## dev_notes/lsmr-and-the-monotone-backward-error.md.
  for (seed in 1:5) {
    f <- lsq_prescribed(45, 15, exp(seq(log(0.5), log(20), length.out = 15)),
                        seed = seed)
    set.seed(600L + seed)
    b <- matrix(stats::rnorm(45), 45, 1)
    mine <- lsmr(linop(f$A), b, tol = 0, maxit = 4L)
    ref <- reference_lsmr(f$A, b, 4L)$x
    expect_equal(mine$iterations, 4L, info = sprintf("seed %d", seed))
    expect_lt(max(Mod(mine$x - ref)) / max(Mod(ref)), 1e-13,
              label = sprintf("seed %d", seed))
  }
})

test_that("the recurrence agrees with the published one on ill-conditioned problems", {
  ## Gate 2 asks for this method specifically. The tolerance tracks the
  ## conditioning rather than being one constant, because the gap between two
  ## implementations of the same recurrence grows with it: measured at four
  ## steps over five seeds it is 3.8e-14 at kappa 1e4, 6.5e-13 at 1e8 and
  ## 1.9e-11 at 1e10.
  for (case in list(list(kappa = 1e4, tol = 1e-12),
                    list(kappa = 1e8, tol = 1e-11),
                    list(kappa = 1e10, tol = 1e-9))) {
    n <- 20
    sigma <- exp(seq(log(1), log(1 / case$kappa), length.out = n))
    for (seed in 1:5) {
      f <- lsq_prescribed(60, n, sigma, seed = seed)
      set.seed(600L + seed)
      b <- matrix(stats::rnorm(60), 60, 1)
      mine <- lsmr(linop(f$A), b, tol = 0, maxit = 4L)$x
      ref <- reference_lsmr(f$A, b, 4L)$x
      expect_lt(max(Mod(mine - ref)) / max(Mod(ref)), case$tol,
                label = sprintf("kappa %.0e, seed %d", case$kappa, seed))
    }
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
    mine <- lsmr(linop(f$A), b, tol = 1e-12, maxit = 500L)
    ## 25 steps on a 15-column operator. In exact arithmetic the Krylov space is
    ## complete at step 15 and the iterate there is the answer; in floating point
    ## the bidiagonal vectors have lost orthogonality by then and the published
    ## recurrence is still wrong in the second digit at step 15, reaching 1e-6
    ## around step 18 and 1e-14 around step 25. Finite termination is the
    ## property this method loses first, which is why the iteration budget
    ## defaults to a multiple of n rather than to n.
    ref <- reference_lsmr(f$A, b, 25L)$x
    expect_lt(max(Mod(mine$x - truth)) / max(Mod(truth)), 1e-9,
              label = sprintf("seed %d, this implementation", seed))
    expect_lt(max(Mod(ref - truth)) / max(Mod(truth)), 1e-9,
              label = sprintf("seed %d, the published one", seed))
  }
})

test_that("columns run in lockstep are the iterates of per-column lsmr", {
  f <- lsq_prescribed(50, 16, exp(seq(log(0.3), log(30), length.out = 16)), seed = 9)
  A <- linop(f$A)
  set.seed(21)
  B <- matrix(stats::rnorm(50 * 4), 50, 4)

  block <- lsmr(A, B, tol = 1e-11)
  for (j in seq_len(4)) {
    one <- lsmr(A, B[, j, drop = FALSE], tol = 1e-11)
    expect_identical(block$x[, j], one$x[, 1L],
                     info = sprintf("column %d diverged from its own solve", j))
  }
})

test_that("a block whose columns converge at different rates still runs in lockstep", {
  f <- lsq_prescribed(40, 12, seq(1, 9, length.out = 12), seed = 22)
  A <- linop(f$A)
  set.seed(23)
  B <- cbind(f$A %*% matrix(stats::rnorm(12), 12, 1),
             matrix(stats::rnorm(40), 40, 1))

  block <- lsmr(A, B, tol = 1e-11)
  for (j in 1:2) {
    one <- lsmr(A, B[, j, drop = FALSE], tol = 1e-11)
    expect_identical(block$x[, j], one$x[, 1L], info = sprintf("column %d", j))
  }
})

## ------------------------------------------- what distinguishes it from lsqr

test_that("the quantity the certificate reports falls monotonically, and lsqr's does not", {
  ## The reason this method is here rather than being a variant of the last one.
  ## LSMR minimises ||A^H r|| over the Krylov space, which is the numerator of
  ## the least-squares backward error the certificate reports, so that number
  ## cannot rise from one step to the next. LSQR minimises ||r|| instead, and
  ## its ||A^H r|| does rise: over 20 steps it rose on every one of 12 seeds at
  ## three conditionings, by up to a factor of 60.
  ##
  ## Both sequences are measured densely from the iterates rather than read off
  ## either recurrence, because neither recurrence's own estimate is a
  ## reportable quantity.
  for (kappa in c(1e2, 1e6)) {
    n <- 16
    sigma <- exp(seq(log(1), log(1 / kappa), length.out = n))
    lsqr_rose <- logical(6)
    for (seed in 1:6) {
      f <- lsq_prescribed(50, n, sigma, seed = seed)
      M <- f$A
      A <- linop(M)
      set.seed(700L + seed)
      b <- matrix(stats::rnorm(50), 50, 1)

      trail <- function(solver) {
        vapply(seq_len(20), function(k) {
          x <- solver(A, b, tol = 0, maxit = k)$x
          sqrt(sum(Mod(crossprod(Conj(M), b - M %*% x))^2))
        }, numeric(1))
      }
      a <- trail(lsmr)
      q <- trail(lsqr)
      ## non-increasing, with a rounding allowance it never needed: the measured
      ## worst rise over all 36 seed and conditioning pairs was exactly zero
      expect_true(all(diff(a) <= 1e-12 * utils::head(a, -1)),
                  info = sprintf("lsmr, kappa %.0e, seed %d", kappa, seed))
      lsqr_rose[seed] <- any(diff(q) > 0)
    }
    expect_true(all(lsqr_rose), info = sprintf("kappa %.0e", kappa))
  }
})

test_that("at a shared budget lsmr reports a smaller backward error than lsqr", {
  ## Same operator, same right-hand side, same budget, same tolerance, and the
  ## comparison is on the certificate's own number rather than on either
  ## recurrence's estimate of it. Measured over 12 seeds the ratio was never
  ## below 8.9 at this budget, so the assertion is that it is smaller rather
  ## than that it is smaller by any particular factor.
  n <- 20
  sigma <- exp(seq(log(1), log(1e-8), length.out = n))
  smaller <- logical(12)
  differs <- logical(12)
  for (seed in 1:12) {
    f <- lsq_prescribed(60, n, sigma, seed = seed)
    A <- linop(f$A)
    set.seed(1200L + seed)
    b <- matrix(stats::rnorm(60), 60, 1)

    a <- lsmr(A, b, tol = 0, maxit = 40L)
    q <- lsqr(A, b, tol = 0, maxit = 40L)
    expect_equal(a$iterations, q$iterations, info = sprintf("seed %d", seed))
    smaller[seed] <- a$certificate$values$backward_error <
                     q$certificate$values$backward_error
    differs[seed] <- !isTRUE(all.equal(a$x, q$x))
  }
  expect_true(all(differs))
  expect_true(all(smaller))
})

## --------------------------------------------------------------- certificate

test_that("an incompatible system does not fail the residual line", {
  m <- 50; n <- 12
  f <- lsq_prescribed(m, n, seq(1, 5, length.out = n), seed = 24)
  set.seed(25)
  b <- matrix(stats::rnorm(m), m, 1)

  cert <- lsmr(linop(f$A), b, tol = 1e-12)$certificate
  expect_equal(cert_status(cert, "residual"), "not_checked")
  expect_equal(cert_status(cert, "backward error"), "pass")
  expect_equal(cert_status(cert, "convergence"), "pass")
  expect_false(cert$overall == "fail")
  detail <- cert$checks$detail[cert$checks$check == "residual"]
  expect_match(detail, "not in the range of A")
})

test_that("the reported backward error is one an exhibited perturbation achieves", {
  ## Stewart's dA = -(r r^H A) / ||r||^2 makes x the exact least-squares solution
  ## of A + dA, and the certificate does not care which method produced x, so
  ## the same check has to hold here as for lsqr.
  m <- 45; n <- 14
  f <- lsq_prescribed(m, n, seq(1, 6, length.out = n), seed = 26)
  M <- f$A
  set.seed(27)
  b <- matrix(stats::rnorm(m), m, 1)

  fit <- lsmr(linop(M), b, tol = 1e-10)
  x <- fit$x
  r <- b - M %*% x
  dA <- -(r %*% (crossprod(Conj(r), M))) / sum(Mod(r)^2)
  Ap <- M + dA

  expect_lt(max(Mod(crossprod(Conj(Ap), b - Ap %*% x))),
            1e-8 * max(Mod(crossprod(Conj(M), b))))
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

    fit <- lsmr(A, b, tol = 1e-11, maxit = 2000L)
    r <- b - M %*% fit$x
    atr <- crossprod(Conj(M), r)
    sv <- svd(M, nu = 0L, nv = 0L)$d

    recovered[i] <- max(Mod(fit$x - lsq_prescribed_solve(f, b))) <=
      1e-6 * max(Mod(lsq_prescribed_solve(f, b)))

    optimality_agrees[i] <-
      abs(fit$certificate$values$optimality -
          sqrt(sum(Mod(atr)^2)) / (fit$certificate$values$norm * sqrt(sum(Mod(r)^2)))) <= 1e-12

    norm_is_a_lower_bound[i] <- fit$certificate$values$norm <= sv[1L] * (1 + 1e-9)

    true_bw <- sqrt(sum(Mod(atr)^2)) / (sv[1L] * sqrt(sum(Mod(r)^2)))
    bounds_true_backward_error[i] <-
      fit$certificate$values$backward_error >= true_bw * (1 - 1e-9)

    est <- lsmr(A, b, tol = 1e-11, maxit = 2000L, norm_control = list(exact_max = 0))
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
  ## The tolerance is set from what the solve delivered rather than guessed,
  ## because the window the floor covers is one c eps wide and a fixed tolerance
  ## either sits inside it or does not depending on the fixture.
  f <- lsq_prescribed(40, 12, seq(1, 4, length.out = 12), seed = 28)
  A <- linop(f$A)
  set.seed(29)
  b <- matrix(stats::rnorm(40), 40, 1)

  achieved <- lsmr(A, b, tol = 0, maxit = 200L)$certificate$values$backward_error
  expect_lt(achieved, 2 * linop:::SOLVE_FLOOR_CONST * .Machine$double.eps)

  with_floor <- lsmr(A, b, tol = achieved / 2, maxit = 200L)
  without <- lsmr(A, b, tol = achieved / 2, maxit = 200L, floor_const = 0)

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

  fit <- lsmr(A, b, tol = 1e-14, maxit = 3L)
  expect_false(fit$converged)
  expect_equal(cert_status(fit$certificate, "convergence"), "fail")
  expect_match(fit$certificate$checks$detail[fit$certificate$checks$check == "convergence"],
               "budget exhausted")
})

## -------------------------------------------------------------- the operator

test_that("lsmr needs an adjoint, as every bidiagonal method does", {
  A <- linop(function(X) rbind(X, 0), dim = c(11L, 10L))
  expect_error(lsmr(A, rep(1, 11)), "adjoint")
})

test_that("a non-conformable right-hand side is refused", {
  f <- lsq_prescribed(20, 6, seq(1, 3, length.out = 6), seed = 32)
  A <- linop(f$A)
  expect_error(lsmr(A, rep(1, 19)), "non-conformable")
  expect_error(lsmr(A, rep(1, 20), x0 = rep(0, 5)), "x0 is")
  expect_error(lsmr(A, rep(1, 20), maxit = 0L), "maxit")
  expect_error(lsmr(matrix(1, 3, 2), rep(1, 3)), "expects a linop")
})

## ---------------------------------------------------------- preconditioning

test_that("left and split preconditioners are refused, and the message says why", {
  f <- lsq_prescribed(20, 6, seq(1, 3, length.out = 6), seed = 34)
  A <- linop(f$A)
  for (sd in c("left", "split")) {
    P <- preconditioner(function(R) R, side = sd, hermitian = TRUE)
    expect_error(lsmr(A, rep(1, 20), preconditioner = P),
                 "acts on the domain of A")
  }
})

test_that("a right preconditioner without M^-H is refused by name", {
  f <- lsq_prescribed(20, 6, seq(1, 3, length.out = 6), seed = 35)
  A <- linop(f$A)
  b <- rep(1, 20)

  bare <- preconditioner(function(R) R, side = "right")
  expect_error(lsmr(A, b, preconditioner = bare), "M\\^-H")

  herm <- preconditioner(function(R) R, side = "right", hermitian = TRUE)
  expect_silent(lsmr(A, b, preconditioner = herm, maxit = 5L))

  supplied <- preconditioner(function(R) R, apply_inverse_adjoint = function(R) R,
                             side = "right")
  expect_silent(lsmr(A, b, preconditioner = supplied, maxit = 5L))
})

test_that("column equilibration converges where the unscaled operator does not", {
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
    plain <- lsmr(A, b, tol = 1e-10, maxit = 40L)
    prec <- lsmr(A, b, tol = 1e-10, preconditioner = P, maxit = 40L)
    helped[s] <- prec$converged && !plain$converged
    differs[s] <- !isTRUE(all.equal(prec$x, plain$x))
  }
  expect_true(all(differs))
  expect_true(all(helped))
})

test_that("a preconditioned solve lands on the same answer as an unpreconditioned one", {
  n <- 14
  d <- 10^seq(-2, 2, length.out = n)
  f <- lsq_prescribed(50, n, seq(1, 5, length.out = n), seed = 37)
  M <- f$A %*% diag(d)
  A <- linop(M)
  P <- preconditioner(function(R) R / d, side = "right", hermitian = TRUE)
  set.seed(38)
  b <- matrix(stats::rnorm(50), 50, 1)

  prec <- lsmr(A, b, tol = 1e-12, preconditioner = P, maxit = 2000L)
  truth <- diag(1 / d) %*% lsq_prescribed_solve(f, b)
  expect_true(prec$converged)
  expect_equal(prec$x, truth, tolerance = 1e-7)
})

## ----------------------------------------------------------- the condition limit

test_that("the condition limit stops a bidiagonal that has stopped meaning anything", {
  n <- 30
  M <- outer(seq(0.1, 2, length.out = 50), seq_len(n) - 1L, "^")
  A <- linop(M)

  for (s in 1:4) {
    set.seed(800L + s)
    b <- matrix(stats::rnorm(50), 50, 1)
    run <- function(cl) lsmr(A, b, tol = 1e-300, maxit = 30L, conlim = cl)
    limited <- run(linop:::KRYLOV_CONDITION_LIMIT)
    unlimited <- run(Inf)
    expect_lte(sqrt(sum(limited$x^2)), sqrt(sum(unlimited$x^2)),
               label = sprintf("seed %d", s))
  }
})

test_that("the condition limit never fires on a well-posed problem", {
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
    a <- lsmr(A, b, tol = 1e-11, maxit = 1200L)
    z <- lsmr(A, b, tol = 1e-11, maxit = 1200L, conlim = Inf)
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
  fit <- lsmr(A, B, tol = 1e-10, history = TRUE)
  expect_equal(ncol(fit$history), 3L)
  expect_equal(nrow(fit$history), fit$iterations)
  expect_null(lsmr(A, B, tol = 1e-10)$history)
})
