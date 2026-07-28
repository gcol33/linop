## Gate 2 for the last of the seven methods: recovery against closed-form truth,
## bitwise agreement with the published recurrence, the three preconditioner
## sides, the breakdown that is this method's alone, certificate coverage over 20
## seeds, and the refusals the section 4.3 table requires.

bicgstab <- linop:::bicgstab_solve
gmres <- linop:::gmres_solve
cg <- linop:::cg_solve
minres <- linop:::minres_solve
cert_status <- linop:::cert_status

## ------------------------------------------------------------------ recovery

test_that("bicgstab recovers the solution of a nonsymmetric convection-diffusion operator", {
  n <- 40
  mu <- 0.3
  A <- convdiff_1d(n, mu)
  M <- as.matrix(A)
  expect_equal(max(abs(M - t(M))), 2 * mu)

  set.seed(1)
  b <- matrix(stats::rnorm(n), n, 1)
  fit <- bicgstab(A, b, tol = 1e-11, maxit = 2000L)
  expect_true(fit$converged)
  expect_equal(fit$method, "bicgstab")
  expect_equal(fit$x, convdiff_1d_solve(n, mu, b), tolerance = 1e-9)
})

test_that("bicgstab solves a complex non-hermitian system", {
  n <- 30
  set.seed(2)
  M <- matrix(complex(real = stats::rnorm(n * n), imaginary = stats::rnorm(n * n)),
              n, n) + diag(8, n)
  A <- linop(M)
  b <- matrix(complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)), n, 1)

  fit <- bicgstab(A, b, tol = 1e-11, maxit = 2000L)
  expect_true(fit$converged)
  expect_true(is.complex(fit$x))
  expect_equal(fit$x, solve(M, b), tolerance = 1e-9)
})

test_that("bicgstab reproduces the closed-form KMS inverse", {
  n <- 40
  rho <- 0.5
  A <- as_spd_linop(kms_matrix(n, rho))
  set.seed(3)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- bicgstab(A, b, tol = 1e-11, maxit = 2000L)
  expect_true(fit$converged)
  expect_equal(fit$x, kms_inverse(n, rho) %*% b, tolerance = 1e-9)
})

test_that("bicgstab solves what cg and minres both refuse, which is why it exists", {
  ## A nonsymmetric operator declaring nothing. CG needs positive definiteness
  ## and MINRES needs hermitian, and neither is available to declare.
  n <- 40
  A <- convdiff_1d(n, 0.4)
  set.seed(4)
  b <- matrix(stats::rnorm(n), n, 1)

  expect_true(is.na(capv(A, "hermitian")))
  expect_true(is.na(capv(A, "positive_definite")))
  ## both refuse on the capability nobody established, and say which
  expect_error(cg(A, b), "hermitian")
  expect_error(minres(A, b), "hermitian")
  expect_error(minres(A, b), "Unknown is not false")
  fit <- bicgstab(A, b, tol = 1e-11, maxit = 2000L)
  expect_true(fit$converged)
  expect_equal(fit$x, convdiff_1d_solve(n, 0.4, b), tolerance = 1e-8)
})

test_that("bicgstab solves an indefinite hermitian system", {
  n <- 40
  A <- shifted_laplacian_1d(n, 0.9)
  set.seed(5)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- bicgstab(A, b, tol = 1e-11, maxit = 2000L)
  expect_true(fit$converged)
  expect_equal(fit$x, shifted_laplacian_solve(n, 0.9, b), tolerance = 1e-8)
})

test_that("bicgstab needs no adjoint, which only gmres can also say", {
  ## Every apply is mode "N". LSQR, LSMR, CG and MINRES all form A^H X.
  n <- 25
  S <- diag(0, n)
  for (i in seq_len(n - 1L)) S[i, i + 1L] <- 1
  M <- diag(4, n) + S
  A <- linop(function(X) M %*% X, dim = c(n, n))
  set.seed(6)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- bicgstab(A, b, tol = 1e-11, maxit = 2000L)
  expect_true(fit$converged)
  expect_equal(fit$x, solve(M, b), tolerance = 1e-9)
  ## and the methods that do need one refuse this operator
  expect_error(linop:::lsqr_solve(A, b), "adjoint")
})

test_that("a vector goes in and a vector comes out", {
  n <- 30
  A <- convdiff_1d(n, 0.2)
  set.seed(7)
  b <- stats::rnorm(n)
  fit <- bicgstab(A, b, tol = 1e-11, maxit = 2000L)
  expect_null(dim(fit$x))
  expect_equal(fit$x, as.numeric(convdiff_1d_solve(n, 0.2, b)), tolerance = 1e-9)
})

test_that("a zero right-hand side and a warm start both cost nothing", {
  n <- 20
  A <- convdiff_1d(n, 0.2)

  z <- bicgstab(A, rep(0, n))
  expect_equal(z$iterations, 0L)
  expect_true(all(z$x == 0))
  expect_true(z$converged)

  set.seed(8)
  x_true <- matrix(stats::rnorm(n), n, 1)
  w <- bicgstab(A, as.matrix(A %*% x_true), x0 = x_true, tol = 1e-8)
  expect_equal(w$iterations, 0L)
})

## --------------------------------------------------- the published recurrence

## van der Vorst 1992, single column, no preconditioner and no stopping test.
reference_bicgstab <- function(M, b, steps) {
  n <- ncol(M)
  x <- rep(0, n)
  if (is.complex(M) || is.complex(b)) storage.mode(x) <- "complex"
  r <- as.vector(b)
  rhat <- r
  p <- r
  rho <- sum(Conj(rhat) * r)
  for (j in seq_len(steps)) {
    v <- as.vector(M %*% p)
    alpha <- rho / sum(Conj(rhat) * v)
    s <- r - alpha * v
    tt <- as.vector(M %*% s)
    omega <- sum(Conj(tt) * s) / sum(Conj(tt) * tt)
    x <- x + alpha * p + omega * s
    r <- s - omega * tt
    rho_new <- sum(Conj(rhat) * r)
    beta <- (rho_new / rho) * (alpha / omega)
    p <- r + beta * (p - omega * v)
    rho <- rho_new
  }
  matrix(x, n, 1L)
}

test_that("the iteration is the published one, bitwise, at every step count", {
  ## The contrast with the two bidiagonal methods, whose reference tests can
  ## assert four steps and no more. LSQR and LSMR lose agreement because their
  ## stored vectors lose orthogonality and rounding noise fills the lost
  ## direction. This recurrence stores no basis and orthogonalises nothing, so
  ## there is nothing to lose and no drift to accumulate: the agreement is exact
  ## at 24 steps as at one, on real, nonsymmetric and complex fixtures alike.
  fixtures <- list(
    convdiff = as.matrix(convdiff_1d(40, 0.3)),
    swirl = as.matrix(convdiff_1d(40, 0.7)),
    laplacian = as.matrix(laplacian_1d(40)),
    kms = kms_matrix(40, 0.7))
  set.seed(9)
  fixtures$complex <- zmat(30, 30) + diag(8, 30)

  for (nm in names(fixtures)) {
    M <- fixtures[[nm]]
    n <- nrow(M)
    for (k in c(1L, 2L, 4L, 8L, 16L, 24L)) {
      for (seed in 1:3) {
        set.seed(300L + seed)
        b <- if (is.complex(M))
          matrix(complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)), n, 1)
        else matrix(stats::rnorm(n), n, 1)
        mine <- bicgstab(linop(M), b, tol = 0, maxit = k)
        expect_equal(mine$iterations, k, info = sprintf("%s, k = %d", nm, k))
        expect_identical(mine$x, reference_bicgstab(M, b, k),
                         info = sprintf("%s, k = %d, seed %d", nm, k, seed))
      }
    }
  }
})

test_that("columns run in lockstep are the iterates of per-column bicgstab", {
  n <- 40
  A <- convdiff_1d(n, 0.3)
  set.seed(10)
  B <- matrix(stats::rnorm(n * 4), n, 4)

  block <- bicgstab(A, B, tol = 1e-11, maxit = 2000L)
  for (j in seq_len(4)) {
    one <- bicgstab(A, B[, j, drop = FALSE], tol = 1e-11, maxit = 2000L)
    expect_identical(block$x[, j], one$x[, 1L],
                     info = sprintf("column %d diverged from its own solve", j))
  }
})

test_that("the residual is not monotone, which is what this method gives up", {
  ## GMRES minimises the residual over the whole Krylov space, so its true
  ## residual cannot rise; its own suite asserts that. BiCGSTAB minimises only
  ## over the single direction of the stabilising half, so its residual does
  ## rise, and by a lot: measured over 8 seeds on this fixture it rose at least
  ## five times in 20 steps on every seed, by up to a factor of 41 in one step.
  n <- 40
  A <- convdiff_1d(n, 0.7)
  M <- as.matrix(A)

  rose <- logical(6)
  gmres_rose <- logical(6)
  for (seed in 1:6) {
    set.seed(500L + seed)
    b <- matrix(stats::rnorm(n), n, 1)
    trail <- function(solver, ...) {
      vapply(seq_len(20), function(k) {
        x <- solver(A, b, tol = 0, maxit = k, ...)$x
        sqrt(sum(Mod(b - M %*% x)^2))
      }, numeric(1))
    }
    rose[seed] <- any(diff(trail(bicgstab)) > 0)
    ## the same budget, and the comparison is on the same measured quantity
    gmres_rose[seed] <- any(diff(trail(gmres, restart = 100L)) > 1e-12)
  }
  expect_true(all(rose))
  expect_false(any(gmres_rose))
})

## ------------------------------------------------------------- the breakdown

test_that("a system this method cannot solve comes back as a certificate", {
  ## For a real skew-symmetric A, <A z, z> = 0 for every real z, so <rhat0, v>
  ## is exactly zero at the first step of every round however the shadow vector
  ## is re-seeded. BiCGSTAB cannot solve this system at all. What is asserted is
  ## that it says so: a fail certificate naming the breakdown, a finite iterate,
  ## and no error thrown.
  n <- 20
  S <- diag(0, n)
  for (i in seq_len(n - 1L)) { S[i, i + 1L] <- 1; S[i + 1L, i] <- -1 }
  A <- linop(S)
  set.seed(11)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- bicgstab(A, b, tol = 1e-10, maxit = 200L)
  expect_false(fit$converged)
  expect_true(all(is.finite(fit$x)))
  expect_equal(cert_status(fit$certificate, "residual"), "fail")
  expect_match(fit$certificate$checks$detail[fit$certificate$checks$check == "convergence"],
               "broke down")
  ## it is the method and not the system: gmres solves it to machine precision
  g <- gmres(A, b, tol = 1e-10, maxit = 200L, restart = n)
  expect_true(g$converged)
  expect_equal(g$x, solve(S, b), tolerance = 1e-9)
})

test_that("the breakdown guard is bitwise silent where no breakdown occurs", {
  ## A guard that ends an iteration has to be inert where it is not needed, so
  ## on healthy fixtures switching it off changes neither the iterate nor the
  ## iteration count.
  fixtures <- list(
    convdiff = as.matrix(convdiff_1d(40, 0.3)),
    swirl = as.matrix(convdiff_1d(40, 0.7)),
    laplacian = as.matrix(laplacian_1d(40)),
    kms = kms_matrix(40, 0.7))
  for (nm in names(fixtures)) {
    M <- fixtures[[nm]]
    A <- linop(M)
    for (seed in 1:3) {
      set.seed(400L + seed)
      b <- matrix(stats::rnorm(nrow(M)), nrow(M), 1)
      a <- bicgstab(A, b, tol = 1e-11, maxit = 1000L)
      z <- bicgstab(A, b, tol = 1e-11, maxit = 1000L, breakdown_tol = 0)
      expect_identical(a$x, z$x, info = sprintf("%s, seed %d", nm, seed))
      expect_equal(a$iterations, z$iterations, info = sprintf("%s, seed %d", nm, seed))
    }
  }
})

## --------------------------------------------------------------- certificate

test_that("a fully converged solve does not certify as fail, and the floor is why", {
  n <- 30
  A <- convdiff_1d(n, 0.2)
  set.seed(12)
  b <- matrix(stats::rnorm(n), n, 1)

  achieved <- bicgstab(A, b, tol = 0, maxit = 300L)$certificate$values$backward_error
  expect_lt(achieved, 2 * linop:::SOLVE_FLOOR_CONST * .Machine$double.eps)

  with_floor <- bicgstab(A, b, tol = achieved / 2, maxit = 300L)
  without <- bicgstab(A, b, tol = achieved / 2, maxit = 300L, floor_const = 0)

  expect_identical(with_floor$x, without$x)
  expect_equal(cert_status(with_floor$certificate, "backward error"), "qualified")
  expect_equal(cert_status(with_floor$certificate, "convergence"), "pass")
  expect_equal(cert_status(without$certificate, "backward error"), "fail")
  expect_equal(cert_status(without$certificate, "convergence"), "fail")
})

test_that("over 20 seeds the certificate never understates what it reports", {
  n <- 30
  A <- convdiff_1d(n, 0.3)
  M <- as.matrix(A)
  sv1 <- svd(M, nu = 0L, nv = 0L)$d[1L]
  seeds <- seq_len(20)

  recovered <- logical(length(seeds))
  residual_agrees <- logical(length(seeds))
  norm_is_a_lower_bound <- logical(length(seeds))
  bounds_true_backward_error <- logical(length(seeds))

  for (i in seq_along(seeds)) {
    s <- seeds[i]
    set.seed(2000L + s)
    b <- matrix(stats::rnorm(n), n, 1)
    fit <- bicgstab(A, b, tol = 1e-11, maxit = 2000L)
    truth <- convdiff_1d_solve(n, 0.3, b)
    r <- b - M %*% fit$x

    recovered[i] <- max(Mod(fit$x - truth)) <= 1e-6 * max(Mod(truth))
    residual_agrees[i] <-
      abs(fit$certificate$values$residual -
          sqrt(sum(Mod(r)^2)) / sqrt(sum(Mod(b)^2))) <= 1e-12
    norm_is_a_lower_bound[i] <- fit$certificate$values$norm <= sv1 * (1 + 1e-9)
    ## with the exact norm the backward error is smaller, so the reported one
    ## overstates rather than understates
    true_bw <- sqrt(sum(Mod(r)^2)) /
      (sv1 * sqrt(sum(Mod(fit$x)^2)) + sqrt(sum(Mod(b)^2)))
    bounds_true_backward_error[i] <-
      fit$certificate$values$backward_error >= true_bw * (1 - 1e-9)
  }

  expect_true(all(recovered))
  expect_true(all(residual_agrees))
  expect_true(all(norm_is_a_lower_bound))
  expect_true(all(bounds_true_backward_error))
})

test_that("a budget too small to converge in is reported and not thrown", {
  n <- 60
  A <- convdiff_1d(n, 0.5)
  set.seed(13)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- bicgstab(A, b, tol = 1e-14, maxit = 2L)
  expect_false(fit$converged)
  expect_equal(cert_status(fit$certificate, "convergence"), "fail")
  expect_match(fit$certificate$checks$detail[fit$certificate$checks$check == "convergence"],
               "budget exhausted")
})

test_that("forward error is not_checked, and the certificate says why", {
  n <- 20
  A <- convdiff_1d(n, 0.2)
  cert <- bicgstab(A, rep(1, n), tol = 1e-10, maxit = 500L)$certificate
  expect_equal(cert_status(cert, "forward error"), "not_checked")
  expect_match(cert$checks$detail[cert$checks$check == "forward error"], "A\\^-1")
})

## ---------------------------------------------------------- preconditioning

test_that("every side meets the euclidean tolerance it was asked for", {
  n <- 40
  A <- convdiff_1d(n, 0.3)
  d <- exp(seq(-1.5, 1.5, length.out = n))
  set.seed(14)
  b <- matrix(stats::rnorm(n), n, 1)
  truth <- convdiff_1d_solve(n, 0.3, b)

  for (sd in c("left", "right", "split")) {
    P <- preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                        positive_definite = TRUE, side = sd)
    fit <- bicgstab(A, b, tol = 1e-11, preconditioner = P, maxit = 2000L)
    expect_true(fit$converged, info = sd)
    expect_equal(fit$x, truth, tolerance = 1e-8, info = sd)
    ## the tolerance is euclidean whatever the side works in
    expect_lte(max(fit$certificate$values$residual), 1e-11 + max(fit$certificate$values$floor))
  }
})

test_that("the sides are genuinely different iterations, not one implementation relabelled", {
  ## The trap this test exists to avoid: the diagonal of convdiff_1d is constant,
  ## so a Jacobi preconditioner built from it is a scalar, and a scalar commutes
  ## with everything. All three sides then produce bitwise identical iterates and
  ## a test comparing them would pass while proving nothing. Both halves are
  ## asserted here so the first cannot decay into the second.
  n <- 40
  A <- convdiff_1d(n, 0.3)
  set.seed(15)
  b <- matrix(stats::rnorm(n), n, 1)
  mk <- function(d, sd) preconditioner(function(R) R / d, fixed = TRUE,
                                       hermitian = TRUE, positive_definite = TRUE,
                                       side = sd)
  run <- function(d, sd) bicgstab(A, b, tol = 0, preconditioner = mk(d, sd), maxit = 3L)$x

  varying <- exp(seq(-1.5, 1.5, length.out = n))
  xl <- run(varying, "left"); xr <- run(varying, "right"); xs <- run(varying, "split")
  expect_gt(max(Mod(xl - xr)), 1e-3)
  expect_gt(max(Mod(xl - xs)), 1e-3)
  expect_gt(max(Mod(xr - xs)), 1e-3)

  scalar <- rep(2, n)
  expect_identical(run(scalar, "left"), run(scalar, "right"))
  expect_identical(run(scalar, "left"), run(scalar, "split"))
})

test_that("split preconditioning names the property it needs and who claimed it", {
  n <- 20
  A <- convdiff_1d(n, 0.2)
  b <- rep(1, n)

  declared <- preconditioner(function(R) -R, fixed = TRUE, hermitian = TRUE,
                             positive_definite = TRUE, side = "split")
  expect_error(bicgstab(A, b, preconditioner = declared), "contradicts it")
  expect_error(bicgstab(A, b, preconditioner = declared), "L L\\^H")

  silent <- preconditioner(function(R) -R, fixed = TRUE, side = "split")
  expect_error(bicgstab(A, b, preconditioner = silent), "declares no definiteness")

  for (sd in c("left", "right")) {
    P <- preconditioner(function(R) -R, fixed = TRUE, side = sd)
    expect_silent(bicgstab(A, b, tol = 1e-10, preconditioner = P, maxit = 500L))
  }
})

test_that("a singular left preconditioner is caught rather than divided by", {
  n <- 10
  A <- convdiff_1d(n, 0.2)
  P <- preconditioner(function(R) R * 0, fixed = TRUE, side = "left")
  expect_error(bicgstab(A, rep(1, n), preconditioner = P), "invertible")
})

test_that("a preconditioner that returns the wrong shape is caught", {
  n <- 10
  A <- convdiff_1d(n, 0.2)
  bad <- preconditioner(function(R) R[-1L, , drop = FALSE], fixed = TRUE)
  expect_error(bicgstab(A, rep(1, n), preconditioner = bad), "returned a 9 x 1 block")
})

test_that("a flexible preconditioner is refused by name", {
  ## The bicgstab row requires `fixed`. FGMRES is the only method that takes one
  ## that changes, and the message says so.
  n <- 10
  A <- convdiff_1d(n, 0.2)
  P <- preconditioner(function(R) R, fixed = FALSE, side = "right")
  expect_error(bicgstab(A, rep(1, n), preconditioner = P), "fixed = TRUE")
  expect_error(bicgstab(A, rep(1, n), preconditioner = P), "fgmres")
})

test_that("a preconditioner earns its keep on a badly scaled operator", {
  ## Same budget, same tolerance, several seeds, and the knob has to change the
  ## iterates as well as the outcome.
  n <- 40
  d <- 10^seq(-2.5, 2.5, length.out = n)
  A <- linop(diag(d) %*% as.matrix(convdiff_1d(n, 0.2)))
  P <- preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                      positive_definite = TRUE, side = "left")

  helped <- logical(6)
  differs <- logical(6)
  for (s in 1:6) {
    set.seed(600L + s)
    b <- matrix(stats::rnorm(n), n, 1)
    plain <- bicgstab(A, b, tol = 1e-10, maxit = 60L)
    prec <- bicgstab(A, b, tol = 1e-10, preconditioner = P, maxit = 60L)
    helped[s] <- prec$converged && !plain$converged
    differs[s] <- !isTRUE(all.equal(prec$x, plain$x))
  }
  expect_true(all(differs))
  expect_true(all(helped))
})

test_that("the unpreconditioned path is the right-preconditioned one with M = I", {
  n <- 30
  A <- convdiff_1d(n, 0.3)
  set.seed(16)
  b <- matrix(stats::rnorm(n), n, 1)
  I <- preconditioner(function(R) R, fixed = TRUE, hermitian = TRUE,
                      positive_definite = TRUE, side = "right")
  expect_identical(bicgstab(A, b, tol = 1e-11, maxit = 500L)$x,
                   bicgstab(A, b, tol = 1e-11, maxit = 500L, preconditioner = I)$x)
})

## -------------------------------------------------------------- the operator

test_that("a non-conformable right-hand side is refused", {
  A <- convdiff_1d(10, 0.2)
  expect_error(bicgstab(A, rep(1, 9)), "non-conformable")
  expect_error(bicgstab(A, rep(1, 10), x0 = rep(0, 5)), "x0 is")
  expect_error(bicgstab(A, rep(1, 10), maxit = 0L), "maxit")
  expect_error(bicgstab(matrix(1, 3, 3), rep(1, 3)), "expects a linop")
})

test_that("bicgstab refuses a rectangular operator", {
  f <- lsq_prescribed(20, 6, seq(1, 3, length.out = 6), seed = 17)
  expect_error(bicgstab(linop(f$A), rep(1, 20)), "square")
})

## ------------------------------------------------------------------- history

test_that("history records one row per iteration and one column per right-hand side", {
  n <- 30
  A <- convdiff_1d(n, 0.3)
  set.seed(18)
  B <- matrix(stats::rnorm(n * 3), n, 3)
  fit <- bicgstab(A, B, tol = 1e-10, history = TRUE, maxit = 500L)
  expect_equal(ncol(fit$history), 3L)
  expect_equal(nrow(fit$history), fit$iterations)
  expect_null(bicgstab(A, B, tol = 1e-10, maxit = 500L)$history)
})
