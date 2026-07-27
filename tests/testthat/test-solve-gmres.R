## Gate 2 for the third and fourth of the seven methods: recovery against
## closed-form truth run to convergence, the minimal-residual property that
## defines GMRES, the nonsymmetric problems neither CG nor MINRES can be given,
## the three preconditioner sides that this is the first method to tell apart,
## FGMRES on a preconditioner that changes between applications, certificate
## coverage over 20 seeds, and the refusals the section 4.3 table requires.

gmres <- linop:::gmres_solve
fgmres <- linop:::fgmres_solve
cg <- linop:::cg_solve
minres <- linop:::minres_solve
cert_status <- linop:::cert_status

## ------------------------------------------------------------------ recovery

test_that("gmres recovers the solution of a nonsymmetric convection-diffusion operator", {
  n <- 40
  mu <- 0.2
  A <- convdiff_1d(n, mu)
  M <- as.matrix(A)
  ## the fixture has to be genuinely nonsymmetric or the test proves nothing;
  ## the off-diagonals differ by exactly 2 mu
  expect_equal(max(abs(M - t(M))), 2 * mu)

  set.seed(1)
  b <- matrix(stats::rnorm(n), n, 1)
  fit <- gmres(A, b, tol = 1e-12, restart = n, maxit = 500L)
  expect_true(fit$converged)
  expect_equal(fit$x, convdiff_1d_solve(n, mu, b), tolerance = 1e-10)
  expect_equal(fit$method, "gmres")
})

test_that("gmres solves a complex non-hermitian system", {
  ## The case that forces the whole inner product rather than its real part. A
  ## hermitian operator makes the Krylov scalars real, which is why CG and MINRES
  ## can work in real arithmetic; here <v_i, A v_j> is genuinely complex and
  ## keeping only Re() would orthogonalise against the wrong direction.
  n <- 30
  set.seed(2)
  M <- matrix(complex(real = stats::rnorm(n * n), imaginary = stats::rnorm(n * n)),
              n, n) + diag(8, n)
  A <- linop(M)
  ## the dense leaf inspects the entries it holds, so this is a proved FALSE
  ## rather than an unknown, and both cg and minres refuse it
  expect_false(capv(A, "hermitian"))

  x_true <- matrix(complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)), n, 1)
  b <- M %*% x_true
  fit <- gmres(A, b, tol = 1e-12, restart = n, maxit = 500L)
  expect_true(fit$converged)
  expect_true(is.complex(fit$x))
  expect_equal(fit$x, x_true, tolerance = 1e-8)
})

test_that("a vector goes in and a vector comes out", {
  n <- 20
  A <- convdiff_1d(n, 0.2)
  set.seed(4)
  b <- stats::rnorm(n)
  fit <- gmres(A, b, tol = 1e-12, restart = n, maxit = 200L)
  expect_null(dim(fit$x))
  expect_equal(fit$x, as.numeric(convdiff_1d_solve(n, 0.2, b)), tolerance = 1e-10)
})

test_that("a zero right-hand side and a warm start both cost nothing", {
  n <- 20
  A <- convdiff_1d(n, 0.2)

  z <- gmres(A, rep(0, n))
  expect_equal(z$iterations, 0L)
  expect_true(all(z$x == 0))
  expect_true(z$converged)

  x_true <- as.numeric(convdiff_1d_solve(n, 0.2, rep(1, n)))
  w <- gmres(A, rep(1, n), x0 = x_true, tol = 1e-8)
  expect_equal(w$iterations, 0L)
})

test_that("gmres needs no adjoint, which no other method in the package can say", {
  ## Every apply below is mode "N". An operator that supplies only a forward
  ## action has no hermitian or definite capability to declare and cannot be
  ## handed to cg or minres at all, and the norm estimate falls back to probes
  ## rather than failing.
  n <- 30
  M <- as.matrix(convdiff_1d(n, 0.2))
  A <- linop(function(X) M %*% X, dim = c(n, n))
  expect_false(linop:::expr_has_adjoint(A))
  expect_error(linop_apply(A, matrix(1, n, 1), "C"))

  set.seed(5)
  b <- matrix(stats::rnorm(n), n, 1)
  fit <- gmres(A, b, tol = 1e-10, restart = n, maxit = 300L)
  expect_true(fit$converged)
  expect_equal(as.numeric(fit$x), solve(M, as.numeric(b)), tolerance = 1e-7)
})

## ------------------------------ the problems cg and minres cannot be given

test_that("gmres solves what cg and minres both refuse, which is why it exists", {
  n <- 40
  mu <- 0.2
  A <- convdiff_1d(n, mu)
  set.seed(6)
  b <- matrix(stats::rnorm(n), n, 1)

  expect_true(is.na(capv(A, "hermitian")))
  expect_true(is.na(capv(A, "positive_definite")))
  ## both refuse on the capability nobody established, and say which
  expect_error(cg(A, b), "hermitian")
  expect_error(minres(A, b), "hermitian")
  expect_error(minres(A, b), "Unknown is not false")

  fit <- gmres(A, b, tol = 1e-12, restart = n, maxit = 500L)
  expect_true(fit$converged)
  expect_equal(fit$x, convdiff_1d_solve(n, mu, b), tolerance = 1e-10)
})

test_that("a declared hermitian operator that is not one is a gmres problem, not a refusal", {
  ## minres contradicts the declaration and throws, because it needs it. gmres
  ## needs nothing of the operator, so the same object is simply solved.
  set.seed(7)
  n <- 30
  M <- matrix(stats::rnorm(n * n), n, n) + diag(10, n)
  A <- linop(M, properties = c(hermitian = TRUE))
  b <- matrix(stats::rnorm(n), n, 1)

  expect_error(minres(A, b, tol = 1e-10), "declares hermitian = TRUE")

  fit <- gmres(A, b, tol = 1e-12, restart = n, maxit = 300L)
  expect_true(fit$converged)
  expect_equal(as.numeric(fit$x), solve(M, as.numeric(b)), tolerance = 1e-8)
})

## ------------------------------------------------ the minimal-residual property

## The exact minimiser of ||b - A x|| over K_m(A, b), by dense least squares on an
## orthonormalised basis. This is the definition of what GMRES returns at step m,
## computed without reference to the recurrence that produces it.
krylov_argmin_gmres <- function(M, b, m) {
  n <- length(b)
  K <- matrix(0, n, m)
  v <- b / sqrt(sum(b^2))
  for (j in seq_len(m)) {
    K[, j] <- v
    v <- as.numeric(M %*% v)
    for (rep in 1:2) for (i in seq_len(j)) v <- v - sum(K[, i] * v) * K[, i]
    nv <- sqrt(sum(v^2))
    if (nv < 1e-12) { K <- K[, seq_len(j), drop = FALSE]; break }
    v <- v / nv
  }
  as.numeric(K %*% qr.solve(M %*% K, b))
}

test_that("the iterate at step m is the exact minimiser over the Krylov space", {
  ## Not a tolerance statement about the answer but an identity about every
  ## intermediate iterate: the rotations solve the projected least-squares problem
  ## exactly, so each step must agree with a dense solve of that problem.
  n <- 40
  A <- convdiff_1d(n, 0.2)
  M <- as.matrix(A)
  set.seed(8)
  b <- stats::rnorm(n)

  for (m in c(1L, 2L, 3L, 5L, 8L, 12L)) {
    fit <- gmres(A, b, tol = 1e-300, restart = m, maxit = m)
    expect_equal(fit$iterations, m)
    expect_equal(fit$x, krylov_argmin_gmres(M, b, m), tolerance = 1e-9,
                 info = sprintf("step %d is not the Krylov minimiser", m))
  }
})

test_that("the true residual never increases, which is what minimal residual means", {
  n <- 40
  A <- convdiff_1d(n, 0.2)
  M <- as.matrix(A)
  set.seed(9)
  b <- stats::rnorm(n)

  res <- vapply(seq_len(20), function(m) {
    x <- gmres(A, b, tol = 1e-300, restart = m, maxit = m)$x
    sqrt(sum((b - M %*% x)^2))
  }, numeric(1))
  expect_true(all(diff(res) <= 1e-10 * res[1L]))
})

test_that("without a preconditioner the recurrence scalar is the true residual", {
  n <- 40
  A <- convdiff_1d(n, 0.2)
  M <- as.matrix(A)
  set.seed(9)
  b <- stats::rnorm(n)

  fit <- gmres(A, b, tol = 1e-300, restart = 15L, maxit = 15L, history = TRUE)
  truth <- vapply(seq_len(15), function(m) {
    sqrt(sum((b - M %*% gmres(A, b, tol = 1e-300, restart = m, maxit = m)$x)^2))
  }, numeric(1))
  expect_equal(as.numeric(fit$history), truth, tolerance = 1e-9)
})

test_that("gmres terminates in as many steps as the operator has distinct eigenvalues", {
  ## The happy breakdown: the Krylov space becomes invariant and the projected
  ## problem is solved exactly, with h_{j+1,j} at rounding level.
  n <- 30
  A <- linop(spd_prescribed(n, rep(c(1, 2, 4), each = 10), seed = 8))
  set.seed(10)
  b <- matrix(stats::rnorm(n), n, 1)
  fit <- gmres(A, b, tol = 1e-12, restart = n, maxit = 300L)
  expect_lte(fit$iterations, 3L)
  expect_true(fit$converged)
})

test_that("full-dimension gmres finishes in at most n steps", {
  ## K_n is the whole space, so the minimiser over it is the solution. Exactly the
  ## termination property a restarted run gives up in exchange for bounded storage.
  n <- 25
  A <- convdiff_1d(n, 0.2)
  set.seed(11)
  b <- matrix(stats::rnorm(n), n, 1)
  fit <- gmres(A, b, tol = 1e-13, restart = n, maxit = 10L * n)
  expect_lte(fit$iterations, n)
  expect_equal(fit$restarts, 0L)
  expect_true(fit$converged)
})

test_that("columns run in lockstep are the iterates of per-column gmres", {
  n <- 60
  A <- convdiff_1d(n, 0.1)
  set.seed(12)
  B <- matrix(stats::rnorm(n * 4), n, 4)

  for (m in c(10L, 60L)) {
    block <- gmres(A, B, tol = 1e-11, restart = m, maxit = 3000L)
    for (j in seq_len(4)) {
      one <- gmres(A, B[, j, drop = FALSE], tol = 1e-11, restart = m, maxit = 3000L)
      expect_identical(block$x[, j], one$x[, 1L],
                       info = sprintf("restart %d, column %d diverged from its own solve", m, j))
    }
  }
})

## ---------------------------------------------------------------- restarting

test_that("a restarted run reaches the same answer as the unrestarted one", {
  n <- 60
  mu <- 0.1
  A <- convdiff_1d(n, mu)
  set.seed(13)
  b <- matrix(stats::rnorm(n), n, 1)
  truth <- convdiff_1d_solve(n, mu, b)

  full <- gmres(A, b, tol = 1e-12, restart = n, maxit = 3000L)
  expect_equal(full$restarts, 0L)

  for (m in c(5L, 10L, 20L)) {
    fit <- gmres(A, b, tol = 1e-12, restart = m, maxit = 3000L)
    expect_true(fit$converged, info = sprintf("restart %d", m))
    expect_gt(fit$restarts, 0L)
    expect_equal(fit$x, truth, tolerance = 1e-9,
                 info = sprintf("restart %d", m))
    ## bounded storage is bought with iterations, which is the whole trade
    expect_gte(fit$iterations, full$iterations)
  }
})

test_that("restart is clamped to the dimension, since a Krylov space cannot exceed it", {
  n <- 12
  A <- convdiff_1d(n, 0.2)
  set.seed(14)
  b <- matrix(stats::rnorm(n), n, 1)
  big <- gmres(A, b, tol = 1e-12, restart = 500L, maxit = 100L)
  exact <- gmres(A, b, tol = 1e-12, restart = n, maxit = 100L)
  expect_identical(big$x, exact$x)
  expect_equal(big$restarts, 0L)
})

## ------------------------------------------------------ the arithmetic floor

test_that("a fully converged solve does not certify as fail, and the floor is what stops it", {
  n <- 30
  A <- convdiff_1d(n, 0.2)
  set.seed(15)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- gmres(A, b, tol = 0, restart = n, maxit = 200L)
  expect_equal(cert_status(fit$certificate, "residual"), "qualified")
  expect_false(fit$certificate$overall == "fail")

  bare <- gmres(A, b, tol = 0, restart = n, maxit = 200L, floor_const = 0)
  expect_equal(cert_status(bare$certificate, "residual"), "fail")
  expect_equal(bare$certificate$overall, "fail")

  ## the two ran the same iteration; only the certificate differs
  expect_identical(fit$x, bare$x)
})

test_that("a line that meets its tolerance only through the floor says estimate", {
  n <- 30
  A <- convdiff_1d(n, 0.2)
  set.seed(15)
  b <- matrix(stats::rnorm(n), n, 1)

  clean <- gmres(A, b, tol = 1e-8, restart = n)$certificate
  row <- clean$checks[clean$checks$check == "residual", ]
  expect_equal(row$status, "pass")
  expect_equal(row$guarantee, "identity")

  floored <- gmres(A, b, tol = 0, restart = n, maxit = 200L)$certificate
  row <- floored$checks[floored$checks$check == "residual", ]
  expect_equal(row$guarantee, "estimate")
  expect_true(is.na(row$confidence))
})

test_that("forward error is not_checked, and the certificate says why", {
  A <- convdiff_1d(20, 0.2)
  set.seed(16)
  cert <- gmres(A, matrix(stats::rnorm(20), 20, 1), tol = 1e-10)$certificate
  row <- cert$checks[cert$checks$check == "forward error", ]
  expect_equal(row$status, "not_checked")
  expect_match(row$detail, "A\\^-1")
  expect_equal(cert$overall, "qualified")
})

## -------------------------------------------------- certificate coverage, 20

test_that("over 20 seeds the certificate never understates what it reports", {
  n <- 40
  seeds <- seq_len(20)

  recovered <- logical(length(seeds))
  residual_agrees <- logical(length(seeds))
  bounds_true_omega <- logical(length(seeds))
  estimate_is_conservative <- logical(length(seeds))

  for (i in seq_along(seeds)) {
    s <- seeds[i]
    set.seed(2000L + s)
    ## nonsymmetric, and kept away from the origin so the system is well posed
    M <- matrix(stats::rnorm(n * n) / sqrt(n), n, n) + diag(4, n)
    A <- linop(M)
    x_true <- matrix(stats::rnorm(n), n, 1)
    b <- M %*% x_true

    fit <- gmres(A, b, tol = 1e-11, restart = n, maxit = 2000L)
    r <- b - M %*% fit$x
    sv <- svd(M, nu = 0L, nv = 0L)$d

    recovered[i] <- max(Mod(fit$x - x_true)) <= 1e-7 * max(Mod(x_true))
    residual_agrees[i] <- abs(fit$certificate$values$residual -
                              sqrt(sum(r^2)) / sqrt(sum(b^2))) <= 1e-12

    true_omega <- sqrt(sum(r^2)) /
      (sv[1L] * sqrt(sum(fit$x^2)) + sqrt(sum(b^2)))
    bounds_true_omega[i] <-
      fit$certificate$values$backward_error >= true_omega * (1 - 1e-9)

    ## Every route to ||A|| returns a lower bound, so an estimated norm shrinks
    ## the denominator and the reported backward error can only grow.
    est <- gmres(A, b, tol = 1e-11, restart = n, maxit = 2000L,
                 norm_control = list(exact_max = 0))
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

## ------------------------------------------------- the three sides, told apart

## A badly scaled nonsymmetric operator, so a diagonal preconditioner has real
## work to do and the sides can be compared on something they change.
badly_scaled <- function(n, mu, spread) {
  d <- exp(seq(0, log(spread), length.out = n))
  M <- as.matrix(convdiff_1d(n, mu))
  d * M * rep(d, each = n)
}

test_that("every side meets the euclidean tolerance it was asked for", {
  ## The currency question, now with three answers. Right preconditioning
  ## minimises the reported quantity; left minimises ||M^-1 r|| and split
  ## minimises ||r||_{M^-1}, and neither is what the caller asked about. The outer
  ## loop converts the target going in and measures b - A x coming out.
  n <- 40
  M <- badly_scaled(n, 0.2, 1e4)
  A <- linop(M)
  set.seed(17)
  b <- matrix(stats::rnorm(n), n, 1)
  tol <- 1e-11
  d <- abs(diag(M))

  for (sd in c("left", "right", "split")) {
    P <- preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                        positive_definite = TRUE, side = sd, dim = c(n, n))
    fit <- gmres(A, b, tol = tol, preconditioner = P, restart = n, maxit = 3000L)
    true_rel <- sqrt(sum((b - M %*% fit$x)^2)) / sqrt(sum(b^2))
    expect_true(fit$converged, info = sd)
    expect_lte(true_rel, tol)
    expect_equal(max(fit$certificate$values$residual), true_rel, tolerance = 1e-12,
                 info = sd)
  }
})

test_that("the sides are genuinely different iterations, not one implementation relabelled", {
  ## Section 4.3 leaves the gmres row unrestricted, and unlike cg that is not a
  ## statement that the three agree. They generate different Krylov spaces and
  ## different iterates, which is what makes `side` select an algorithm here.
  n <- 40
  M <- badly_scaled(n, 0.2, 1e4)
  A <- linop(M)
  set.seed(18)
  b <- matrix(stats::rnorm(n), n, 1)
  d <- abs(diag(M))

  x <- lapply(c("left", "right", "split"), function(sd) {
    P <- preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                        positive_definite = TRUE, side = sd, dim = c(n, n))
    gmres(A, b, tol = 1e-300, restart = 6L, maxit = 6L, preconditioner = P)$x
  })
  rel <- function(u, v) max(abs(u - v)) / max(abs(u))
  expect_gt(rel(x[[1L]], x[[2L]]), 1e-6)
  expect_gt(rel(x[[1L]], x[[3L]]), 1e-6)
  expect_gt(rel(x[[2L]], x[[3L]]), 1e-6)
})

test_that("a preconditioner earns its keep on a badly scaled operator", {
  ## Both sides get the same budget and the same restart, or the comparison is
  ## between budgets rather than between preconditioners. Both are also required
  ## to certify `pass` rather than `qualified`, so the comparison is between two
  ## converged solves and not between two readings of the arithmetic floor: on a
  ## badly enough scaled operator the floor swallows the tolerance and every run
  ## certifies as met.
  n <- 60
  M <- badly_scaled(n, 0.1, 1e5)
  A <- linop(M)
  set.seed(19)
  b <- matrix(stats::rnorm(n), n, 1)
  tol <- 1e-11
  budget <- 1200L

  plain <- gmres(A, b, tol = tol, restart = n, maxit = budget)
  d <- abs(diag(M))
  P <- preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                      positive_definite = TRUE, side = "right", dim = c(n, n))
  prec <- gmres(A, b, tol = tol, restart = n, maxit = budget, preconditioner = P)

  expect_equal(cert_status(plain$certificate, "residual"), "pass")
  expect_equal(cert_status(prec$certificate, "residual"), "pass")
  expect_lt(prec$iterations, plain$iterations)
})

test_that("a budget too small to converge in is reported and not thrown", {
  n <- 60
  M <- badly_scaled(n, 0.1, 1e5)
  A <- linop(M)
  set.seed(19)
  b <- matrix(stats::rnorm(n), n, 1)

  starved <- gmres(A, b, tol = 1e-11, restart = 20L, maxit = 200L)
  expect_false(starved$converged)
  expect_equal(starved$iterations, 200L)
  expect_equal(cert_status(starved$certificate, "residual"), "fail")
  expect_match(starved$certificate$checks$detail[
    starved$certificate$checks$check == "convergence"], "budget exhausted")
})

test_that("the unpreconditioned path is the right-preconditioned one with M = I", {
  n <- 30
  A <- convdiff_1d(n, 0.2)
  set.seed(20)
  b <- matrix(stats::rnorm(n), n, 1)
  P <- preconditioner(function(R) R, fixed = TRUE, side = "right", dim = c(n, n))
  expect_identical(gmres(A, b, tol = 1e-11, restart = 10L, maxit = 300L)$x,
                   gmres(A, b, tol = 1e-11, restart = 10L, maxit = 300L,
                         preconditioner = P)$x)
})

## ------------------------------------------------------------------- fgmres

test_that("fgmres accepts an inner solve that changes between applications", {
  ## The canonical flexible preconditioner, and the only v0.1 consumer of the
  ## variable_inexact half of the section 4.2 fidelity lattice.
  n <- 40
  M <- badly_scaled(n, 0.2, 1e4)
  A <- linop(M)
  set.seed(21)
  b <- matrix(stats::rnorm(n), n, 1)

  ## a loose inner gmres, which is a different linear map every time it is called
  ## because the Krylov space it builds depends on the block it is given
  inner <- preconditioner(function(R) gmres(A, R, tol = 1e-2, restart = 5L,
                                            maxit = 5L)$x,
                          fixed = FALSE, side = "right", dim = c(n, n))

  fit <- fgmres(A, b, tol = 1e-10, preconditioner = inner, restart = n, maxit = 500L)
  expect_equal(fit$method, "fgmres")
  expect_true(fit$converged)
  true_rel <- sqrt(sum((b - M %*% fit$x)^2)) / sqrt(sum(b^2))
  expect_lte(true_rel, 1e-10)
})

test_that("gmres refuses the flexible preconditioner fgmres accepts, and says which to use", {
  n <- 20
  A <- convdiff_1d(n, 0.2)
  b <- rep(1, n)
  flexible <- preconditioner(function(R) R / 2, fixed = FALSE, side = "right")

  expect_error(gmres(A, b, preconditioner = flexible), "fixed")
  expect_error(gmres(A, b, preconditioner = flexible), "fgmres")
  expect_silent(fgmres(A, b, tol = 1e-10, preconditioner = flexible))
})

test_that("fgmres is defined for right preconditioning only, and says why", {
  n <- 20
  A <- convdiff_1d(n, 0.2)
  b <- rep(1, n)
  for (sd in c("left", "split")) {
    P <- preconditioner(function(R) R / 2, fixed = FALSE, side = sd)
    expect_error(fgmres(A, b, preconditioner = P), "does not accept")
    expect_error(fgmres(A, b, preconditioner = P), "right-preconditioned basis")
  }
})

test_that("with a fixed preconditioner fgmres and gmres agree", {
  ## FGMRES is a flag on GMRES, not a second solver. The two assemble the same
  ## update by different routes, M^-1 (V_m y) against sum_i (M^-1 v_i) y_i, which
  ## are equal in exact arithmetic and not bit for bit, so this is a tolerance
  ## statement and the iteration counts have to match exactly.
  n <- 40
  M <- badly_scaled(n, 0.2, 1e2)
  A <- linop(M)
  set.seed(22)
  b <- matrix(stats::rnorm(n), n, 1)
  d <- abs(diag(M))
  mk <- function(fixed) preconditioner(function(R) R / d, fixed = fixed,
                                       side = "right", dim = c(n, n))

  g <- gmres(A, b, tol = 1e-11, restart = n, maxit = 500L, preconditioner = mk(TRUE))
  f <- fgmres(A, b, tol = 1e-11, restart = n, maxit = 500L, preconditioner = mk(FALSE))
  expect_true(g$converged && f$converged)
  expect_equal(g$iterations, f$iterations)
  expect_equal(g$x, f$x, tolerance = 1e-12)
})

## -------------------------------------------- what split preconditioning needs

test_that("split preconditioning names the property it needs and who claimed it", {
  ## The gmres row asks only for `fixed`, so a split preconditioner arrives either
  ## declaring definiteness or declaring nothing. Those are different facts: one
  ## is a declaration the iteration has contradicted, the other is a requirement
  ## of the method the caller never claimed to meet.
  n <- 20
  A <- convdiff_1d(n, 0.2)
  b <- rep(1, n)

  declared <- preconditioner(function(R) -R, fixed = TRUE, hermitian = TRUE,
                             positive_definite = TRUE, side = "split")
  expect_error(gmres(A, b, preconditioner = declared), "contradicts it")
  expect_error(gmres(A, b, preconditioner = declared), "L L\\^H")

  silent <- preconditioner(function(R) -R, fixed = TRUE, side = "split")
  expect_error(gmres(A, b, preconditioner = silent), "declares no definiteness")

  ## and the same indefinite map is fine on the sides that ask nothing of M
  for (sd in c("left", "right")) {
    P <- preconditioner(function(R) -R, fixed = TRUE, side = sd)
    expect_silent(gmres(A, b, tol = 1e-10, restart = n, preconditioner = P))
  }
})

test_that("a singular left preconditioner is caught rather than divided by", {
  n <- 10
  A <- convdiff_1d(n, 0.2)
  P <- preconditioner(function(R) R * 0, fixed = TRUE, side = "left")
  expect_error(gmres(A, rep(1, n), preconditioner = P), "invertible")
})

test_that("a preconditioner that returns the wrong shape is caught", {
  n <- 10
  A <- convdiff_1d(n, 0.2)
  bad <- preconditioner(function(R) R[-1L, , drop = FALSE], fixed = TRUE)
  expect_error(gmres(A, rep(1, n), preconditioner = bad), "returned a 9 x 1 block")
})

## ------------------------------------------------------------------- refusals

test_that("gmres refuses a rectangular operator and says what a rectangular system is", {
  R <- linop(rmat(5, 4, seed = 1))
  expect_error(gmres(R, rep(1, 4)), "square")
  expect_error(gmres(R, rep(1, 4)), "least-squares")
  expect_error(fgmres(R, rep(1, 4)), "square")
})

test_that("a non-conformable right-hand side is refused", {
  A <- convdiff_1d(6, 0.2)
  expect_error(gmres(A, rep(1, 5)), "non-conformable")
  expect_error(gmres(A, rep(1, 6), x0 = rep(0, 5)), "x0 is")
})

test_that("a nonsensical budget is refused by name", {
  A <- convdiff_1d(6, 0.2)
  expect_error(gmres(A, rep(1, 6), maxit = 0L), "maxit")
  expect_error(gmres(A, rep(1, 6), restart = 0L), "restart")
})

test_that("a non-convergent solve comes back rather than throwing", {
  ## A singular operator with an inconsistent right-hand side has no solution at
  ## all. That is a fact about the problem, so it is reported; a contradicted
  ## declaration would be a fact about the operator, and gmres has none to
  ## contradict.
  n <- 30
  M <- spd_prescribed(n, c(rep(0, 4), seq(1, 5, length.out = 26)), seed = 8)
  ev <- eigen(M, symmetric = TRUE)
  null_dir <- ev$vectors[, which.min(abs(ev$values))]
  A <- linop(M)

  set.seed(23)
  b <- as.numeric(M %*% stats::rnorm(n)) + null_dir
  fit <- gmres(A, b, tol = 1e-12, restart = 10L, maxit = 200L)
  expect_false(fit$converged)
  expect_equal(cert_status(fit$certificate, "residual"), "fail")
  expect_true(is.finite(max(fit$certificate$values$residual)))
})

## ------------------------------------------ condition of the projected problem

## A Vandermonde matrix on 40 nodes, kappa about 3.6e25. Its Krylov basis goes
## numerically rank deficient early, which is the case the condition estimate
## exists for.
vandermonde <- function(n) outer(seq(0.1, 2, length.out = n), seq_len(n) - 1L, "^")

test_that("the condition limit stops a Krylov space that has stopped meaning anything", {
  ## The failure it prevents is not a slow solve but a divergent one: past the
  ## point where the projected problem is numerically singular, the recurrence
  ## residual keeps falling while the true residual climbs and ||x|| runs away.
  ## The certificate reports that honestly either way, so what this buys is the
  ## difference between a bad answer and a much worse one.
  M <- vandermonde(40)
  n <- nrow(M)
  A <- linop(M)

  for (s in 1:4) {
    set.seed(s)
    b <- matrix(stats::rnorm(n), n, 1)
    bn <- sqrt(sum(b^2))
    rel <- function(cl) {
      f <- gmres(A, b, tol = 1e-300, restart = 30L, maxit = 30L, conlim = cl)
      sqrt(sum((b - M %*% f$x)^2)) / bn
    }
    limited <- rel(linop:::GMRES_CONDITION_LIMIT)
    unlimited <- rel(Inf)
    expect_lt(limited * 20, unlimited, label = sprintf("seed %d", s))
    ## and neither certifies as converged, because neither has converged
    expect_gt(limited, 1e-2)
  }
})

test_that("the condition limit never fires on a well-posed problem", {
  ## A stopping rule that ends an iteration has to be silent where it is not
  ## needed, so on operators whose Krylov basis stays healthy the two settings
  ## have to agree bit for bit rather than merely closely.
  fixtures <- list(convdiff = as.matrix(convdiff_1d(60, 0.1)),
                   scaled = badly_scaled(60, 0.1, 1e8),
                   kms = kms_matrix(80, 0.999))
  for (nm in names(fixtures)) {
    M <- fixtures[[nm]]
    n <- nrow(M)
    A <- linop(M)
    set.seed(26)
    b <- matrix(stats::rnorm(n), n, 1)
    a <- gmres(A, b, tol = 1e-11, restart = 20L, maxit = 1200L)
    z <- gmres(A, b, tol = 1e-11, restart = 20L, maxit = 1200L, conlim = Inf)
    expect_identical(a$x, z$x, info = nm)
    expect_equal(a$iterations, z$iterations, info = nm)
  }
})

## ------------------------------------------------------- orthogonalisation

test_that("turning reorthogonalisation off keeps every column deterministic", {
  ## The Daniel-Gragg-Kaufman-Stewart criterion is per column, and the mask that
  ## applies it is exact, so a column that does not need the second pass is left
  ## bitwise untouched by it. That is what keeps a data-dependent test on the
  ## orthogonalisation from costing the lockstep identity, on either setting.
  n <- 40
  A <- convdiff_1d(n, 0.1)
  set.seed(25)
  B <- matrix(stats::rnorm(n * 3), n, 3)
  for (re in c(TRUE, FALSE)) {
    block <- gmres(A, B, tol = 1e-10, restart = 12L, maxit = 1000L, reorth = re)
    for (j in seq_len(3)) {
      one <- gmres(A, B[, j, drop = FALSE], tol = 1e-10, restart = 12L,
                   maxit = 1000L, reorth = re)
      expect_identical(block$x[, j], one$x[, 1L],
                       info = sprintf("reorth = %s, column %d", re, j))
    }
  }
})

test_that("reorthogonalisation changes the iterates it is asked about", {
  ## What is asserted here is only that the setting does something, which is the
  ## precondition for any claim about it. It is deliberately not a claim that
  ## reorthogonalisation improves the answer: over 12 seeds on two ill-conditioned
  ## fixtures the drift of the recurrence scalar from the true residual improves
  ## by a median factor of 1.6 and gets worse on 4 of 12, so the accuracy claim
  ## that would read well here is not one the measurements support. The reasoning
  ## and the numbers are in dev_notes/gmres-and-the-second-pass.md.
  n <- 50
  M <- badly_scaled(n, 0.2, 1e6)
  A <- linop(M)
  set.seed(27)
  b <- matrix(stats::rnorm(n), n, 1)
  on <- gmres(A, b, tol = 1e-300, restart = 40L, maxit = 40L, reorth = TRUE)
  off <- gmres(A, b, tol = 1e-300, restart = 40L, maxit = 40L, reorth = FALSE)
  expect_false(isTRUE(all.equal(on$x, off$x, tolerance = 1e-14)))
})

test_that("turning reorthogonalisation off keeps every column deterministic", {
  ## The selective form, a second pass taken only where the first one cancelled,
  ## would make a column's arithmetic depend on the other columns sharing the
  ## block. This flat one does not, on either setting.
  n <- 40
  A <- convdiff_1d(n, 0.1)
  set.seed(25)
  B <- matrix(stats::rnorm(n * 3), n, 3)
  for (re in c(TRUE, FALSE)) {
    block <- gmres(A, B, tol = 1e-10, restart = 12L, maxit = 1000L, reorth = re)
    for (j in seq_len(3)) {
      one <- gmres(A, B[, j, drop = FALSE], tol = 1e-10, restart = 12L,
                   maxit = 1000L, reorth = re)
      expect_identical(block$x[, j], one$x[, 1L],
                       info = sprintf("reorth = %s, column %d", re, j))
    }
  }
})
