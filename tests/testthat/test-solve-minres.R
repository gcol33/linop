## Gate 2 for the second of the seven methods: recovery against closed-form truth
## run to convergence, the minimal-residual property that defines MINRES, the
## indefinite problems CG cannot be given at all, certificate coverage over 20
## seeds, and the refusals the section 4.3 table requires.

minres <- linop:::minres_solve
cg <- linop:::cg_solve
cert_status <- linop:::cert_status

## ------------------------------------------------------------------ recovery

test_that("minres recovers the solution of an indefinite shifted Laplacian", {
  n <- 60
  sigma <- 1.5
  lambda <- laplacian_1d_eigenvalues(n) - sigma
  ## the fixture has to be genuinely indefinite or the test proves nothing
  expect_true(any(lambda < 0) && any(lambda > 0))

  A <- shifted_laplacian_1d(n, sigma)
  set.seed(1)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- minres(A, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$x, shifted_laplacian_solve(n, sigma, b), tolerance = 1e-8)
  expect_equal(fit$method, "minres")
})

test_that("minres reproduces the closed-form KMS inverse, all columns at once", {
  n <- 50
  rho <- 0.7
  A <- as_spd_linop(kms_matrix(n, rho))
  fit <- minres(A, diag(1, n), tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$x, kms_inverse(n, rho), tolerance = 1e-7)
})

test_that("minres solves a complex hermitian indefinite system", {
  n <- 25
  A <- linop(hpd_prescribed(n, c(seq(-6, -1, length.out = 10),
                                 seq(1, 9, length.out = 15)), seed = 4),
             properties = c(hermitian = TRUE))
  set.seed(3)
  x_true <- matrix(complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)), n, 1)
  b <- A %*% x_true

  fit <- minres(A, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_true(is.complex(fit$x))
  expect_equal(fit$x, x_true, tolerance = 1e-8)
})

test_that("a vector goes in and a vector comes out", {
  A <- shifted_laplacian_1d(20, 1.5)
  set.seed(4)
  b <- stats::rnorm(20)
  fit <- minres(A, b, tol = 1e-12)
  expect_null(dim(fit$x))
  expect_equal(fit$x, as.numeric(shifted_laplacian_solve(20, 1.5, b)), tolerance = 1e-8)
})

test_that("a zero right-hand side and a warm start both cost nothing", {
  A <- shifted_laplacian_1d(20, 1.5)

  z <- minres(A, rep(0, 20))
  expect_equal(z$iterations, 0L)
  expect_true(all(z$x == 0))
  expect_true(z$converged)

  x_true <- as.numeric(shifted_laplacian_solve(20, 1.5, rep(1, 20)))
  w <- minres(A, rep(1, 20), x0 = x_true, tol = 1e-8)
  expect_equal(w$iterations, 0L)
})

## ------------------------------------------- the problems cg cannot be given

test_that("minres solves what cg refuses, which is the whole reason it exists", {
  ## An indefinite operator cannot carry positive_definite = TRUE, so cg refuses
  ## it by capability before running. This is the dispatch decision method =
  ## "auto" will make once there is more than one method to choose between.
  n <- 40
  sigma <- 1.5
  A <- shifted_laplacian_1d(n, sigma)
  set.seed(5)
  b <- matrix(stats::rnorm(n), n, 1)

  expect_true(is.na(capv(A, "positive_definite")))
  expect_error(cg(A, b), "positive_definite")

  fit <- minres(A, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$x, shifted_laplacian_solve(n, sigma, b), tolerance = 1e-8)
})

test_that("unknown definiteness is a minres problem rather than a refusal", {
  ## The three-valued rule in section 5.3 doing visible work: nobody has proved
  ## this operator definite, and it happens to be, but minres does not need to
  ## know either way.
  n <- 30
  M <- spd_prescribed(n, seq(1, 5, length.out = n), seed = 2)
  A <- linop(M)
  expect_true(is.na(capv(A, "positive_definite")))
  expect_true(isTRUE(capv(A, "hermitian")))

  set.seed(6)
  b <- matrix(stats::rnorm(n), n, 1)
  fit <- minres(A, b, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(as.numeric(fit$x), solve(M, as.numeric(b)), tolerance = 1e-8)
})

## ------------------------------------------------ the minimal-residual property

cdot <- function(x, y) sum(Conj(x) * y)
vnorm <- function(x) sqrt(Re(cdot(x, x)))

## An orthonormal basis for K_m(A, b), by modified Gram-Schmidt with a second
## pass. The threshold on an exhausted space is relative, so scaling the operator
## does not move where the space is declared to end.
krylov_basis <- function(M, b, m) {
  n <- length(b)
  K <- matrix(if (is.complex(M) || is.complex(b)) complex(1) else 0, n, m)
  v <- b / vnorm(b)
  used <- 0L
  for (j in seq_len(m)) {
    K[, j] <- v
    used <- j
    v <- as.vector(M %*% v)
    n0 <- vnorm(v)
    for (pass in 1:2) for (i in seq_len(j)) v <- v - cdot(K[, i], v) * K[, i]
    nv <- vnorm(v)
    if (n0 <= 0 || nv <= 1e-12 * n0) break
    v <- v / nv
  }
  K[, seq_len(used), drop = FALSE]
}

## The exact minimiser of ||b - A x|| over K_m(A, b), by dense least squares over
## that basis. This is the definition of what MINRES returns at step m, computed
## without reference to the recurrence that produces it. The least-squares step is
## a truncated SVD rather than a QR so that it degrades gracefully where the basis
## has lost rank, and so that one routine serves the real and complex cases.
krylov_argmin <- function(M, b, m) {
  K <- krylov_basis(M, b, m)
  G <- M %*% K
  s <- svd(G)
  keep <- s$d > max(dim(G)) * .Machine$double.eps * max(s$d)
  y <- rep(if (is.complex(G) || is.complex(b)) complex(1) else 0, ncol(G))
  if (any(keep)) {
    y[keep] <- as.vector(s$v[, keep, drop = FALSE] %*%
                         (as.vector(crossprod(Conj(s$u[, keep, drop = FALSE]), b)) /
                          s$d[keep]))
  }
  as.vector(K %*% y)
}

test_that("the iterate at step m is the exact minimiser over the Krylov space", {
  ## Not a tolerance statement about the answer but an identity about every
  ## intermediate iterate: the rotations solve the projected least-squares problem
  ## exactly, so each step must agree with a dense solve of that problem.
  n <- 40
  M <- spd_prescribed(n, c(seq(-5, -1, length.out = 15), seq(1, 8, length.out = 25)),
                      seed = 1)
  A <- linop(M, properties = c(hermitian = TRUE))
  set.seed(1)
  b <- stats::rnorm(n)

  for (m in c(1L, 2L, 3L, 5L, 8L, 12L)) {
    fit <- minres(A, b, tol = 1e-300, maxit = m)
    expect_equal(fit$iterations, m)
    expect_equal(fit$x, krylov_argmin(M, b, m), tolerance = 1e-10,
                 info = sprintf("step %d is not the Krylov minimiser", m))
  }
})

test_that("the true residual never increases, which is what minimal residual means", {
  ## CG's residual is not monotone; MINRES's is, by construction. Measured on the
  ## recomputed b - A x rather than on the recurrence, so it is a statement about
  ## the iterates and not about the scalar the rotations carry.
  n <- 40
  M <- spd_prescribed(n, c(seq(-4, -1, length.out = 18), seq(1, 6, length.out = 22)),
                      seed = 3)
  A <- linop(M, properties = c(hermitian = TRUE))
  set.seed(8)
  b <- stats::rnorm(n)

  res <- vapply(seq_len(20), function(m) {
    x <- minres(A, b, tol = 1e-300, maxit = m)$x
    sqrt(sum((b - M %*% x)^2))
  }, numeric(1))
  expect_true(all(diff(res) <= 1e-10 * res[1L]))
})

test_that("without a preconditioner the recurrence scalar is the true residual", {
  n <- 40
  M <- spd_prescribed(n, c(seq(-4, -1, length.out = 18), seq(1, 6, length.out = 22)),
                      seed = 3)
  A <- linop(M, properties = c(hermitian = TRUE))
  set.seed(8)
  b <- stats::rnorm(n)

  fit <- minres(A, b, tol = 1e-300, maxit = 15L, history = TRUE)
  truth <- vapply(seq_len(15), function(m) {
    sqrt(sum((b - M %*% minres(A, b, tol = 1e-300, maxit = m)$x)^2))
  }, numeric(1))
  expect_equal(as.numeric(fit$history), truth, tolerance = 1e-10)
})

test_that("minres terminates in as many steps as the operator has distinct eigenvalues", {
  n <- 30
  A <- linop(spd_prescribed(n, rep(c(-2, 1, 4), each = 10), seed = 8),
             properties = c(hermitian = TRUE))
  set.seed(9)
  b <- matrix(stats::rnorm(n), n, 1)
  fit <- minres(A, b, tol = 1e-12)
  expect_lte(fit$iterations, 3L)
  expect_true(fit$converged)
})

test_that("columns run in lockstep are the iterates of per-column minres", {
  n <- 40
  A <- linop(spd_prescribed(n, c(seq(-9, -1, length.out = 16),
                                 seq(1, 12, length.out = 24)), seed = 9),
             properties = c(hermitian = TRUE))
  set.seed(11)
  B <- matrix(stats::rnorm(n * 4), n, 4)

  block <- minres(A, B, tol = 1e-11)
  for (j in seq_len(4)) {
    one <- minres(A, B[, j, drop = FALSE], tol = 1e-11)
    expect_identical(block$x[, j], one$x[, 1L],
                     info = sprintf("column %d diverged from its own solve", j))
  }
})

## ------------------------------------ agreement with a reference, Gate 2 --

## The published recurrence, transcribed from the Paige and Saunders form with
## the preconditioner set to the identity and no shift. It shares no code with
## the implementation and no idea with krylov_argmin above: one is a short
## recurrence, the other a variational characterisation of the same iterate.
##
## The transcription is checked against SciPy 1.17.1's minres by
## dev_notes/spikes/minres_export.R, which is what makes it a reference rather
## than a second copy of the same beliefs. Over four conditionings and both
## breakdown fixtures it tracks SciPy to the accuracy this implementation does:
## 6.9e-15 at eight steps at kappa 1e2, 4.3e-12 at kappa 1e10.
reference_minres <- function(M, b, steps) {
  n <- length(b)
  cx <- is.complex(M) || is.complex(b)
  zero <- if (cx) complex(n) else numeric(n)
  x <- zero
  r1 <- as.vector(b)
  y <- r1
  beta1 <- sqrt(Re(cdot(r1, y)))
  if (beta1 <= 0) return(x)
  oldb <- 0; beta <- beta1; dbar <- 0; epsln <- 0
  phibar <- beta1
  cs <- -1; sn <- 0
  w <- zero; w2 <- zero
  r2 <- r1
  for (itn in seq_len(steps)) {
    if (beta <= 0) break
    v <- y / beta
    y <- as.vector(M %*% v)
    if (itn >= 2) y <- y - (beta / oldb) * r1
    ## alpha is real for a hermitian operator; taking the real part states that
    ## rather than repairs it.
    alfa <- Re(cdot(v, y))
    y <- y - (alfa / beta) * r2
    r1 <- r2; r2 <- y
    oldb <- beta
    beta <- sqrt(max(Re(cdot(r2, y)), 0))
    oldeps <- epsln
    delta <- cs * dbar + sn * alfa
    gbar <- sn * dbar - cs * alfa
    epsln <- sn * beta
    dbar <- -cs * beta
    gamma <- max(sqrt(gbar^2 + beta^2), .Machine$double.eps)
    cs <- gbar / gamma
    sn <- beta / gamma
    phi <- cs * phibar
    phibar <- sn * phibar
    w1 <- w2; w2 <- w
    w <- (v - oldeps * w1 - delta * w2) / gamma
    x <- x + phi * w
  }
  x
}

## Half below zero and half above, magnitudes from 1 down to 1/kappa, so kappa_2
## is exactly kappa and no definite method may be given the operator at all.
indef_spectrum <- function(n, kappa) {
  rep(c(-1, 1), length.out = n) * exp(seq(log(1), log(1 / kappa), length.out = n))
}

test_that("the iterate agrees with both references on ill-conditioned problems", {
  ## Gate 2 asks for this method specifically. tol = 0 removes every stopping
  ## test, so both sides run the same fixed number of steps and the gap is the
  ## recurrence rather than the budget.
  ##
  ## The tolerance tracks the conditioning rather than being one constant,
  ## because the gap between two computations of the same iterate grows with it.
  ## Measured over five seeds at eight steps in
  ## dev_notes/spikes/minres_reference.R: 3.7e-15 at kappa 1e2, 2.2e-14 at 1e4,
  ## 3.0e-13 at 1e6, 1.1e-12 at 1e8 and 1.7e-11 at 1e10.
  cases <- list(list(kappa = 1e2,  tol = 1e-13),
                list(kappa = 1e4,  tol = 1e-12),
                list(kappa = 1e6,  tol = 1e-11),
                list(kappa = 1e8,  tol = 1e-10),
                list(kappa = 1e10, tol = 1e-9))
  n <- 40
  for (case in cases) {
    for (seed in 1:5) {
      M <- spd_prescribed(n, indef_spectrum(n, case$kappa), seed = seed)
      A <- linop(M, properties = c(hermitian = TRUE))
      set.seed(700L + seed)
      b <- matrix(stats::rnorm(n), n, 1L)
      for (m in c(1L, 2L, 4L, 6L, 8L)) {
        mine <- as.vector(minres(A, b, tol = 0, maxit = m)$x)
        defn <- krylov_argmin(M, b, m)
        publ <- reference_minres(M, as.vector(b), m)
        expect_lt(max(Mod(mine - defn)) / max(Mod(defn)), case$tol,
                  label = sprintf("kappa %.0e, seed %d, step %d, definition",
                                  case$kappa, seed, m))
        expect_lt(max(Mod(mine - publ)) / max(Mod(publ)), case$tol,
                  label = sprintf("kappa %.0e, seed %d, step %d, published",
                                  case$kappa, seed, m))
      }
    }
  }
})

test_that("the two references part company where neither is truth for the other", {
  ## What the ceiling is a property of. The definition and the published
  ## recurrence share no code, and at kappa 1e6 they agree to 1.4e-15 at four
  ## steps, 2.8e-13 at eight, 2.0e-08 at twelve and not at all at sixteen. That
  ## schedule is the same one this implementation follows against each of them,
  ## so the divergence belongs to the method rather than to any one program, and
  ## a reference test cannot be extended past it by tightening anything.
  n <- 40
  M <- spd_prescribed(n, indef_spectrum(n, 1e6), seed = 5L)
  set.seed(705L)
  b <- matrix(stats::rnorm(n), n, 1L)

  gap <- function(m) {
    a <- krylov_argmin(M, b, m)
    max(Mod(a - reference_minres(M, as.vector(b), m))) / max(Mod(a))
  }
  expect_lt(gap(4L), 1e-13)
  expect_lt(gap(8L), 1e-11)
  expect_gt(gap(16L), 1e-3)
})

test_that("an ill-conditioned solve run to convergence lands where kappa allows", {
  ## The other half of the same fact. What the implementations agree on is the
  ## solution, whatever they did on the way to it, and how close that can be is
  ## set by the conditioning: the measured forward error over these fixtures is
  ## between 0.5 and 22 times kappa * eps, so the bound is written that way
  ## rather than as a constant that happens to pass.
  n <- 40
  for (kappa in c(1e2, 1e4, 1e6, 1e8)) {
    for (seed in 1:3) {
      M <- spd_prescribed(n, indef_spectrum(n, kappa), seed = seed)
      A <- linop(M, properties = c(hermitian = TRUE))
      set.seed(800L + seed)
      x_true <- matrix(stats::rnorm(n), n, 1L)
      b <- M %*% x_true
      fit <- minres(A, b, tol = 1e-13, maxit = 40L * n)
      err <- max(Mod(as.vector(fit$x) - as.vector(x_true))) / max(Mod(x_true))
      expect_lt(err, 200 * kappa * .Machine$double.eps,
                label = sprintf("kappa %.0e, seed %d", kappa, seed))
    }
  }
})

test_that("the complex case agrees with both references, where no external one exists", {
  ## SciPy's minres casts a complex operator to real and solves a different
  ## problem, so this case has only the two references carried here. Measured at
  ## eight steps over three seeds: 3.8e-15 at kappa 1e2 and 5.5e-12 at 1e6.
  for (case in list(list(kappa = 1e2, tol = 1e-13),
                    list(kappa = 1e6, tol = 1e-10))) {
    n <- 30
    for (seed in 1:3) {
      M <- hpd_prescribed(n, indef_spectrum(n, case$kappa), seed = seed)
      A <- linop(M, properties = c(hermitian = TRUE))
      set.seed(950L + seed)
      b <- matrix(complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)), n, 1L)
      for (m in c(1L, 2L, 4L, 8L)) {
        mine <- as.vector(minres(A, b, tol = 0, maxit = m)$x)
        defn <- krylov_argmin(M, b, m)
        publ <- reference_minres(M, as.vector(b), m)
        expect_true(is.complex(mine))
        expect_lt(max(Mod(mine - defn)) / max(Mod(defn)), case$tol,
                  label = sprintf("kappa %.0e, seed %d, step %d, definition",
                                  case$kappa, seed, m))
        expect_lt(max(Mod(mine - publ)) / max(Mod(publ)), case$tol,
                  label = sprintf("kappa %.0e, seed %d, step %d, published",
                                  case$kappa, seed, m))
      }
    }
  }
})

## --------------------------------------------------- breakdown, and past it --

## The off-diagonal sequence of the recurrence MINRES runs, by a dense Lanczos
## with no reorthogonalisation, scaled by ||A v_j|| so the row is scale free.
beta_sequence <- function(M, b, steps) {
  out <- rep(NA_real_, steps)
  v <- as.vector(b) / vnorm(b)
  v_old <- rep(0, length(v))
  beta <- 0
  for (j in seq_len(steps)) {
    p <- as.vector(M %*% v)
    w <- p - Re(cdot(v, p)) * v - beta * v_old
    bn <- vnorm(w)
    out[j] <- bn / vnorm(p)
    if (bn <= 0) break
    v_old <- v; v <- w / bn; beta <- bn
  }
  out
}

## b inside a d-dimensional invariant subspace, so the Krylov space is exhausted
## at step d and beta_{d+1} = 0 in exact arithmetic. The shifted Laplacian's
## eigenvectors are closed form, so the subspace is constructed rather than found.
invariant_rhs <- function(n, d, seed) {
  set.seed(seed)
  matrix(laplacian_1d_eigenvectors(n)[, seq_len(d) * 5L, drop = FALSE] %*%
           stats::rnorm(d), n, 1L)
}

test_that("an exact breakdown returns the exact solution at the dimension of the subspace", {
  n <- 60
  sigma <- 1.5
  A <- shifted_laplacian_1d(n, sigma)
  for (d in c(1L, 2L, 3L)) {
    b <- invariant_rhs(n, d, 900L + d)
    fit <- minres(A, b, tol = 1e-12, maxit = 10L * n)
    expect_equal(fit$iterations, d, info = sprintf("d = %d", d))
    expect_equal(cert_status(fit$certificate, "residual"), "pass")
    expect_lt(max(Mod(as.vector(fit$x) -
                      as.vector(shifted_laplacian_solve(n, sigma, b)))) /
              max(Mod(shifted_laplacian_solve(n, sigma, b))), 1e-12,
              label = sprintf("d = %d", d))
    ## and the definition finds the space exhausted at the same step
    expect_equal(ncol(krylov_basis(as.matrix(A), b, n)), d,
                 info = sprintf("d = %d", d))
  }
})

test_that("the off-diagonal that collapses at a breakdown recovers on the next step", {
  ## What an exact breakdown looks like in floating point is one collapsed entry
  ## rather than the end of the sequence: rounding in A v re-seeds the space with
  ## a component the exact iteration does not have, and the next off-diagonal
  ## comes back at O(1). Measured depth of the collapse, which degrades with d:
  ## 8.3e-16 at d = 1, 7.3e-13 at d = 3, 9.6e-11 at d = 8. Measured size of the
  ## recovery: 0.36, 0.14 and 0.042 on the same three.
  ##
  ## This is why nothing in the solver thresholds beta to detect a breakdown, and
  ## why neither reference stops at d either.
  n <- 60
  M <- as.matrix(shifted_laplacian_1d(n, 1.5))
  for (d in c(1L, 2L, 3L, 5L, 8L)) {
    bs <- beta_sequence(M, invariant_rhs(n, d, 900L + d), d + 2L)
    expect_lt(bs[d], 1e-9, label = sprintf("d = %d, the collapse", d))
    expect_gt(bs[d + 1L], 1e-2, label = sprintf("d = %d, the recovery", d))
  }
})

test_that("iterating past a breakdown improves the iterate rather than damaging it", {
  ## The consequence of the previous test, and the reason the residual is the
  ## right thing to stop on. Finite termination is exact only in exact
  ## arithmetic: at d = 8 the iterate at the breakdown step is wrong by 3.9e-12
  ## and the steps taken past it bring that to 9.7e-15, so a solver that
  ## detected the breakdown and stopped would return the worse answer.
  n <- 60
  sigma <- 1.5
  A <- shifted_laplacian_1d(n, sigma)
  for (d in c(1L, 3L, 5L, 8L)) {
    b <- invariant_rhs(n, d, 900L + d)
    truth <- shifted_laplacian_solve(n, sigma, b)
    errs <- vapply(d + 0:8, function(m) {
      max(Mod(as.vector(minres(A, b, tol = 0, maxit = m)$x) - as.vector(truth))) /
        max(Mod(truth))
    }, numeric(1))
    expect_true(all(errs <= 2 * errs[1L]),
                info = sprintf("d = %d: %s", d, paste(sprintf("%.1e", errs),
                                                      collapse = " ")))
    expect_true(all(errs < 1e-11), info = sprintf("d = %d", d))
    if (d >= 5L) expect_lt(errs[length(errs)], errs[1L], label = sprintf("d = %d", d))
  }
})

test_that("the near-breakdown window is crossed without incident", {
  ## b = v_i + delta v_j sits inside a one-dimensional invariant subspace to
  ## relative accuracy delta, so beta_2 is O(delta). The sweep crosses
  ## sqrt(eps) = 1.5e-08, where the solver's exhaustion guard sits, and crosses
  ## the point where one step already meets the tolerance: at delta >= 1e-10 the
  ## solve takes two steps and at delta <= 1e-12 it takes one. Everything in the
  ## window converges, and the worst forward error over it is 4.4e-12, at the
  ## delta where a single step is just good enough to stop on.
  n <- 60
  sigma <- 1.5
  A <- shifted_laplacian_1d(n, sigma)
  V <- laplacian_1d_eigenvectors(n)
  tol <- 1e-12

  for (delta in c(1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16, 0)) {
    b <- matrix(V[, 7L] + delta * V[, 23L], n, 1L)
    truth <- shifted_laplacian_solve(n, sigma, b)
    fit <- minres(A, b, tol = tol, maxit = 10L * n)
    lab <- sprintf("delta %.0e", delta)
    expect_true(fit$converged, info = lab)
    expect_equal(cert_status(fit$certificate, "residual"), "pass", info = lab)
    expect_lte(max(fit$certificate$values$residual), tol, label = lab)
    expect_lt(fit$iterations, 3L, label = lab)
    expect_lt(max(Mod(as.vector(fit$x) - as.vector(truth))) / max(Mod(truth)),
              1e-10, label = lab)
  }
})

test_that("the near-breakdown iterates agree with both references step for step", {
  n <- 60
  A <- shifted_laplacian_1d(n, 1.5)
  M <- as.matrix(A)
  V <- laplacian_1d_eigenvectors(n)
  for (delta in c(1e-4, 1e-8, 1e-12, 1e-16)) {
    b <- matrix(V[, 7L] + delta * V[, 23L], n, 1L)
    for (m in 1:3) {
      mine <- as.vector(minres(A, b, tol = 0, maxit = m)$x)
      defn <- krylov_argmin(M, b, m)
      publ <- reference_minres(M, as.vector(b), m)
      lab <- sprintf("delta %.0e, step %d", delta, m)
      expect_lt(max(Mod(mine - defn)) / max(Mod(defn)), 1e-10, label = lab)
      expect_lt(max(Mod(mine - publ)) / max(Mod(publ)), 1e-10, label = lab)
    }
  }
})

## ------------------------------------------------------ the arithmetic floor

test_that("a fully converged solve does not certify as fail, and the floor is what stops it", {
  n <- 30
  A <- linop(spd_prescribed(n, c(seq(-4, -1, length.out = 12),
                                 seq(1, 4, length.out = 18)), seed = 6),
             properties = c(hermitian = TRUE))
  set.seed(7)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- minres(A, b, tol = 0, maxit = 200L)
  expect_equal(cert_status(fit$certificate, "residual"), "qualified")
  expect_false(fit$certificate$overall == "fail")

  bare <- minres(A, b, tol = 0, maxit = 200L, floor_const = 0)
  expect_equal(cert_status(bare$certificate, "residual"), "fail")
  expect_equal(bare$certificate$overall, "fail")

  ## the two ran the same iteration; only the certificate differs
  expect_identical(fit$x, bare$x)
})

test_that("a line that meets its tolerance only through the floor says estimate", {
  n <- 30
  A <- linop(spd_prescribed(n, c(seq(-4, -1, length.out = 12),
                                 seq(1, 4, length.out = 18)), seed = 6),
             properties = c(hermitian = TRUE))
  set.seed(7)
  b <- matrix(stats::rnorm(n), n, 1)

  clean <- minres(A, b, tol = 1e-8)$certificate
  row <- clean$checks[clean$checks$check == "residual", ]
  expect_equal(row$status, "pass")
  expect_equal(row$guarantee, "identity")

  floored <- minres(A, b, tol = 0, maxit = 200L)$certificate
  row <- floored$checks[floored$checks$check == "residual", ]
  expect_equal(row$guarantee, "estimate")
  expect_true(is.na(row$confidence))
})

test_that("forward error is not_checked, and the certificate says why", {
  A <- shifted_laplacian_1d(20, 1.5)
  set.seed(8)
  cert <- minres(A, matrix(stats::rnorm(20), 20, 1), tol = 1e-10)$certificate
  row <- cert$checks[cert$checks$check == "forward error", ]
  expect_equal(row$status, "not_checked")
  expect_match(row$detail, "A\\^-1")
  expect_equal(cert$overall, "qualified")
})

## -------------------------------------------------- certificate coverage, 20

test_that("over 20 seeds the certificate never understates what it reports", {
  n <- 40
  lambda <- c(-exp(seq(log(1e-1), log(1e2), length.out = 18)),
              exp(seq(log(1e-1), log(1e2), length.out = 22)))
  seeds <- seq_len(20)

  recovered <- logical(length(seeds))
  residual_agrees <- logical(length(seeds))
  bounds_true_omega <- logical(length(seeds))
  estimate_is_conservative <- logical(length(seeds))

  for (i in seq_along(seeds)) {
    s <- seeds[i]
    M <- spd_prescribed(n, lambda, seed = s)
    A <- linop(M, properties = c(hermitian = TRUE))
    set.seed(1000L + s)
    x_true <- matrix(stats::rnorm(n), n, 1)
    b <- M %*% x_true

    fit <- minres(A, b, tol = 1e-11, maxit = 2000L)
    r <- b - M %*% fit$x
    sv <- svd(M, nu = 0L, nv = 0L)$d

    ## kappa is 1e3 here and the tolerance 1e-11, so kappa * tol is what the
    ## residual entitles the solution to and the bound is that rather than a
    ## number chosen to pass.
    recovered[i] <- max(Mod(fit$x - x_true)) <= 1e-8 * max(Mod(x_true))
    residual_agrees[i] <- abs(fit$certificate$values$residual -
                              sqrt(sum(r^2)) / sqrt(sum(b^2))) <= 1e-12

    true_omega <- sqrt(sum(r^2)) /
      (sv[1L] * sqrt(sum(fit$x^2)) + sqrt(sum(b^2)))
    bounds_true_omega[i] <-
      fit$certificate$values$backward_error >= true_omega * (1 - 1e-9)

    ## Every route to ||A|| returns a lower bound, so an estimated norm shrinks
    ## the denominator and the reported backward error can only grow.
    est <- minres(A, b, tol = 1e-11, maxit = 2000L, norm_control = list(exact_max = 0))
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

## ------------------------------------------------------------- preconditioner

test_that("a preconditioned solve meets the euclidean tolerance it was asked for", {
  ## The currency question. MINRES minimises ||r||_{M^-1} and the caller asked in
  ## ||r||_2; the two differ by the conditioning of M, so a solver that reported
  ## its own recurrence scalar would claim a tolerance it had not reached. The
  ## outer loop converts the target going in and measures b - A x coming out.
  n <- 60
  A <- linop(spd_prescribed(n, c(seq(-9, -1, length.out = 25),
                                 seq(1, 9, length.out = 35)), seed = 4),
             properties = c(hermitian = TRUE))
  M <- as.matrix(A)
  set.seed(5)
  b <- matrix(stats::rnorm(n), n, 1)
  tol <- 1e-12

  for (spread in c(1e3, 1e6)) {
    d <- exp(seq(0, log(spread), length.out = n))
    P <- preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                        positive_definite = TRUE, dim = c(n, n))
    fit <- minres(A, b, tol = tol, preconditioner = P, maxit = 5000L)
    true_rel <- sqrt(sum((b - M %*% fit$x)^2)) / sqrt(sum(b^2))
    expect_true(fit$converged, info = sprintf("spread %g", spread))
    expect_lte(true_rel, tol)
    expect_equal(max(fit$certificate$values$residual), true_rel, tolerance = 1e-12)
  }
})

test_that("a jacobi preconditioner turns an indefinite system minres cannot finish into one it can", {
  ## Badly scaled and indefinite. The preconditioner takes the absolute value of
  ## the diagonal, which the indefinite case forces: a preconditioner for minres
  ## has to stay positive definite where the operator does not, or the M^-1 norm
  ## it minimises is not a norm.
  n <- 60
  s <- 10^seq(-2, 2, length.out = n)
  ev <- eigen(kms_matrix(n, 0.5), symmetric = TRUE)
  lam <- ev$values
  lam[seq(1, n, by = 3)] <- -lam[seq(1, n, by = 3)]
  M <- diag(s) %*% (ev$vectors %*% diag(lam) %*% t(ev$vectors)) %*% diag(s)
  M <- (M + t(M)) / 2

  spectrum <- eigen(M, symmetric = TRUE)$values
  expect_true(any(spectrum < 0) && any(spectrum > 0))
  sv <- svd(M, nu = 0L, nv = 0L)$d
  expect_gt(sv[1L] / sv[n], 1e7)

  A <- linop(M, properties = c(hermitian = TRUE))
  set.seed(21)
  b <- matrix(stats::rnorm(n), n, 1)
  tol <- 1e-12
  budget <- 10L * n

  ## Both sides get the same budget, or the comparison is between budgets rather
  ## than between preconditioners.
  plain <- minres(A, b, tol = tol, maxit = budget)
  expect_false(plain$converged)
  expect_equal(plain$iterations, budget)
  expect_equal(cert_status(plain$certificate, "residual"), "fail")
  expect_match(plain$certificate$checks$detail[
    plain$certificate$checks$check == "convergence"], "budget exhausted")

  d <- abs(diag(M))
  P <- preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                      positive_definite = TRUE, dim = c(n, n))
  prec <- minres(A, b, tol = tol, preconditioner = P, maxit = budget)
  expect_true(prec$converged)
  expect_lt(prec$iterations, budget %/% 4L)
  expect_lte(max(prec$certificate$values$residual), tol)
})

test_that("minres accepts a preconditioner on any of the three sides", {
  ## Section 4.3 leaves the MINRES row unrestricted, pending the sourced pass
  ## dev_notes/fgmres-and-preconditioner-sides.md records as open. What this
  ## asserts is the contract as it currently stands, so narrowing the row is a
  ## deliberate edit that changes a passing test.
  n <- 20
  A <- shifted_laplacian_1d(n, 1.5)
  set.seed(22)
  b <- matrix(stats::rnorm(n), n, 1)
  for (sd in c("left", "right", "split")) {
    P <- preconditioner(function(R) R / 2, fixed = TRUE, hermitian = TRUE,
                        positive_definite = TRUE, side = sd, dim = c(n, n))
    expect_silent(minres(A, b, tol = 1e-10, preconditioner = P))
  }
})

test_that("minres refuses a preconditioner whose declared properties it requires and lacks", {
  n <- 20
  A <- shifted_laplacian_1d(n, 1.5)
  b <- rep(1, n)

  flexible <- preconditioner(function(R) R, fixed = FALSE, hermitian = TRUE,
                             positive_definite = TRUE)
  expect_error(minres(A, b, preconditioner = flexible), "fgmres")

  undeclared <- preconditioner(function(R) R)
  expect_error(minres(A, b, preconditioner = undeclared), "hermitian")

  half <- preconditioner(function(R) R, hermitian = TRUE)
  expect_error(minres(A, b, preconditioner = half), "positive_definite")
})

test_that("a preconditioner that returns the wrong shape is caught", {
  n <- 10
  A <- shifted_laplacian_1d(n, 1.5)
  bad <- preconditioner(function(R) R[-1L, , drop = FALSE], fixed = TRUE,
                        hermitian = TRUE, positive_definite = TRUE)
  expect_error(minres(A, rep(1, n), preconditioner = bad), "returned a 9 x 1 block")
})

## ------------------------------------------------------------------- refusals

test_that("minres names the capability it needed", {
  n <- 10
  A <- linop(function(X) X * 2, adjoint = function(X) X * 2, dim = c(n, n))
  expect_true(is.na(capv(A, "hermitian")))
  expect_error(minres(A, rep(1, n)), "hermitian")
  expect_error(minres(A, rep(1, n)), "Unknown is not false")
})

test_that("minres refuses a rectangular operator and says what a rectangular system is", {
  R <- linop(rmat(5, 4, seed = 1))
  expect_error(minres(R, rep(1, 4)), "square")
  expect_error(minres(R, rep(1, 4)), "least-squares")
})

test_that("minres contradicts a false hermitian declaration by name", {
  set.seed(2)
  n <- 30
  N <- matrix(stats::rnorm(n * n), n, n)
  A <- linop(N, properties = c(hermitian = TRUE))
  expect_error(minres(A, rep(1, n), tol = 1e-10), "declares hermitian = TRUE")
  expect_error(minres(A, rep(1, n), tol = 1e-10), "contradicts it")
})

test_that("minres contradicts a false declaration on the preconditioner by name", {
  A <- shifted_laplacian_1d(4, 1.5)
  P <- preconditioner(function(R) -R, fixed = TRUE, hermitian = TRUE,
                      positive_definite = TRUE)
  expect_error(minres(A, c(1, 0, 0, 0), preconditioner = P), "preconditioner declares")
})

test_that("a non-conformable right-hand side is refused", {
  A <- shifted_laplacian_1d(6, 1.5)
  expect_error(minres(A, rep(1, 5)), "non-conformable")
  expect_error(minres(A, rep(1, 6), x0 = rep(0, 5)), "x0 is")
})

test_that("a non-convergent solve comes back rather than throwing", {
  ## An inconsistent singular system has no solution at all. The distinction the
  ## solver draws is between that, which is a fact about the problem and is
  ## reported, and a contradicted declaration, which is a fact about the operator
  ## and is raised.
  n <- 30
  lambda <- c(rep(0, 5), seq(-3, 4, length.out = 25))
  M <- spd_prescribed(n, lambda, seed = 8)
  ev <- eigen(M, symmetric = TRUE)
  null_dir <- ev$vectors[, which.min(abs(ev$values))]
  A <- linop(M, properties = c(hermitian = TRUE))

  set.seed(10)
  b <- as.numeric(M %*% stats::rnorm(n)) + null_dir

  fit <- minres(A, b, tol = 1e-10, maxit = 200L)
  expect_false(fit$converged)
  expect_equal(cert_status(fit$certificate, "residual"), "fail")
  expect_true(is.finite(max(fit$certificate$values$residual)))
})

## ---------------------------------------------------------- evidence for auto

test_that("the requirement method = auto would apply admits proof and refuses declaration", {
  n <- 8
  req <- linop:::MINRES_HERMITIAN_REQUIREMENT

  ## a bare declaration is the caller's claim, not an argument
  declared <- linop(function(X) X, adjoint = function(X) X, dim = c(n, n),
                    properties = c(hermitian = TRUE))
  expect_false(evidence_satisfies(linop:::cape(declared, "hermitian"), req))

  ## the dense leaf establishes it by inspecting the entries it holds
  proved <- linop(spd_prescribed(n, seq(1, 3, length.out = n), seed = 1))
  expect_true(evidence_satisfies(linop:::cape(proved, "hermitian"), req))

  expect_false(evidence_satisfies(linop:::ev_probe(), req))

  ## naming the method accepts the declaration; the requirement is for auto
  expect_silent(minres(declared, rep(1, n), tol = 1e-10))
})

## ------------------------------------------- what the contradiction test is on

test_that("the contradiction test rests on the adjoint identity, not the recurrence", {
  ## Measured, not assumed. On a clustered spectrum classical Lanczos loses
  ## orthogonality, and with it the identity v_{j-1}^H A v_j = beta_j, for a
  ## perfectly hermitian operator. Anyone tempted to test symmetry that way, as a
  ## quantity the recurrence already carries, gets a solver that refuses correct
  ## operators. <x, A y> = <A x, y> does not depend on orthogonality and stays at
  ## rounding level across the same fixture.
  n <- 50
  M <- spd_prescribed(n, c(rep(1, n - 3), 1e6, -1e6, 1e-3), seed = 1)
  set.seed(101)
  b <- stats::rnorm(n)

  worst <- c(recurrence = 0, adjoint = 0)
  beta1 <- sqrt(sum(b^2)); v <- b / beta1; u <- v
  uo <- rep(0, n); vo <- rep(0, n); po <- rep(0, n); beta <- 0
  for (j in seq_len(250)) {
    p <- as.numeric(M %*% v)
    if (j >= 2) {
      sc <- sqrt(sum(vo^2)) * sqrt(sum(p^2))
      if (sc > 0) {
        worst["recurrence"] <- max(worst["recurrence"], abs(sum(vo * p) - beta) / sc)
        worst["adjoint"] <- max(worst["adjoint"], abs(sum(vo * p) - sum(po * v)) / sc)
      }
    }
    al <- sum(v * p)
    w <- p - al * u - beta * uo
    bn <- sqrt(sum(w^2))
    if (bn <= 1e-14 * sqrt(sum(p^2))) break
    uo <- u; vo <- v; po <- p; u <- w / bn; v <- w / bn; beta <- bn
  }

  tol <- linop:::MINRES_SYMMETRY_TOL
  expect_gt(worst[["recurrence"]], tol)   # would refuse a hermitian operator
  expect_lt(worst[["adjoint"]], tol)      # the identity the solver actually uses

  ## and the solver itself is untroubled by the same fixture
  A <- linop(M, properties = c(hermitian = TRUE))
  expect_silent(minres(A, b, tol = 1e-8, maxit = 500L))
})

test_that("no hermitian fixture in the suite trips the contradiction test", {
  spectra <- list(
    spread    = c(-exp(seq(0, log(1e5), length.out = 25)),
                   exp(seq(0, log(1e5), length.out = 25))),
    clustered = rep(c(-1, 1, -1.0000001, 1.0000001), length.out = 50),
    gapped    = c(rep(1, 47), 1e6, -1e6, 1e-3),
    tiny      = c(1e-10, -1e-9, seq(-3, 5, length.out = 48)),
    repeats   = rep(c(-2, 3), length.out = 50)
  )
  for (nm in names(spectra)) {
    for (s in 1:4) {
      A <- linop(spd_prescribed(50, spectra[[nm]], seed = s),
                 properties = c(hermitian = TRUE))
      set.seed(100L + s)
      expect_error(minres(A, stats::rnorm(50), tol = 1e-10, maxit = 400L), NA,
                   info = sprintf("%s, seed %d", nm, s))
    }
  }
})
