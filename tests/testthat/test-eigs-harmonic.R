harmonic_residual <- function(A, X, values) {
  AX <- linop:::linop_apply(A, X, "N")
  sqrt(colSums(Mod(AX - sweep(X, 2L, values, "*"))^2)) /
    sqrt(colSums(Mod(X)^2))
}

harmonic_expect_pass <- function(fit, k, info = NULL) {
  expect_equal(fit$nconv, k, info = info)
  for (row in c("residual", "backward error", "convergence")) {
    expect_equal(cert_status(fit$certificate, row), "pass", info = info)
  }
  expect_equal(cert_status(fit$certificate, "target identity"), "not_checked", info = info)
}

test_that("harmonic thick and Schur restarts recover interior Laplacian pairs", {
  n <- 80L
  S <- laplacian_1d_eigenvectors(n)
  for (shift in c(0, 2)) for (method in c("lanczos", "arnoldi")) {
    A <- shifted_laplacian_1d(n, shift)
    truth <- laplacian_1d_eigenvalues(n) - shift
    target <- 2.02 - shift
    take <- order(abs(truth - target))[1:3]
    expect_lt(max(abs(linop:::linop_apply(A, S, "N") -
                        sweep(S, 2L, truth, "*"))), 1e-13)
    fit <- eigs(A, 3L, target = target, extraction = "harmonic", method = method,
                ncv = 24L, maxit = 1000L)
    harmonic_expect_pass(fit, 3L, paste(method, shift))
    expect_gt(fit$restarts, 1L)
    expect_lt(fit$iterations, fit$maxit)
    expect_equal(fit$values, truth[take], tolerance = 1e-11)
    expect_lt(max(harmonic_residual(A, fit$vectors, truth[take])), 5e-10)
    for (j in 1:3) {
      expect_gt(abs(sum(Conj(S[, take[j]]) * fit$vectors[, j])), 1 - 1e-10)
    }
    expect_equal(cert_status(fit$certificate, "orthogonality"), "pass")
    expect_equal(cert_status(fit$certificate, "forward error"), "pass")
    expect_true(is.double(fit$vectors))
  }
})

test_that("harmonic Arnoldi recovers convection-diffusion inside its truth budget", {
  n <- 80L
  mu <- 0.04
  rho <- sqrt((1 + mu) / (1 - mu))
  expect_lt(rho^n, 1e4)
  A <- convdiff_1d(n, mu)
  truth <- convdiff_1d_eigenvalues(n, mu)
  S <- rho^seq_len(n) * laplacian_1d_eigenvectors(n)
  S <- sweep(S, 2L, sqrt(colSums(S^2)), "/")
  expect_lt(max(abs(linop:::linop_apply(A, S, "N") -
                      sweep(S, 2L, truth, "*"))), 1e-13)
  take <- order(abs(truth - 2.02))[1:3]
  for (seed in 1:4) {
    fit <- eigs(A, 3L, target = 2.02, extraction = "harmonic", seed = seed,
                ncv = 24L, maxit = 1200L)
    harmonic_expect_pass(fit, 3L, sprintf("seed %d", seed))
    expect_equal(fit$method, "arnoldi harmonic")
    expect_gt(fit$restarts, 0L)
    expect_equal(fit$values, truth[take], tolerance = 1e-8)
    expect_lt(max(harmonic_residual(A, fit$vectors, truth[take])), 2e-9)
    for (j in 1:3) {
      expect_gt(abs(sum(S[, take[j]] * fit$vectors[, j])), 1 - 1e-10)
    }
    expect_equal(cert_status(fit$certificate, "forward error"), "not_checked")
  }
})

test_that("harmonic extraction improves the paired budget across twelve starts", {
  n <- 120L
  for (method in c("lanczos", "arnoldi")) {
    steps <- matrix(0, 12L, 2L, dimnames = list(NULL, c("ritz", "harmonic")))
    for (seed in 1:12) {
      v0 <- linop:::eigen_start_vector(n, "double", NULL, seed)
      for (extraction in colnames(steps)) {
        fit <- eigs(laplacian_1d(n), 3L, target = 2.02, extraction = extraction,
                    method = method, v0 = v0, seed = seed, ncv = 24L,
                    maxit = 900L, tol = 1e-10)
        harmonic_expect_pass(fit, 3L, sprintf("%s %s seed %d", method, extraction, seed))
        steps[seed, extraction] <- fit$iterations
      }
      traces <- lapply(colnames(steps), function(extraction) {
        directions <- list()
        act <- function(X) {
          if (ncol(X) == 1L) directions[[length(directions) + 1L]] <<- X
          laplacian_1d_apply(X)
        }
        A <- linop(act, adjoint = act, dim = c(n, n), properties = c(hermitian = TRUE))
        eigs(A, 3L, target = 2.02, extraction = extraction, method = method,
             v0 = v0, seed = seed, ncv = 24L, maxit = 36L, tol = 1e-10)
        directions
      })
      expect_equal(traces[[1L]][1:24], traces[[2L]][1:24], tolerance = 1e-14)
      expect_lt(abs(sum(Conj(traces[[1L]][[25L]]) * traces[[2L]][[25L]])), 0.999)
    }
    expect_lt(stats::median(steps[, "harmonic"]), stats::median(steps[, "ritz"]))
    expect_lt(mean(steps[, "harmonic"]), mean(steps[, "ritz"]))
    expect_gt(sum(steps[, "harmonic"] < steps[, "ritz"]), nrow(steps) / 2)
  }
})

test_that("the harmonic backward bound covers the true eigenvalue error over twenty starts", {
  n <- 60L
  shift <- 0.37
  A <- shifted_laplacian_1d(n, shift)
  truth <- laplacian_1d_eigenvalues(n) - shift
  for (method in c("lanczos", "arnoldi")) for (seed in 1:20) {
    fit <- eigs(A, 3L, target = 1.65, extraction = "harmonic", method = method,
                seed = seed, ncv = 24L, maxit = 900L)
    info <- sprintf("%s seed %d", method, seed)
    harmonic_expect_pass(fit, 3L, info)
    cv <- fit$certificate$values
    err <- vapply(fit$values, function(z) min(abs(truth - z)), numeric(1))
    expect_true(all(err / cv$norm <= cv$backward_error + cv$floor), info = info)
    expect_true(all(err <= cv$forward_bound + cv$norm * cv$floor), info = info)
  }
})

test_that("target applies A directly while sigma invokes the solver", {
  solves <- 0L
  original <- linop:::solve.linop
  testthat::local_mocked_bindings(solve.linop = function(...) {
    solves <<- solves + 1L
    original(...)
  }, .package = "linop")
  n <- 40L
  applies <- 0L
  act <- function(X) {
    applies <<- applies + ncol(X)
    laplacian_1d_apply(X)
  }
  A <- linop(act, adjoint = act, dim = c(n, n), properties = c(hermitian = TRUE))
  for (method in c("lanczos", "arnoldi")) {
    applies <- 0L
    fit <- eigs(A, 2L, target = 2.02, extraction = "harmonic", method = method,
                ncv = 20L, maxit = 600L)
    harmonic_expect_pass(fit, 2L)
    expect_equal(solves, 0L)
    expect_gte(applies, fit$iterations)
    expect_match(fit$certificate$dispatch$reason, "harmonic extraction near 2.02 on A", fixed = TRUE)
  }
  for (method in c("lanczos", "arnoldi")) {
    shifted <- eigs(A, 2L, sigma = 2.02, method = method, ncv = 20L, maxit = 600L,
                    inner = list(method = "minres"))
    harmonic_expect_pass(shifted, 2L)
  }
  expect_gt(solves, 0L)
  expect_match(shifted$method, "shift-invert", fixed = TRUE)
})

test_that("the reported harmonic value is its measured Rayleigh quotient", {
  A <- laplacian_1d(80)
  for (method in c("lanczos", "arnoldi")) {
    fit <- eigs(A, 2L, target = 2.02, extraction = "harmonic", method = method,
                ncv = 16L, maxit = 16L)
    AX <- laplacian_1d_apply(fit$vectors)
    rq <- colSums(Conj(fit$vectors) * AX) / colSums(Mod(fit$vectors)^2)
    expect_equal(fit$values, Re(rq), tolerance = 1e-14)
    expect_equal(cert_status(fit$certificate, "convergence"), "fail")
    v <- linop:::eigen_start_vector(80, "double", NULL, 1L)
    run <- linop:::lanczos_run(A, v, NULL, 16L, 16L)
    hp <- linop:::harmonic_projection(run$Tm, run$beta_last, 2.02)
    ht <- linop:::harmonic_values(eigen(hp$M)$values, 2.02, hp$reciprocal)
    ht <- ht[order(abs(ht - 2.02))][1:2]
    expect_gt(max(abs(fit$values - ht)), 1e-3)
  }
})

test_that("harmonic projection and recovery preserve the defining equations", {
  A <- laplacian_1d(60)
  v <- linop:::eigen_start_vector(60, "double", NULL, 11L)
  ex <- linop:::lanczos_run(A, v, NULL, 20L, 20L)
  H <- ex$Tm
  tau <- 2.02
  hp <- linop:::harmonic_projection(H, ex$beta_last, tau)
  ee <- eigen(hp$M)
  take <- order(-abs(ee$values))[1:4]
  Z <- ee$vectors[, take]
  delta <- 1 / ee$values[take]
  C <- rbind(H - diag(tau, 20L), c(rep(0, 19L), ex$beta_last))
  expect_lt(max(abs(crossprod(C) %*% Z -
                      sweep(crossprod(H - diag(tau, 20L), Z), 2L, delta, "*"))), 1e-12)
  Q <- linop:::harmonic_keep_basis(hp$M, tau, 8L, FALSE, hp$reciprocal)
  hr <- linop:::harmonic_restart(ex$V, H, ex$beta_last, ex$vnext, Q)
  defect <- laplacian_1d_apply(hr$Y) - hr$Y %*% hr$H - hr$v %*% t(hr$brow)
  expect_lt(max(abs(defect)), 1e-13)
  expect_lt(max(abs(crossprod(hr$Y, hr$v))), 1e-14)
  expect_equal(crossprod(hr$Y), diag(8), tolerance = 1e-13)
  expect_true(is.double(hr$Y) && is.double(hr$H) && is.double(hr$v))
  expect_lt(max(abs(hr$H - t(hr$H))), 1e-13)
  old <- ex$vnext - hr$Y %*% crossprod(hr$Y, ex$vnext)
  expect_gt(sqrt(sum((hr$v - old)^2)), 1e-3)

  H <- matrix(c(0, 1, 1, 0), 2L)
  hp <- linop:::harmonic_projection(H, 0.4, 1)
  expect_true(all(is.finite(hp$M)))
  ee <- eigen(hp$M)
  j <- which.max(abs(ee$values))
  C <- rbind(H - diag(2), c(0, 0.4))
  z <- ee$vectors[, j, drop = FALSE]
  expect_lt(max(abs(crossprod(C) %*% z - crossprod(H - diag(2), z) / ee$values[j])), 1e-13)
})

test_that("an invariant space admits an exact target and a full spectrum request", {
  for (method in c("lanczos", "arnoldi")) for (n in c(1L, 2L, 10L)) {
    fit <- eigs(linop(diag(seq_len(n), n)), min(n, 2L), target = 1,
                extraction = "harmonic", method = method, ncv = n)
    harmonic_expect_pass(fit, min(n, 2L))
    expect_equal(fit$values, as.numeric(seq_len(min(n, 2L))), tolerance = 1e-10)
  }
})

test_that("complex hermitian harmonic extraction recovers the interior", {
  truth <- seq(-2, 3, length.out = 40)
  A <- linop(hpd_prescribed(40, truth, seed = 7L), properties = c(hermitian = TRUE))
  take <- order(abs(truth - 0.37))[1:3]
  for (method in c("lanczos", "arnoldi")) {
    fit <- eigs(A, 3L, target = 0.37, extraction = "harmonic", method = method,
                ncv = 24L, maxit = 1200L)
    harmonic_expect_pass(fit, 3L)
    expect_equal(fit$values, truth[take], tolerance = 1e-10)
    expect_lt(max(harmonic_residual(A, fit$vectors, truth[take])), 1e-10)
    expect_equal(cert_status(fit$certificate, "orthogonality"), "pass")
  }
})

test_that("harmonic extraction converges at an exact interior target", {
  n <- 61L
  truth <- laplacian_1d_eigenvalues(n)
  for (method in c("lanczos", "arnoldi")) {
    fit <- eigs(laplacian_1d(n), 3L, target = 2, extraction = "harmonic",
                method = method, ncv = 24L, maxit = 900L)
    harmonic_expect_pass(fit, 3L)
    expect_equal(fit$values[1L], 2, tolerance = 1e-12)
    take <- vapply(fit$values, function(v) which.min(abs(v - truth)), integer(1))
    expect_setequal(take, 30:32)
    expect_lt(max(harmonic_residual(laplacian_1d(n), fit$vectors, truth[take])), 3e-10)
    expect_equal(cert_status(fit$certificate, "orthogonality"), "pass")
  }
})

test_that("a complex target selects complex harmonic eigenpairs", {
  n <- 40L
  Q <- qr.Q(qr(zmat(n, n, seed = 53L)))
  truth <- seq(-2, 3, length.out = n) + 0.3i * sin(seq_len(n))
  M <- Q %*% (truth * Conj(t(Q)))
  target <- 0.37 + 0.15i
  take <- order(Mod(truth - target))[1:3]
  expect_lt(max(Mod(M %*% Q - sweep(Q, 2L, truth, "*"))), 1e-13)
  fit <- eigs(linop(M), 3L, target = target, extraction = "harmonic",
              ncv = 24L, maxit = 1000L)
  harmonic_expect_pass(fit, 3L)
  expect_equal(fit$values, truth[take], tolerance = 1e-9)
  expect_lt(max(harmonic_residual(linop(M), fit$vectors, truth[take])), 1e-9)
  expect_true(is.complex(fit$vectors))
})

test_that("real harmonic Krylov-Schur keeps conjugate blocks together", {
  p <- 12L
  D <- matrix(0, 2L * p, 2L * p)
  lam <- seq(1, 3, length.out = p) + 1i * seq(0.2, 0.8, length.out = p)
  for (i in seq_len(p)) {
    j <- 2L * i + (-1:0)
    D[j, j] <- matrix(c(Re(lam[i]), -Im(lam[i]), Im(lam[i]), Re(lam[i])), 2L)
  }
  truth <- c(lam, Conj(lam))
  fit <- eigs(linop(D), 3L, target = 1.87, extraction = "harmonic",
              ncv = 16L, maxit = 1000L)
  harmonic_expect_pass(fit, 3L)
  expect_equal(Mod(fit$values - 1.87),
               sort(Mod(truth - 1.87))[1:3], tolerance = 1e-9)
  take <- vapply(fit$values, function(v) which.min(Mod(v - truth)), integer(1))
  expect_lt(max(harmonic_residual(linop(D), fit$vectors, truth[take])), 1e-9)
  expect_gt(fit$restarts, 0L)
})

test_that("harmonic requests refuse ambiguous and unsupported targets by name", {
  A <- laplacian_1d(20)
  expect_error(eigs(A, 2, extraction = "harmonic"), "harmonic.*requires target")
  expect_error(eigs(A, 2, sigma = 2, extraction = "harmonic"), "requires target")
  expect_error(eigs(A, 2, target = 2, sigma = 2), "target and sigma")
  expect_error(eigs(A, 2, target = 2, which = "smallest"), "target and which")
  expect_error(eigs(A, 2, target = NA_real_), "target must be")
  expect_error(eigs(A, 2, target = c(1, 2)), "target must be")
  expect_error(eigs(A, 2, target = 2 + 1i), "target must be real")
  expect_error(eigs(A, 2, target = 2, B = linop_eye(20)), "target with B")
  expect_error(eigs(A, 2, target = 2, inner = list(tol = 1e-6)), "inner.*target")
  expect_error(eigs(linop_jacobi(diagonal = 1), 1, target = 0), "finite dimension")
  expect_error(eigs(finite_section(linop_jacobi(diagonal = 1), 10), 1, target = 0),
               "standard eigenpair certificate")
  for (method in c("primme", "rspectra")) {
    expect_error(eigs(A, 2, target = 2, extraction = "harmonic", method = method),
                 "does not take extraction = 'harmonic'", fixed = TRUE)
  }
})
