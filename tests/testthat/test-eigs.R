## Section 7.2, eigs(). Recovery against closed-form spectra, certificate
## coverage across seeds, and the refusals.
##
## Nothing here compares sorted spectra. The convdiff fixture in helper-linop.R
## records why: a closed form whose eigenvalue set is right and whose eigenvectors
## are wrong passes a sorted comparison and fails everything that matters. Every
## recovery test below checks A x - lambda x, or the eigenvector itself.

## max ||A x - lambda x|| / (||A|| ||x||), recomputed from the operator rather
## than read off the certificate, against a ||A|| the fixture knows in closed
## form. That is the quantity tol is expressed in: the run stops on
## ||r|| <= tol * anorm_lb with anorm_lb a lower bound on ||A||, so the relative
## backward error is at most tol and comparing an absolute residual against tol
## would be off by a factor of ||A||.
pair_backward <- function(A, values, vectors, anorm) {
  R <- linop_apply(A, vectors, "N") - vectors * rep(values, each = nrow(vectors))
  max(sqrt(colSums(Mod(R)^2)) / sqrt(colSums(Mod(vectors)^2))) / anorm
}

## The distance from each reported value to the nearest true eigenvalue, which is
## what the forward-error row bounds.
nearest_gap <- function(values, truth) {
  vapply(values, function(v) min(abs(v - truth)), numeric(1))
}

test_that("eigs recovers both ends of the 1-D Laplacian, values and vectors", {
  n <- 60
  A <- laplacian_1d(n)
  truth <- laplacian_1d_eigenvalues(n)
  Vtrue <- laplacian_1d_eigenvectors(n)

  ## Smallest algebraic are eigenvalues 1..3 of the closed form, largest are n..n-2.
  lo <- eigs(A, k = 3, which = "smallest_algebraic")
  expect_equal(lo$nconv, 3L)
  expect_equal(lo$values, truth[1:3], tolerance = 1e-12)
  expect_lte(pair_backward(A, lo$values, lo$vectors, max(truth)), 1e-10)

  ## The eigenvector, not only the eigenvalue: each column is the closed-form one
  ## up to a sign, which is the check a sorted spectrum cannot make.
  for (j in 1:3) {
    got <- lo$vectors[, j]
    want <- Vtrue[, j]
    expect_lt(min(max(abs(got - want)), max(abs(got + want))), 1e-8)
  }

  hi <- eigs(A, k = 3, which = "largest_algebraic")
  expect_equal(hi$nconv, 3L)
  expect_equal(sort(hi$values), sort(truth[(n - 2):n]), tolerance = 1e-12)
  expect_lte(pair_backward(A, hi$values, hi$vectors, max(truth)), 1e-10)
})

test_that("largest magnitude is not largest algebraic on an indefinite operator", {
  n <- 40
  sigma <- 2
  A <- shifted_laplacian_1d(n, sigma)
  truth <- laplacian_1d_eigenvalues(n) - sigma
  expect_true(any(truth < 0) && any(truth > 0))

  lm <- eigs(A, k = 2, which = "largest")
  la <- eigs(A, k = 2, which = "largest_algebraic")
  expect_equal(lm$nconv, 2L)
  expect_equal(la$nconv, 2L)
  expect_equal(sort(lm$values), sort(truth[order(-abs(truth))][1:2]), tolerance = 1e-10)
  expect_equal(sort(la$values), sort(truth[order(-truth)][1:2]), tolerance = 1e-10)
  ## and the two requests genuinely differ here, so neither is the other relabelled
  expect_false(isTRUE(all.equal(sort(lm$values), sort(la$values))))
})

test_that("eigs recovers a prescribed spectrum with a dialled-in gap", {
  n <- 40
  lambda <- c(100, 99.5, 40, 10, seq(1, 0.1, length.out = n - 4L))
  M <- spd_prescribed(n, lambda, seed = 7L)
  A <- as_spd_linop(M)

  fit <- eigs(A, k = 3, which = "largest_algebraic")
  expect_equal(fit$nconv, 3L)
  expect_equal(fit$values, lambda[1:3], tolerance = 1e-9)
  expect_lte(pair_backward(A, fit$values, fit$vectors, max(lambda)), 1e-10)
})

test_that("eigs recovers a complex hermitian spectrum and keeps complex vectors", {
  n <- 30
  lambda <- c(50, 20, 8, seq(2, 0.5, length.out = n - 3L))
  M <- hpd_prescribed(n, lambda, seed = 3L)
  A <- linop(M, properties = c(hermitian = TRUE, positive_definite = TRUE))

  fit <- eigs(A, k = 3, which = "largest_algebraic")
  expect_equal(fit$nconv, 3L)
  expect_equal(fit$values, lambda[1:3], tolerance = 1e-9)
  expect_true(is.complex(fit$vectors))
  ## the reported value is real, because a Rayleigh quotient of a hermitian
  ## operator is
  expect_false(is.complex(fit$values))
  expect_lte(pair_backward(A, fit$values, fit$vectors, max(lambda)), 1e-10)
})

test_that("a repeated eigenvalue comes back as an invariant subspace, not one vector twice", {
  n <- 30
  lambda <- c(10, 10, 10, seq(2, 0.5, length.out = n - 3L))
  M <- spd_prescribed(n, lambda, seed = 11L)
  A <- as_spd_linop(M)

  fit <- eigs(A, k = 3, which = "largest_algebraic")
  expect_equal(fit$nconv, 3L)
  expect_equal(fit$values, rep(10, 3), tolerance = 1e-9)
  ## Three distinct directions. Locking deflates rather than excludes, so a
  ## second eigenvector for the same eigenvalue is still reachable.
  expect_lt(max(abs(crossprod(fit$vectors) - diag(3))), 1e-10)
  expect_lte(pair_backward(A, fit$values, fit$vectors, max(lambda)), 1e-10)
})

test_that("shift-invert finds the eigenvalues nearest sigma, which the plain run does not", {
  n <- 60
  A <- laplacian_1d(n)
  truth <- laplacian_1d_eigenvalues(n)
  sigma <- 2
  want <- truth[order(abs(truth - sigma))][1:3]

  fs <- eigs(A, k = 3, sigma = sigma)
  expect_equal(fs$nconv, 3L)
  expect_equal(sort(fs$values), sort(want), tolerance = 1e-9)
  expect_lte(pair_backward(A, fs$values, fs$vectors, max(truth)), 1e-10)

  ## The knob did what is claimed: an unshifted run at the same budget reaches
  ## the ends of the spectrum instead, so these are not the same pairs relabelled.
  fp <- eigs(A, k = 3, which = "largest_algebraic")
  expect_gt(min(abs(outer(fp$values, want, "-"))), 1e-3)
})

test_that("eigs runs on an operator supplying no adjoint", {
  n <- 40
  A <- linop(laplacian_1d_apply, dim = c(n, n),
             properties = c(hermitian = TRUE, positive_definite = TRUE))
  expect_error(linop_apply(A, matrix(1, n, 1), "C"), "no adjoint")

  fit <- eigs(A, k = 2, which = "largest_algebraic")
  expect_equal(fit$nconv, 2L)
  expect_equal(sort(fit$values), sort(laplacian_1d_eigenvalues(n)[(n - 1):n]),
               tolerance = 1e-10)
})

test_that("the thick restart is what makes a small subspace converge", {
  ## This exact configuration, with the restart taken from one vector instead,
  ## spends 240 of its 300 iterations over 9 rounds and converges none of the
  ## four, stalling at a backward error of 5.7e-4. The shipped restart reaches
  ## 5.1e-13 in 96. Measured in dev_notes/spikes/restart-comparison.R, which
  ## swaps only the block that chooses what the next round starts from.
  n <- 60
  A <- laplacian_1d(n)
  fit <- eigs(A, k = 4, which = "smallest_algebraic", ncv = 24L, maxit = 300L)
  expect_equal(fit$nconv, 4L)
  expect_gt(fit$restarts, 0L)
  expect_lt(fit$iterations, 200L)
  expect_equal(fit$values, laplacian_1d_eigenvalues(n)[1:4], tolerance = 1e-10)
})

test_that("the forward-error bound contains the true error over 20 seeds", {
  n <- 30
  lambda <- c(30, 12, 11.8, 5, seq(2, 0.4, length.out = n - 4L))
  ## The certificate's forward line is a deterministic bound, so it holds every
  ## time rather than at a nominal rate; the seeds vary the operator and the
  ## starting vector together.
  worst_slack <- 0
  for (s in 1:20) {
    M <- spd_prescribed(n, lambda, seed = s)
    A <- as_spd_linop(M)
    fit <- eigs(A, k = 3, which = "largest_algebraic", seed = s)
    truth <- eigen(M, symmetric = TRUE, only.values = TRUE)$values

    bound <- fit$certificate$values$forward_bound
    err <- nearest_gap(fit$values, truth)
    nA <- fit$certificate$values$norm
    expect_true(all(err <= bound + 8 * .Machine$double.eps * nA),
                info = sprintf("seed %d: worst err %.3e against bound %.3e", s,
                               max(err), max(bound)))
    worst_slack <- max(worst_slack, max(err - bound))
    expect_equal(cert_status(fit$certificate, "forward error"), "pass")
  }
  ## The bound is not vacuous: it is not orders of magnitude above the error it
  ## bounds on every seed, or it would contain anything.
  expect_lt(worst_slack, 1e-12)
})

test_that("the certificate is qualified and never pass, because target identity is not checked", {
  A <- laplacian_1d(40)
  fit <- eigs(A, k = 2, which = "largest_algebraic")
  expect_equal(fit$certificate$overall, "qualified")
  expect_equal(cert_status(fit$certificate, "target identity"), "not_checked")
  expect_true("target identity" %in% fit$certificate$without_deterministic_bound)
})

test_that("the forward bound counts as a deterministic bound in the summary line", {
  ## Before eigs() no row carried guarantee 'deterministic_bound', so the summary
  ## line's test could read != 'identity' without contradiction. It cannot now.
  M <- spd_prescribed(20, seq(10, 1, length.out = 20), seed = 2L)
  fit <- eigs(as_spd_linop(M), k = 2, which = "largest_algebraic")
  row <- fit$certificate$checks
  expect_equal(row$guarantee[row$check == "forward error"], "deterministic_bound")
  expect_false("forward error" %in% fit$certificate$without_deterministic_bound)
})

test_that("the forward bound records the declaration it rests on", {
  ## A dense symmetric leaf establishes hermitian by an exact check on data it
  ## already holds, so the bound clears a requirement that excludes declarations.
  M <- spd_prescribed(20, seq(10, 1, length.out = 20), seed = 5L)
  computed <- eigs(linop(M), k = 2, which = "largest_algebraic")
  ev_computed <- computed$certificate$evidence[["forward error"]]
  strict <- requirement(sources = c("theorem", "construction", "computation"),
                        guarantees = c("identity", "deterministic_bound"))
  expect_true(evidence_satisfies(ev_computed, strict))

  ## The same operator behind a callback declares it instead, and the bound is
  ## then only as good as the declaration. This is the section 5.3 laundering
  ## case reaching the certificate.
  n <- 20
  declared <- eigs(linop(function(X) M %*% X, adjoint = function(X) M %*% X,
                         dim = c(n, n), properties = c(hermitian = TRUE)),
                   k = 2, which = "largest_algebraic")
  ev_declared <- declared$certificate$evidence[["forward error"]]
  expect_false(evidence_satisfies(ev_declared, strict))
  expect_true(evidence_satisfies(
    ev_declared, requirement(sources = c("theorem", "user_declaration"))))
})

test_that("the arithmetic floor is load-bearing, on byte-identical iterates", {
  ## A tolerance below what the arithmetic can deliver. The iteration never sees
  ## floor_const, so both runs produce the same pairs and only the certificate
  ## differs, which is the whole claim.
  M <- spd_prescribed(20, seq(10, 1, length.out = 20), seed = 4L)
  A <- as_spd_linop(M)
  with_floor <- eigs(A, k = 2, which = "largest_algebraic", tol = 1e-17, maxit = 60L)
  without <- eigs(A, k = 2, which = "largest_algebraic", tol = 1e-17, maxit = 60L,
                  floor_const = 0)

  expect_identical(with_floor$values, without$values)
  expect_identical(with_floor$vectors, without$vectors)
  expect_equal(cert_status(with_floor$certificate, "backward error"), "qualified")
  expect_equal(cert_status(without$certificate, "backward error"), "fail")
})

test_that("a budget too small comes back as a fail certificate, not an error", {
  A <- laplacian_1d(80)
  fit <- eigs(A, k = 4, which = "smallest_algebraic", ncv = 8L, maxit = 6L)
  expect_lt(fit$nconv, 4L)
  expect_equal(cert_status(fit$certificate, "convergence"), "fail")
  expect_equal(fit$certificate$overall, "fail")
  ## and it still returns finite pairs with their measured residuals, at most k
  ## of them and never padded to k with anything no measurement produced
  expect_lte(length(fit$values), 4L)
  expect_gte(length(fit$values), 1L)
  expect_true(all(is.finite(fit$values)))
  detail <- fit$certificate$checks$detail[fit$certificate$checks$check == "convergence"]
  expect_match(detail, "of 4 requested pairs")
})

test_that("eigs refuses what it is not defined for", {
  n <- 10
  herm <- laplacian_1d(n)

  expect_error(eigs(linop(rmat(6, 4)), k = 2), "square")
  expect_error(eigs(linop(rmat(5, 5)), k = 2), "requires the operator to be hermitian")
  expect_error(eigs(herm, k = 2, B = linop_eye(n)), "does not take B")
  expect_error(eigs(herm, k = 2, sigma = 1, which = "smallest"), "cannot both be given")
  expect_error(eigs(herm, k = 2, sigma = 1i), "single finite real")
  expect_error(eigs(herm, k = n + 1L), "exceeds")
  expect_error(eigs(herm, k = 0), "positive integer")
  expect_error(eigs(herm, k = 2, ncv = 2L), "must exceed k")
  expect_error(eigs(herm, k = 2, which = "biggest"), "does not take which")
  expect_error(eigs(herm, k = 2, method = "arnoldi"), "'auto' or 'lanczos'")
  expect_error(eigs(herm, k = 2, v0 = rep(1, n + 1L)), "it has to be")
})

test_that("RSpectra's short forms name the same four requests", {
  A <- shifted_laplacian_1d(30, 1)
  for (pair in list(c("largest", "LM"), c("smallest", "SM"),
                    c("largest_algebraic", "LA"), c("smallest_algebraic", "SA"))) {
    long <- eigs(A, k = 2, which = pair[1L])
    short <- eigs(A, k = 2, which = pair[2L])
    expect_identical(long$values, short$values, info = pair[1L])
  }
})

test_that("the caller's random stream is not moved", {
  set.seed(99)
  before <- .Random.seed
  eigs(laplacian_1d(20), k = 2, which = "largest_algebraic")
  expect_identical(.Random.seed, before)
})
