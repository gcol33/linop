## Section 7.2, svds(). Recovery against closed-form singular values, the two
## halves of a triplet, certificate coverage across seeds, and the refusals.

## The backward error of the augmented eigenpair, recomputed from the operator
## against a ||A|| the fixture knows in closed form. The run stops on that
## quantity relative to a lower bound on ||A||, so it is what tol means; an
## absolute residual compared against tol would be off by a factor of ||A||.
triplet_backward <- function(A, d, u, v, anorm) {
  ru <- linop_apply(A, v, "N") - u * rep(d, each = nrow(u))
  rv <- linop_apply(A, u, "C") - v * rep(d, each = nrow(v))
  max(sqrt(colSums(Mod(ru)^2) + colSums(Mod(rv)^2)) / sqrt(2)) / anorm
}

test_that("svds recovers a prescribed singular spectrum, both halves of each triplet", {
  n <- 20
  sigma <- c(40, 20, 12, seq(4, 0.5, length.out = n - 3L))
  f <- lsq_prescribed(30, n, sigma, seed = 2L)
  A <- linop(f$A)

  fit <- svds(A, k = 3)
  expect_equal(fit$nconv, 3L)
  expect_equal(fit$d, sigma[1:3], tolerance = 1e-10)
  ## A v = sigma u and A^H u = sigma v both, which is what makes it a triplet
  ## rather than an eigenpair of the Gram operator.
  expect_lte(triplet_backward(A, fit$d, fit$u, fit$v, max(sigma)), 1e-10)

  ## and the singular vectors themselves, up to a sign
  for (j in 1:3) {
    expect_lt(min(max(abs(fit$v[, j] - f$V[, j])), max(abs(fit$v[, j] + f$V[, j]))),
              1e-7)
  }
})

test_that("svds recovers the matrix-free rank-deficient fixture", {
  n <- 50
  A <- diff_1d(n)
  truth <- rev(sort(diff_1d_singular_values(n)))

  fit <- svds(A, k = 3)
  expect_equal(fit$nconv, 3L)
  expect_equal(fit$d, truth[1:3], tolerance = 1e-10)
  expect_lte(triplet_backward(A, fit$d, fit$u, fit$v, max(truth)), 1e-10)
})

test_that("the two bases are orthonormal separately, not only after augmentation", {
  ## The certificate's orthogonality row is measured on [u; v] / sqrt(2), where a
  ## deviation in U could in principle be cancelled by one in V. That the two are
  ## orthonormal on their own is asserted here instead.
  n <- 20
  sigma <- seq(10, 1, length.out = n)
  f <- lsq_prescribed(28, n, sigma, seed = 6L)
  fit <- svds(linop(f$A), k = 4)
  expect_lt(max(abs(crossprod(fit$u) - diag(4))), 1e-12)
  expect_lt(max(abs(crossprod(fit$v) - diag(4))), 1e-12)
  expect_equal(cert_status(fit$certificate, "orthogonality"), "pass")
})

test_that("svds recovers a complex spectrum and keeps complex vectors", {
  n <- 12
  sigma <- c(20, 9, 5, seq(2, 0.4, length.out = n - 3L))
  f <- lsq_prescribed(18, n, sigma, seed = 4L, dtype = "complex")
  A <- linop(f$A)

  fit <- svds(A, k = 3)
  expect_equal(fit$nconv, 3L)
  expect_equal(fit$d, sigma[1:3], tolerance = 1e-10)
  expect_true(is.complex(fit$u) && is.complex(fit$v))
  ## a singular value is a non-negative real whatever the operator is
  expect_false(is.complex(fit$d))
  expect_true(all(fit$d >= 0))
  expect_lte(triplet_backward(A, fit$d, fit$u, fit$v, max(sigma)), 1e-10)
})

test_that("the smallest singular values are reachable and are not the largest relabelled", {
  n <- 16
  sigma <- seq(20, 1, length.out = n)
  f <- lsq_prescribed(24, n, sigma, seed = 8L)
  A <- linop(f$A)

  lo <- svds(A, k = 2, which = "smallest", ncv = 16L, maxit = 2000L)
  hi <- svds(A, k = 2, which = "largest")
  expect_equal(lo$nconv, 2L)
  expect_equal(sort(lo$d), sort(sigma[(n - 1):n]), tolerance = 1e-8)
  expect_equal(sort(hi$d), sort(sigma[1:2]), tolerance = 1e-10)
  expect_gt(min(hi$d) - max(lo$d), 1)
})

test_that("the certificate is the eigen certificate of the augmented operator", {
  n <- 15
  sigma <- c(9, 4, 2, seq(1, 0.3, length.out = n - 3L))
  f <- lsq_prescribed(20, n, sigma, seed = 9L)
  fit <- svds(linop(f$A), k = 2)

  expect_equal(fit$certificate$subject, "svd")
  expect_setequal(fit$certificate$checks$check,
                  c("arithmetic floor", "residual", "orthogonality",
                    "backward error", "target identity", "convergence",
                    "forward error"))
  expect_equal(cert_status(fit$certificate, "target identity"), "not_checked")

  ## The augmented operator is hermitian by construction with nothing under it,
  ## so the forward bound on a singular value rests on no declaration at all --
  ## which the same bound out of eigs() cannot say unless the operator's own
  ## symmetry was established rather than declared.
  ev <- fit$certificate$evidence[["forward error"]]
  expect_true(evidence_satisfies(
    ev, requirement(sources = c("theorem", "construction"),
                    guarantees = c("identity", "deterministic_bound"),
                    min_confidence = 1)))
})

test_that("the forward-error bound contains the true error over 20 seeds", {
  n <- 12
  sigma <- c(15, 7, 6.8, seq(3, 0.5, length.out = n - 3L))
  worst_slack <- 0
  for (s in 1:20) {
    f <- lsq_prescribed(18, n, sigma, seed = s)
    fit <- svds(linop(f$A), k = 3, seed = s)
    bound <- fit$certificate$values$forward_bound
    err <- vapply(fit$d, function(x) min(abs(x - sigma)), numeric(1))
    nA <- fit$certificate$values$norm
    expect_true(all(err <= bound + 8 * .Machine$double.eps * nA),
                info = sprintf("seed %d: worst err %.3e against bound %.3e",
                               s, max(err), max(bound)))
    worst_slack <- max(worst_slack, max(err - bound))
  }
  expect_lt(worst_slack, 1e-12)
})

test_that("||A|| for the augmented operator is inherited structurally, not re-estimated", {
  ## Its spectrum is the singular values of A with both signs, so the norm is
  ## exactly ||A||. A stored matrix small enough to factor gives that exactly, and
  ## the structural rule over it has to stay exact.
  f <- lsq_prescribed(12, 8, seq(6, 1, length.out = 8), seed = 3L)
  fit <- svds(linop(f$A), k = 2)
  expect_equal(fit$certificate$values$norm, max(f$sigma), tolerance = 1e-12)
  row <- fit$certificate$checks
  expect_equal(row$source[row$check == "arithmetic floor"], "construction")
  expect_equal(row$guarantee[row$check == "arithmetic floor"], "identity")
})

test_that("the arithmetic floor is load-bearing, on byte-identical iterates", {
  f <- lsq_prescribed(16, 8, seq(6, 1, length.out = 8), seed = 5L)
  A <- linop(f$A)
  with_floor <- svds(A, k = 2, tol = 1e-17, maxit = 40L)
  without <- svds(A, k = 2, tol = 1e-17, maxit = 40L, floor_const = 0)

  expect_identical(with_floor$d, without$d)
  expect_identical(with_floor$u, without$u)
  expect_identical(with_floor$v, without$v)
  expect_equal(cert_status(with_floor$certificate, "backward error"), "qualified")
  expect_equal(cert_status(without$certificate, "backward error"), "fail")
})

test_that("a budget too small comes back as a fail certificate, not an error", {
  ## A run whose subspace never reached k wide returns what it found and nothing
  ## else. Padding to k would put a number in the result that no measurement
  ## produced, and the convergence line counts against the request rather than
  ## against what survived, so a short result still reads as incomplete.
  A <- diff_1d(80)
  fit <- svds(A, k = 4, ncv = 8L, maxit = 6L)
  expect_lt(fit$nconv, 4L)
  expect_equal(cert_status(fit$certificate, "convergence"), "fail")
  expect_lte(length(fit$d), 4L)
  expect_gte(length(fit$d), 1L)
  expect_true(all(is.finite(fit$d)))
  detail <- fit$certificate$checks$detail[fit$certificate$checks$check == "convergence"]
  expect_match(detail, "of 4 requested pairs")
})

test_that("svds refuses what it is not defined for", {
  n <- 12
  A <- linop(rmat(16, n, seed = 1L))
  no_adjoint <- linop(function(X) X[-1L, , drop = FALSE], dim = c(n - 1L, n))

  expect_error(svds(no_adjoint, k = 2), "needs the adjoint")
  expect_error(svds(A, k = n + 1L), "exceeds")
  expect_error(svds(A, k = 2, which = "largest_algebraic"), "does not take which")
  expect_error(svds(A, k = 2, ncv = 2L), "must exceed k")
  expect_error(svds(A, k = 2, method = "lanczos"), "'auto' or 'golub-kahan'")
})

test_that("svds is not eigs on the Gram operator", {
  ## The separation the bidiagonalisation exists for: a singular value at 1e-9
  ## relative to the largest sits at 1e-18 in A^H A and is below the rounding
  ## level of its largest eigenvalue, so the squared route cannot resolve it and
  ## this one does.
  n <- 8
  sigma <- c(1, 1e-9, seq(0.9, 0.2, length.out = n - 2L))
  f <- lsq_prescribed(12, n, sigma, seed = 12L)
  A <- linop(f$A)

  fit <- svds(A, k = n, which = "smallest", ncv = 8L, maxit = 4000L)
  expect_equal(min(fit$d), 1e-9, tolerance = 1e-6)

  gram <- eigen(crossprod(f$A), symmetric = TRUE, only.values = TRUE)$values
  expect_gt(abs(sqrt(max(gram[gram > 0])) - 1) + abs(min(gram)), 0)
  ## the squared route loses it: the smallest Gram eigenvalue is at the rounding
  ## level of the largest rather than at 1e-18
  expect_gt(min(abs(gram)), 1e-24)
})

test_that("the caller's random stream is not moved", {
  A <- linop(rmat(10, 6, seed = 2L))
  set.seed(101)
  before <- .Random.seed
  svds(A, k = 2)
  expect_identical(.Random.seed, before)
})
