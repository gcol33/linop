## ||A|| exists for one reason: the arithmetic floor S0.6 requires. What matters
## about it is therefore not only the number but the evidence attached, since a
## certificate line inherits its guarantee from whichever route produced it.

norm_of <- function(A, ...) linop:::norm2(A, ...)

test_that("structural rules give the spectral norm exactly", {
  expect_equal(norm_of(linop_eye(7))$value, 1)
  expect_equal(norm_of(linop_scaling(c(2, -5, 3)))$value, 5)
  expect_equal(norm_of(linop_scaling(c(1 + 1i, 3)))$value, 3)

  ## and say so: an exact route carries an identity guarantee all the way down
  for (A in list(linop_eye(7), linop_scaling(c(2, -5, 3)))) {
    expect_true(evidence_satisfies(norm_of(A)$evidence,
                                   requirement(guarantees = "identity")))
  }
})

test_that("a scale multiplies the norm and a view preserves it", {
  D <- linop_scaling(c(2, -5, 3))
  expect_equal(norm_of(3 * D)$value, 15)
  expect_equal(norm_of((-2) * D)$value, 10)

  M <- rmat(6, 4, seed = 11)
  A <- linop(M)
  base <- norm_of(A)$value
  expect_equal(norm_of(t(A))$value, base)
  expect_equal(norm_of(adjoint(A))$value, base)
  expect_equal(norm_of(Conj(A))$value, base)
  expect_equal(base, max(svd(M)$d))
})

test_that("an estimate underneath a structural rule stays visible at the top", {
  ## The scale rule is exact, but it is exact about an estimate. A composite that
  ## reported a bare identity here would launder the estimate, which is the same
  ## failure evidence_satisfies() recursion exists to prevent for capabilities.
  n <- 60
  A <- laplacian_1d(n)
  est <- norm_of(A)
  expect_equal(est$method, "power")
  expect_false(evidence_satisfies(est$evidence, requirement(guarantees = "identity")))

  scaled <- norm_of(4 * A)
  expect_equal(scaled$value, 4 * est$value)
  expect_equal(scaled$evidence$source, "construction")
  expect_false(evidence_satisfies(scaled$evidence, requirement(guarantees = "identity")))
})

test_that("a stored matrix small enough to factor is exact, and a large one is not", {
  M <- rmat(20, 20, seed = 3)
  A <- linop(M)
  est <- norm_of(A)
  expect_equal(est$method, "svd")
  expect_equal(est$value, max(svd(M)$d))
  expect_true(evidence_satisfies(est$evidence, requirement(guarantees = "identity")))

  ## the same operator above the guard falls through to the iteration
  loose <- norm_of(A, exact_max = 10)
  expect_equal(loose$method, "power")
  expect_false(evidence_satisfies(loose$evidence, requirement(guarantees = "identity")))
})

test_that("the power iteration is a lower bound and converges to the norm", {
  n <- 80
  A <- laplacian_1d(n)
  truth <- max(laplacian_1d_eigenvalues(n))

  ## every iterate is ||A v|| for a unit v, so it can never exceed the norm
  for (m in c(1L, 3L, 10L, 50L)) {
    expect_lte(norm_of(A, tol = 0, maxit = m)$value, truth * (1 + 1e-12))
  }
  ## and run out it gets there
  expect_equal(norm_of(A, tol = 1e-12, maxit = 2000)$value, truth, tolerance = 1e-6)
})

test_that("an operator with no adjoint falls back to probes and says so", {
  n <- 30
  A <- linop(function(X) laplacian_1d_apply(X), dim = c(n, n))
  est <- norm_of(A)
  expect_equal(est$method, "probe")
  expect_match(est$detail, "declares no adjoint")
  expect_lte(est$value, max(laplacian_1d_eigenvalues(n)) * (1 + 1e-12))
  expect_gt(est$value, 0)
})

test_that("estimating a norm does not move the caller's random stream", {
  A <- laplacian_1d(50)
  set.seed(99); expected <- stats::runif(3)
  set.seed(99)
  invisible(norm_of(A))
  invisible(norm_of(linop(function(X) X, dim = c(20, 20))))
  expect_equal(stats::runif(3), expected)

  ## including when the caller had no stream to begin with
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  invisible(norm_of(A))
  expect_false(exists(".Random.seed", envir = globalenv()))
})

test_that("an empty operator has norm zero", {
  A <- linop(matrix(0, 0, 0))
  expect_equal(norm_of(A)$value, 0)
})
