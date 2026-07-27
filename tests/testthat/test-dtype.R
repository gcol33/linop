## Gate 1: the dtype promotion table, exhaustively. Section 5.5.

test_that("the promotion table is exhaustive and correct", {
  grid <- expand.grid(a = c("double", "complex"), b = c("double", "complex"),
                      stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    a <- grid$a[i]; b <- grid$b[i]
    want <- if (a == "complex" || b == "complex") "complex" else "double"
    expect_equal(linop:::promote(a, b), want, info = sprintf("promote(%s, %s)", a, b))
  }
  expect_equal(nrow(grid), 4L)
})

test_that("promotion is commutative and idempotent", {
  for (a in c("double", "complex")) for (b in c("double", "complex")) {
    expect_equal(linop:::promote(a, b), linop:::promote(b, a))
  }
  for (a in c("double", "complex")) expect_equal(linop:::promote(a, a), a)
})

test_that("promote rejects unknown types", {
  expect_error(linop:::promote("single", "double"), "dtype must be one of")
  expect_error(linop:::promote("double", "integer"), "dtype must be one of")
})

test_that("a real operator applied to a complex block promotes the result", {
  A <- linop(rmat(4, 3, seed = 1))
  expect_equal(dtype(A), "double")
  X <- zmat(3, 2, seed = 2)
  Y <- A %*% X
  expect_true(is.complex(Y))
  expect_equal(Y, as.matrix(A) %*% X, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("no path silently downcasts", {
  expect_error(linop:::cast_block(zmat(3, 1, seed = 3), "double"), "refusing to downcast")
})

test_that("expression dtype follows the lattice", {
  R <- linop(rmat(3, 3, seed = 4))
  Z <- linop(zmat(3, 3, seed = 5))
  expect_equal(dtype(R + R), "double")
  expect_equal(dtype(R + Z), "complex")
  expect_equal(dtype(Z %*% R), "complex")
  expect_equal(dtype(2 * R), "double")
  expect_equal(dtype((1 + 1i) * R), "complex")
  expect_equal(dtype(t(Z)), "complex")
})

test_that("a complex operator refuses to lose its imaginary part on materialisation", {
  Z <- linop(zmat(3, 3, seed = 6))
  M <- as.matrix(Z)
  expect_true(is.complex(M))
})
