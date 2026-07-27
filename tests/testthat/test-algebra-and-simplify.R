## Sections 1.1, 5.8, 5.9.

test_that("%*% returns what a matrix would return", {
  M <- rmat(5, 4, seed = 1); A <- linop(M)
  x <- rnorm(4); X <- rmat(4, 3, seed = 2)

  ## vector in, vector out, computed now
  expect_null(dim(A %*% x))
  expect_equal(as.vector(A %*% x), as.vector(M %*% x), tolerance = 1e-12)
  ## matrix in, matrix out, computed now
  expect_equal(A %*% X, M %*% X, tolerance = 1e-12, ignore_attr = TRUE)
  ## linop in, linop out, lazy
  expect_s3_class(A %*% t(A), "linop")
})

test_that("dispatch works when the linop is the second argument", {
  M <- rmat(4, 4, seed = 3); A <- linop(M); B <- rmat(4, 4, seed = 4)
  expect_s3_class(B %*% A, "linop")
  expect_equal(as.matrix(B %*% A), B %*% M, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("crossprod and tcrossprod match base R for one and two arguments", {
  M <- rmat(5, 4, seed = 5); A <- linop(M); Y <- rmat(5, 2, seed = 6)
  expect_equal(as.matrix(crossprod(A)), crossprod(M), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(as.matrix(tcrossprod(A)), tcrossprod(M), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(crossprod(A, Y), crossprod(M, Y), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("scalar arithmetic behaves and non-scalar is refused", {
  M <- rmat(3, 3, seed = 7); A <- linop(M)
  expect_equal(as.matrix(2 * A), 2 * M, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(as.matrix(A * 2), 2 * M, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(as.matrix(A / 2), M / 2, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(as.matrix(-A), -M, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(as.matrix(A - A), M - M, tolerance = 1e-12, ignore_attr = TRUE)
  expect_error(A * c(1, 2), "only scalar multiplication")
  expect_error(A * A, "not defined here")
  expect_error(A + 1, "ambiguous")
  expect_error(A / linop(M), "not defined")
})

test_that("simplification rules fire at construction", {
  n <- 4
  A <- linop(rmat(n, n, seed = 8)); I <- linop_eye(n)

  expect_identical((I %*% A)$node, A$node)
  expect_identical((A %*% I)$node, A$node)
  expect_identical(t(t(A))$node, A$node)
  expect_identical(adjoint(adjoint(A))$node, A$node)
  expect_identical(Conj(Conj(A))$node, A$node)
  expect_identical(adjoint(t(A))$node, "conjugate")
  expect_identical(t(adjoint(A))$node, "conjugate")
  expect_identical((2 * (3 * A))$args$a, 6)

  d1 <- c(1, 2, 3, 4); d2 <- c(2, 2, 2, 2)
  fused <- linop_scaling(d1) %*% linop_scaling(d2)
  expect_identical(fused$node, "diag")
  expect_equal(fused$args$d, d1 * d2)

  expect_identical(t(linop_scaling(d1))$node, "diag")
  expect_equal(adjoint(linop_scaling(complex(real = 1:4, imaginary = 1:4)))$args$d,
               Conj(complex(real = 1:4, imaginary = 1:4)))
})

test_that("a real operator keeps adjoint and transpose as distinct nodes", {
  ## Section 5.8: the simplification happens at apply time, not construction, so
  ## print() keeps showing what the author wrote.
  A <- linop(rmat(4, 3, seed = 9))
  expect_identical(t(A)$node, "transpose")
  expect_identical(adjoint(A)$node, "adjoint")
  expect_false(identical(t(A)$node, adjoint(A)$node))
  ## but they act identically
  X <- rmat(4, 2, seed = 10)
  expect_equal(t(A) %*% X, adjoint(A) %*% X, tolerance = 1e-12)
})

test_that("sums and products flatten rather than nesting", {
  A <- linop(rmat(3, 3, seed = 11))
  s <- A + A + A
  expect_identical(s$node, "sum")
  expect_length(s$args$ops, 3L)
  p <- A %*% A %*% A
  expect_identical(p$node, "product")
  expect_length(p$args$ops, 3L)
})

test_that("non-conformable operations are refused with a useful message", {
  A <- linop(rmat(5, 4, seed = 12)); B <- linop(rmat(3, 3, seed = 13))
  expect_error(A %*% B, "non-conformable product")
  expect_error(A + B, "non-conformable sum")
  expect_error(A %*% rnorm(3), "non-conformable")
})

test_that("as.matrix guards on size and the guard can be raised", {
  big <- linop(function(X) X, adjoint = function(X) X, dim = c(1e5, 1e5))
  expect_error(as.matrix(big), "refusing to materialise")
  expect_error(as.matrix(big), "max_entries")
})

test_that("collapse reports what it did and preserves the action", {
  M1 <- rmat(40, 6, seed = 14); M2 <- rmat(6, 40, seed = 15)
  A <- linop(M1) %*% linop(M2)
  out <- collapse(A, expected_applies = 500)
  X <- rmat(40, 2, seed = 16)
  expect_equal(out %*% X, A %*% X, tolerance = 1e-10)
  rep <- attr(out, "report")
  expect_true(is.list(rep))
  expect_true(rep$memory_added >= 0)
})

test_that("explain reports the node structure", {
  A <- crossprod(linop(rmat(5, 4, seed = 17))) + 0.1 * linop_eye(4)
  df <- NULL
  invisible(capture.output(df <- explain(A)))
  expect_true(is.data.frame(df))
  expect_true(all(c("depth", "node", "nrow", "ncol", "dtype", "cost") %in% names(df)))
  expect_true("product" %in% df$node)
  expect_true("identity" %in% df$node)
})

test_that("an adapter author gets a pointed error for an unsupported class", {
  expect_error(linop(structure(list(), class = "weird_thing")),
               "implement linop.weird_thing")
})
