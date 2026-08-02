## Section 5.7 and the primary acceptance criterion: a third party implements
## linop.myclass() for a storage format this project has never seen, passes
## verify(), and gets the ecosystem working with no changes to linop.

test_that("an external adapter needs only a linop() method and a node registration", {
  ## A storage format linop has never seen: a circulant, held as its first column.
  circ <- structure(list(col = c(4, 1, 0, 2)), class = "my_circulant")

  linop:::linop_register_node(
    "my_circulant",
    apply = function(op, X, mode) {
      cc <- op$args$col
      n <- length(cc)
      ## build the dense circulant once per apply; correctness, not speed
      M <- outer(seq_len(n), seq_len(n), function(i, j) cc[((i - j) %% n) + 1L])
      switch(mode, N = M %*% X, T = t(M) %*% X, C = Conj(t(M)) %*% X, R = Conj(M) %*% X)
    },
    materialize = function(op) {
      cc <- op$args$col; n <- length(cc)
      outer(seq_len(n), seq_len(n), function(i, j) cc[((i - j) %% n) + 1L])
    },
    overwrite = TRUE)

  registerS3method("linop", "my_circulant", function(x, ...) {
    n <- length(x$col)
    linop:::new_linop("my_circulant", c(n, n), "double",
                      linop:::new_caps(), list(col = x$col))
  }, envir = globalenv())

  A <- linop(circ)
  expect_s3_class(A, "linop")

  ## the conformance suite passes
  cert <- verify(A)
  expect_false(cert$overall == "fail")
  expect_length(cert$checks$check[cert$checks$status == "fail"], 0)

  ## and the whole algebra works on it, unchanged
  n <- 4
  expect_equal(dim(A), c(n, n))
  X <- rmat(n, 2, seed = 1)
  M <- as.matrix(A)
  expect_equal(A %*% X, M %*% X, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(as.matrix(crossprod(A)), crossprod(M), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(as.matrix(t(A) %*% A + 2 * linop_eye(n)), t(M) %*% M + 2 * diag(n),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_true(isTRUE(capv(adjoint(A) %*% A, "hermitian")))
})

test_that("the registry refuses a duplicate without overwrite", {
  expect_error(linop:::linop_register_node("dense", function(op, X, mode) X),
               "already registered")
  expect_silent(linop:::linop_register_node("dense", linop:::dense_apply,
                                            linop:::dense_materialize,
                                            overwrite = TRUE))
})

test_that("all v0.1 node types are registered and none of the deferred ones are", {
  have <- linop:::linop_nodes()
  v01 <- c("dense", "sparse", "fun", "identity", "diag",
           "transpose", "adjoint", "conjugate", "scale", "sum", "product")
  expect_true(all(v01 %in% have))

  ## Section 5.7 defers these until an external adapter has exercised the
  ## registry. If one appears, it was added from inside, which is what the
  ## deferral exists to prevent.
  deferred <- c("lowrank", "perm", "zero", "power", "inverse",
                "hstack", "vstack", "blockdiag", "kron")
  expect_equal(intersect(deferred, have), character(0))
})

test_that("a node returning the wrong shape is caught", {
  linop:::linop_register_node("bad_shape",
    apply = function(op, X, mode) X[1, , drop = FALSE], overwrite = TRUE)
  A <- linop:::new_linop("bad_shape", c(4, 4), "double", linop:::new_caps(), list())
  expect_error(A %*% rmat(4, 1, seed = 2), "returned a 1 x 1 block")
})

test_that("the matrix adapter reads structure off the data as computation", {
  S <- diag(c(1, 2, 3))
  A <- linop(S)
  expect_true(isTRUE(capv(A, "symmetric")))
  expect_true(isTRUE(capv(A, "diagonal")))
  expect_equal(linop:::cape(A, "symmetric")$source, "computation")
  expect_equal(linop:::cape(A, "symmetric")$guarantee, "identity")
})

test_that("the Matrix adapter reads structure off the class as an adapter contract", {
  skip_if_not_installed("Matrix")
  M <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 2, 3),
                            dims = c(3, 3), symmetric = TRUE)
  A <- linop(M)
  expect_true(isTRUE(capv(A, "symmetric")))
  expect_equal(linop:::cape(A, "symmetric")$source, "adapter_contract")
})

test_that("complex sparse is refused at construction rather than silently wrong", {
  skip_if_not_installed("Matrix")
  ## Matrix 1.7.5 has no concrete complex sparse class; see S0.1 section 7.
  expect_error(Matrix::Matrix(zmat(3, 3, seed = 3)))
})

test_that("a callback operator without an adjoint refuses adjoint modes clearly", {
  n <- 4
  A <- linop(function(X) 2 * X, dim = c(n, n))
  expect_error(t(A) %*% rmat(n, 1, seed = 4), "no adjoint")
  expect_error(adjoint(A) %*% rmat(n, 1, seed = 5), "no adjoint")
  ## but the forward direction is fine
  expect_equal(A %*% rep(1, n), rep(2, n), tolerance = 1e-12)
})
