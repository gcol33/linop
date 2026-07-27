## Gate 1: the conformance suite passes on every node type, real and complex
## where the backing storage admits it.

test_that("verify() passes on every node type", {
  for (f in all_node_fixtures()) {
    cert <- verify(f$A)
    expect_false(cert$overall == "fail",
                 info = sprintf("node fixture '%s' failed conformance", f$label))
    failed <- cert$checks$check[cert$checks$status == "fail"]
    expect_length(failed, 0)
  }
})

test_that("every node type materialises to its dense reference", {
  for (f in all_node_fixtures()) {
    if (is.null(f$M)) next
    expect_equal(as.matrix(f$A), unname(as.matrix(f$M)), tolerance = 1e-12,
                 ignore_attr = TRUE,
                 info = sprintf("node fixture '%s'", f$label))
  }
})

test_that("apply agrees with the dense reference in all four modes", {
  for (f in all_node_fixtures()) {
    if (is.null(f$M)) next
    A <- f$A; M <- as.matrix(f$M)
    m <- nrow(M); n <- ncol(M)
    Xn <- zmat(n, 3, seed = 7); Xm <- zmat(m, 3, seed = 8)
    expect_equal(linop:::linop_apply(A, Xn, "N"), M %*% Xn, tolerance = 1e-10,
                 ignore_attr = TRUE, info = paste(f$label, "mode N"))
    expect_equal(linop:::linop_apply(A, Xm, "T"), t(M) %*% Xm, tolerance = 1e-10,
                 ignore_attr = TRUE, info = paste(f$label, "mode T"))
    expect_equal(linop:::linop_apply(A, Xm, "C"), Conj(t(M)) %*% Xm, tolerance = 1e-10,
                 ignore_attr = TRUE, info = paste(f$label, "mode C"))
    expect_equal(linop:::linop_apply(A, Xn, "R"), Conj(M) %*% Xn, tolerance = 1e-10,
                 ignore_attr = TRUE, info = paste(f$label, "mode R"))
  }
})

test_that("verify() catches a wrong adjoint", {
  n <- 5
  ## missing conjugation: the single most common third-party bug
  bad <- linop(function(X) (1 + 2i) * X,
               adjoint = function(X) (1 + 2i) * X,
               dim = c(n, n), dtype = "complex")
  cert <- verify(bad)
  expect_equal(cert$overall, "fail")
  expect_true("adjoint consistency" %in% cert$checks$check[cert$checks$status == "fail"])
})

test_that("verify() catches a nonlinear operator wrapped as a linop", {
  n <- 4
  ## an adaptive solve dressed as an operator: scaling depends on the input
  bad <- linop(function(X) X / max(1e-12, max(abs(X))), dim = c(n, n))
  cert <- verify(bad)
  expect_true("linearity" %in% cert$checks$check[cert$checks$status == "fail"])
})

test_that("verify() catches a block apply that disagrees with single columns", {
  n <- 4
  bad <- linop(function(X) { X[, 1L] <- 0; X }, adjoint = function(X) X, dim = c(n, n))
  cert <- verify(bad)
  expect_true("block consistency" %in% cert$checks$check[cert$checks$status == "fail"])
})

test_that("verify() catches an impure operator", {
  n <- 4
  counter <- 0
  bad <- linop(function(X) { counter <<- counter + 1; X * counter },
               adjoint = function(X) X, dim = c(n, n))
  cert <- verify(bad)
  expect_true("purity" %in% cert$checks$check[cert$checks$status == "fail"])
})

test_that("a contradicted declaration is a failure, not a warning", {
  M <- rmat(4, 4, seed = 3)
  M[1, 2] <- M[2, 1] + 1        # definitely not symmetric
  bad <- linop_dense(M, properties = c(hermitian = TRUE), check = FALSE)
  cert <- verify(bad)
  expect_equal(cert$overall, "fail")
  expect_true("declared capabilities" %in% cert$checks$check[cert$checks$status == "fail"])
})

test_that("an operator without an adjoint reports not_checked rather than failing", {
  n <- 4
  A <- linop(function(X) 2 * X, dim = c(n, n))
  cert <- verify(A)
  expect_false(cert$overall == "fail")
  expect_equal(cert$checks$status[cert$checks$check == "adjoint consistency"], "not_checked")
})
