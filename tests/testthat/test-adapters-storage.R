## Phase 4's adapters: out-of-core and device storage through the linop()
## generic, each skipped when its package is absent. Every mode against the
## dense reference, the refusal of single precision, and a solve and a spectrum
## through each.

modes_vs_dense <- function(A, M) {
  cplx <- is.complex(M)
  gen <- function(r, c) if (cplx) zmat(r, c) else rmat(r, c)
  X <- gen(ncol(M), 2); Y <- gen(nrow(M), 2)
  max(Mod(linop:::linop_apply(A, X, "N") - M %*% X),
      Mod(linop:::linop_apply(A, Y, "T") - t(M) %*% Y),
      Mod(linop:::linop_apply(A, Y, "C") - Conj(t(M)) %*% Y),
      Mod(linop:::linop_apply(A, X, "R") - Conj(M) %*% X))
}

test_that("a DelayedMatrix applies by DelayedArray's product in every mode", {
  skip_if_not_installed("DelayedArray")
  for (cplx in c(FALSE, TRUE)) {
    M <- if (cplx) zmat(12, 7, seed = 1) else rmat(12, 7, seed = 1)
    A <- linop(DelayedArray::DelayedArray(M))
    expect_equal(A$dtype, if (cplx) "complex" else "double")
    expect_lt(modes_vs_dense(A, M), 1e-13)
    expect_false(verify(A)$overall == "fail")
  }
  ## nothing is computed from the entries; the caller declares
  S <- crossprod(rmat(8, 8, seed = 2)) + diag(8)
  D <- linop(DelayedArray::DelayedArray(S), properties = c(hermitian = TRUE))
  expect_true(is.na(capv(linop(DelayedArray::DelayedArray(S)), "hermitian")))
  fit <- eigs(D, k = 2, which = "largest_algebraic", method = "lanczos")
  expect_equal(fit$values, eigen(S, symmetric = TRUE)$values[1:2], tolerance = 1e-10)
})

test_that("a torch tensor applies on its device in every mode, double precision only", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "libtorch is not installed")
  for (cplx in c(FALSE, TRUE)) {
    M <- if (cplx) zmat(9, 5, seed = 3) else rmat(9, 5, seed = 3)
    Tm <- torch::torch_tensor(M, dtype = if (cplx) torch::torch_cdouble() else torch::torch_float64())
    A <- linop(Tm)
    expect_lt(modes_vs_dense(A, M), 1e-13)
    expect_false(verify(A)$overall == "fail")
  }
  expect_error(linop(torch::torch_tensor(rmat(3, 3, seed = 4))), "float64 or complex128")
  expect_error(linop(torch::torch_tensor(array(1, c(2, 2, 2)), dtype = torch::torch_float64())),
               "two-dimensional")
  S <- crossprod(rmat(6, 6, seed = 5)) + diag(6)
  A <- linop(torch::torch_tensor(S, dtype = torch::torch_float64()),
             properties = c(hermitian = TRUE, positive_definite = TRUE))
  x <- solve(A, rep(1, 6), method = "cg", tol = 1e-12)
  expect_equal(as.numeric(x), solve(S, rep(1, 6)), tolerance = 1e-9)
})

test_that("a big.matrix, file-backed too, applies a block of rows at a time", {
  skip_if_not_installed("bigmemory")
  M <- rmat(37, 6, seed = 6)
  for (backed in c(FALSE, TRUE)) {
    B <- if (backed) {
      dir <- tempfile(); dir.create(dir)
      bigmemory::as.big.matrix(M, backingfile = "m.bin", descriptorfile = "m.desc",
                               backingpath = dir)
    } else bigmemory::as.big.matrix(M)
    A <- linop(B)
    expect_lt(modes_vs_dense(A, M), 1e-13)
    expect_lt(modes_vs_dense(A, M + 0i), 1e-13)
  }
  ## several row blocks, by shrinking the block size
  old <- linop:::BIGMATRIX_BLOCK_ENTRIES
  ns <- asNamespace("linop")
  unlockBinding("BIGMATRIX_BLOCK_ENTRIES", ns)
  assign("BIGMATRIX_BLOCK_ENTRIES", 30, envir = ns)
  on.exit({ assign("BIGMATRIX_BLOCK_ENTRIES", old, envir = ns)
            lockBinding("BIGMATRIX_BLOCK_ENTRIES", ns) }, add = TRUE)
  A <- linop(bigmemory::as.big.matrix(M))
  expect_lt(modes_vs_dense(A, M), 1e-13)
  fit <- svds(linop(bigmemory::as.big.matrix(M)), k = 2)
  expect_equal(fit$d, svd(M)$d[1:2], tolerance = 1e-9)
})
