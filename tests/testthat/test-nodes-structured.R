## The structured node types: every apply mode against base R, the
## shapes they refuse, what they compute with, and recovery of closed forms
## through the solvers and eigensolvers when an operator is built from them.

modes_agree <- function(A, M) {
  cplx <- is.complex(M) || A$dtype == "complex"
  gen <- function(r, c) if (cplx) zmat(r, c) else rmat(r, c)
  X <- gen(ncol(M), 3); Y <- gen(nrow(M), 3)
  max(Mod(as.matrix(A) - M),
      Mod(linop:::linop_apply(A, X, "N") - M %*% X),
      Mod(linop:::linop_apply(A, Y, "T") - t(M) %*% Y),
      Mod(linop:::linop_apply(A, Y, "C") - Conj(t(M)) %*% Y),
      Mod(linop:::linop_apply(A, X, "R") - Conj(M) %*% X)) / max(Mod(M), 1)
}

test_that("every structured node agrees with base R in all four modes, over seeds", {
  for (s in 1:20) {
    set.seed(s)
    A <- zmat(3, 2); B <- rmat(2, 4); S <- zmat(4, 4) + diag(6, 4)
    expect_lt(modes_agree(linop_kron(A, B), kronecker(A, B)), 1e-13)
    M <- matrix(0i, 5, 6); M[1:3, 1:2] <- A; M[4:5, 3:6] <- B
    expect_lt(modes_agree(linop_blockdiag(A, B), M), 1e-13)
    C3 <- zmat(3, 3); R1 <- rmat(1, 4)
    expect_lt(modes_agree(cbind(linop(A), C3), cbind(A, C3)), 1e-13)
    expect_lt(modes_agree(rbind(linop(B), R1), rbind(B, R1)), 1e-13)
    U <- zmat(5, 2); V <- rmat(4, 2); d <- c(2, 1i)
    expect_lt(modes_agree(linop_lowrank(U, V, d), U %*% diag(d) %*% t(V)), 1e-13)
    p <- sample.int(6)
    expect_lt(modes_agree(linop_perm(p), diag(6)[p, ]), 1e-15)
    expect_lt(modes_agree(linop_power(S, 3), S %*% S %*% S), 1e-12)
    expect_lt(modes_agree(solve(linop(S)), solve(S)), 1e-12)
    expect_lt(modes_agree(linop_power(S, -2), solve(S %*% S)), 1e-11)
  }
})

test_that("a Kronecker product applies without forming either factor", {
  n <- 30
  calls <- 0
  L1 <- linop(function(X) { calls <<- calls + 1; laplacian_1d_apply(X) },
              adjoint = laplacian_1d_apply, dim = c(n, n))
  K <- linop_kron(L1, linop_eye(n))
  X <- rmat(n * n, 4, seed = 1)
  Y <- K %*% X
  ## one call of the factor for the whole block of four columns
  expect_equal(calls, 1)
  L <- as.matrix(laplacian_1d(n))
  expect_lt(max(abs(Y - kronecker(L, diag(n)) %*% X)), 1e-12)
})

test_that("the 2-D Laplacian built as a Kronecker sum has its closed-form spectrum", {
  ## L_n (x) I_m + I_n (x) L_m has eigenvalues lambda_i + mu_j and eigenvectors
  ## the products of the 1-D sine vectors. On a square grid lambda_i + lambda_j
  ## repeats for i != j, and a Krylov space reaches one vector of a repeated
  ## eigenvalue; a 20 x 13 grid has no such repeat at the top. Both factors are
  ## established hermitian by computation, and the sum by construction.
  n <- 20; m <- 13
  L2 <- linop_kron(linop(as.matrix(laplacian_1d(n))), linop_eye(m)) +
        linop_kron(linop_eye(n), linop(as.matrix(laplacian_1d(m))))
  expect_true(isTRUE(capv(L2, "hermitian")))
  expect_equal(linop:::cape(L2, "hermitian")$source, "construction")
  truth <- sort(outer(laplacian_1d_eigenvalues(n), laplacian_1d_eigenvalues(m), "+"),
                decreasing = TRUE)
  expect_gt(min(-diff(truth[1:5])), 1e-3)
  fit <- eigs(L2, k = 4, which = "largest_algebraic")
  expect_equal(fit$method, "lanczos")
  expect_equal(fit$values, truth[1:4], tolerance = 1e-9)
  expect_equal(certificate(fit)$overall %in% c("pass", "qualified"), TRUE)
})

test_that("a block operator solves through the package's own methods", {
  set.seed(4)
  A <- crossprod(rmat(6, 6)) + diag(6)
  B <- linop_blockdiag(linop(A, properties = c(positive_definite = TRUE)),
                       linop_scaling(c(2, 3)))
  b <- rnorm(8)
  x <- solve(B, b, tol = 1e-12)
  ref <- c(solve(A, b[1:6]), b[7:8] / c(2, 3))
  expect_equal(as.numeric(x), ref, tolerance = 1e-9)
})

test_that("solve(a) forms an exact inverse or refuses with solve(a, b) named", {
  S <- rmat(4, 4, seed = 2) + diag(5, 4)
  Si <- solve(linop(S))
  expect_equal(Si$node, "inverse")
  expect_lt(max(abs(as.matrix(Si) - solve(S))), 1e-13)
  ## structure gives exact inverses without a factorisation
  expect_equal(solve(linop_scaling(c(2, 4)))$node, "diag")
  expect_equal(as.matrix(solve(linop_scaling(c(2, 4)))), diag(c(0.5, 0.25)))
  P <- linop_perm(c(2, 3, 1))
  expect_equal(as.matrix(solve(P) %*% P), diag(3))
  expect_equal(solve(solve(linop(S)))$node, "dense")
  K <- linop_kron(linop(S), linop_scaling(c(1, 2)))
  expect_equal(solve(K)$node, "kron")
  expect_lt(max(abs(as.matrix(solve(K)) - solve(as.matrix(K)))), 1e-12)
  expect_error(solve(laplacian_1d(5)), "an iterative solve is not a linear operator")
  expect_error(solve(linop_scaling(c(1, 0))), "zero entry")
  expect_error(solve(linop(matrix(1, 3, 3))), "singular")
  expect_error(solve(linop(rmat(3, 2))), "square")
  ## a unitary declaration is the caller's word, and solve(a) takes it
  Q <- qr.Q(qr(rmat(4, 4, seed = 3)))
  U <- linop(function(X) Q %*% X, adjoint = function(X) crossprod(Q, X), dim = c(4, 4),
             properties = c(unitary = TRUE))
  expect_lt(max(abs(as.matrix(solve(U)) - t(Q))), 1e-14)
})

test_that("the norm is structural where the structure gives it", {
  A <- rmat(3, 3, seed = 1); B <- rmat(2, 2, seed = 2)
  nk <- linop:::norm2(linop_kron(linop(A), linop(B)))
  expect_equal(nk$value, norm(kronecker(A, B), "2"), tolerance = 1e-12)
  nb <- linop:::norm2(linop_blockdiag(linop(A), linop(B)))
  expect_equal(nb$value, max(norm(A, "2"), norm(B, "2")), tolerance = 1e-12)
  U <- rmat(10, 3, seed = 3); V <- rmat(8, 3, seed = 4)
  nl <- linop:::norm2(linop_lowrank(U, V, c(3, 2, 1)))
  expect_equal(nl$value, norm(U %*% diag(c(3, 2, 1)) %*% t(V), "2"), tolerance = 1e-12)
  expect_equal(linop:::norm2(linop_perm(c(2, 1, 3)))$value, 1)
  H <- linop(crossprod(A))
  expect_equal(linop:::norm2(linop_power(H, 3))$value, norm(crossprod(A), "2")^3,
               tolerance = 1e-10)
})

test_that("construction simplifies what is not a node", {
  expect_equal(linop_power(linop(rmat(3, 3, seed = 1)), 0)$node, "identity")
  A <- linop(rmat(3, 3, seed = 1))
  expect_identical(linop_power(A, 1), A)
  expect_equal(linop_power(linop_power(A, 2), 3)$args$k, 6L)
  expect_equal(linop_perm(1:4)$node, "identity")
  expect_equal(linop_kron(linop_scaling(1:2), linop_scaling(3:4))$node, "diag")
  expect_equal(linop_blockdiag(linop_scaling(1:2), linop_eye(2))$node, "diag")
  expect_equal(length(cbind(cbind(A, A), A)$args$ops), 3L)
  expect_equal(length(linop_blockdiag(linop_blockdiag(A, A), A)$args$ops), 3L)
})

test_that("what the constructors refuse", {
  expect_error(linop_perm(c(1, 1, 2)), "permutation")
  expect_error(linop_power(linop(rmat(2, 3)), 2), "square")
  expect_error(linop_power(linop(rmat(2, 2)), 1.5), "whole number")
  expect_error(linop_lowrank(rmat(4, 2), rmat(3, 3)), "same number of columns")
  expect_error(linop_lowrank(rmat(4, 2), d = 1), "one weight per column")
  expect_error(cbind(linop(rmat(2, 2)), rmat(3, 2)), "same number of rows")
  expect_error(rbind(linop(rmat(2, 2)), rmat(2, 3)), "same number of columns")
  expect_error(linop_kron(linop_jacobi(diagonal = 1), linop_eye(2)), "finite dimensions")
  expect_error(linop_blockdiag(), "at least one block")
  expect_error(linop_kron("a", linop_eye(2)), "linops or matrices")
})

test_that("the structured nodes pass the conformance suite", {
  ops <- list(linop_kron(zmat(2, 3, seed = 1), rmat(2, 2, seed = 2)),
              linop_blockdiag(rmat(2, 2, seed = 3), zmat(1, 3, seed = 4)),
              cbind(linop(rmat(3, 2, seed = 5)), rmat(3, 1, seed = 6)),
              rbind(linop(zmat(2, 3, seed = 7)), rmat(1, 3, seed = 8)),
              linop_lowrank(zmat(4, 2, seed = 9), rmat(3, 2, seed = 10)),
              linop_perm(c(4, 2, 1, 3)),
              linop_power(linop(zmat(3, 3, seed = 11)), 2),
              solve(linop(rmat(3, 3, seed = 12) + diag(4, 3))))
  for (A in ops) {
    v <- verify(A)
    expect_false(v$overall == "fail", info = A$node)
  }
})

test_that("print and explain show the structure", {
  K <- linop_kron(linop(rmat(2, 2, seed = 1)), linop_perm(c(2, 1)))
  out <- capture.output(print(K))
  expect_true(any(grepl("kron", out)))
  expect_true(any(grepl("perm\\[2\\]", out)))
  df <- explain(rbind(K, cbind(linop(rmat(1, 2)), linop(rmat(1, 2)))))
  expect_true(all(c("vstack", "kron", "hstack", "perm") %in% df$node))
})
