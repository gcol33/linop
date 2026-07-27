## Shared fixtures. Every node type, real and complex where the backing storage
## admits it. Sparse is real-only: Matrix 1.7.5 declares the zMatrix virtual
## classes but not the concrete ones (see dev_notes/S0.1-dispatch.md section 7).

rmat <- function(m, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  matrix(stats::rnorm(m * n), m, n)
}

zmat <- function(m, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  matrix(complex(real = stats::rnorm(m * n), imaginary = stats::rnorm(m * n)), m, n)
}

## A shift operator as a callback leaf, with a correct adjoint.
shift_op <- function(n, dtype = "double") {
  linop(function(X) rbind(X[-1L, , drop = FALSE], 0),
        adjoint = function(X) rbind(0, X[-n, , drop = FALSE]),
        dim = c(n, n), dtype = dtype)
}

## ------------------------------------------------- closed-form ground truth --
## Section 10. Each fixture below has its spectrum, or its inverse, in closed
## form, so a solver is checked against truth rather than against itself.

## 1-D Dirichlet Laplacian, tridiag(-1, 2, -1). Eigenvalues
## 4 sin^2(k pi / (2(n+1))), exact for any n and k, so the condition number is
## known and CG's predicted rate is computable.
laplacian_1d_apply <- function(X) {
  n <- nrow(X)
  Y <- 2 * X
  Y[seq_len(n - 1L), ] <- Y[seq_len(n - 1L), , drop = FALSE] - X[2:n, , drop = FALSE]
  Y[2:n, ] <- Y[2:n, , drop = FALSE] - X[seq_len(n - 1L), , drop = FALSE]
  Y
}

laplacian_1d <- function(n) {
  linop(laplacian_1d_apply, adjoint = laplacian_1d_apply, dim = c(n, n),
        properties = c(hermitian = TRUE, positive_definite = TRUE))
}

laplacian_1d_eigenvalues <- function(n) 4 * sin(seq_len(n) * pi / (2 * (n + 1)))^2

## KMS, rho^|i-j|. Positive definite for |rho| < 1, with a tridiagonal inverse in
## closed form.
kms_matrix <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))

kms_inverse <- function(n, rho) {
  T <- diag(c(1, rep(1 + rho^2, n - 2L), 1), n)
  for (i in seq_len(n - 1L)) {
    T[i, i + 1L] <- -rho
    T[i + 1L, i] <- -rho
  }
  T / (1 - rho^2)
}

## Prescribed spectrum, Q diag(lambda) Q^H. Clustering and gaps dialled in
## deliberately; the condition number is exactly max(lambda)/min(lambda).
spd_prescribed <- function(n, lambda, seed = 1L) {
  set.seed(seed)
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n, n)))
  M <- Q %*% diag(lambda, n) %*% t(Q)
  (M + t(M)) / 2
}

## Complex hermitian positive definite, by the same route.
hpd_prescribed <- function(n, lambda, seed = 1L) {
  set.seed(seed)
  Q <- qr.Q(qr(matrix(complex(real = stats::rnorm(n * n),
                              imaginary = stats::rnorm(n * n)), n, n)))
  M <- Q %*% diag(lambda, n) %*% Conj(t(Q))
  (M + Conj(t(M))) / 2
}

as_spd_linop <- function(M) {
  linop(M, properties = c(hermitian = TRUE, positive_definite = TRUE))
}

## Every node type, as a list of (label, linop, dense reference).
all_node_fixtures <- function(seed = 42L) {
  set.seed(seed)
  out <- list()
  add <- function(label, A, M) out[[length(out) + 1L]] <<- list(label = label, A = A, M = M)

  Mr <- rmat(5, 4); add("dense real", linop(Mr), Mr)
  Mz <- zmat(4, 4); add("dense complex", linop(Mz), Mz)

  add("identity", linop_eye(4), diag(1, 4))

  dr <- c(2, -1, 3, 0.5); add("diag real", linop_scaling(dr), diag(dr))
  dz <- complex(real = c(1, 2, -1, 0.5), imaginary = c(0.5, -1, 2, 1))
  add("diag complex", linop_scaling(dz), diag(dz))

  n <- 4
  S <- diag(0, n); for (i in seq_len(n - 1L)) S[i, i + 1L] <- 1
  add("fun real", shift_op(n), S)

  fz <- function(X) (1 + 2i) * X
  gz <- function(X) Conj(1 + 2i) * X
  add("fun complex", linop(fz, adjoint = gz, dim = c(n, n), dtype = "complex"),
      diag((1 + 2i), n))

  add("transpose", t(linop(Mr)), t(Mr))
  add("adjoint", adjoint(linop(Mz)), Conj(t(Mz)))
  add("conjugate", Conj(linop(Mz)), Conj(Mz))
  add("scale real", 2.5 * linop(Mr), 2.5 * Mr)
  add("scale complex", (1 + 1i) * linop(Mz), (1 + 1i) * Mz)

  M2 <- rmat(5, 4)
  add("sum", linop(Mr) + linop(M2), Mr + M2)
  add("product", linop(Mr) %*% linop(rmat(4, 3)), NULL)  # reference filled below
  out[[length(out)]]$M <- as.matrix(out[[length(out)]]$A)

  if (requireNamespace("Matrix", quietly = TRUE)) {
    Sp <- Matrix::sparseMatrix(i = c(1, 2, 3, 4), j = c(1, 1, 3, 4),
                               x = c(1, 2, 3, 4), dims = c(4, 4))
    add("sparse real", linop(Sp), as.matrix(Sp))
  }
  out
}
