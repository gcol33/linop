## Gate 1: the flag propagation table checked by brute force, including the
## A^H A versus A^T A distinction on complex inputs.
##
## The table is the corrected one from dev_notes/S0.1-dispatch.md section 6:
## base R's crossprod(A) is A^T A, which is symmetric and generally not
## hermitian. The plan's original table had these swapped.

is_herm <- function(M) isTRUE(all.equal(M, Conj(t(M)), check.attributes = FALSE))
is_symm <- function(M) isTRUE(all.equal(M, t(M), check.attributes = FALSE))

## A declared flag must never be TRUE when the materialised matrix contradicts
## it. NA is always allowed: the lattice is conservative, not complete.
expect_sound <- function(A, label) {
  M <- as.matrix(A)
  if (isTRUE(capv(A, "hermitian"))) expect_true(is_herm(M), info = paste(label, "hermitian"))
  if (isTRUE(capv(A, "symmetric"))) expect_true(is_symm(M), info = paste(label, "symmetric"))
  if (isTRUE(capv(A, "real"))) expect_false(is.complex(M) && any(Im(M) != 0),
                                            info = paste(label, "real"))
  if (isTRUE(capv(A, "diagonal"))) {
    off <- M; diag(off) <- 0
    expect_true(all(off == 0), info = paste(label, "diagonal"))
  }
  if (isTRUE(capv(A, "positive_semidefinite"))) {
    ev <- eigen((M + Conj(t(M))) / 2, only.values = TRUE)$values
    expect_gt(min(Re(ev)), -1e-8)
  }
}

test_that("crossprod(A) is A^T A: symmetric always, hermitian only when A is real", {
  Z <- zmat(4, 3, seed = 1)
  A <- linop(Z)
  cp <- crossprod(A)

  ## base R agreement is the binding constraint
  expect_equal(as.matrix(cp), crossprod(Z), tolerance = 1e-12, ignore_attr = TRUE)
  expect_true(is_symm(crossprod(Z)))
  expect_false(is_herm(crossprod(Z)))

  expect_true(isTRUE(capv(cp, "symmetric")))
  expect_false(isTRUE(capv(cp, "hermitian")))

  ## and for a real operator the same expression is both
  R <- rmat(4, 3, seed = 2)
  cpr <- crossprod(linop(R))
  expect_true(isTRUE(capv(cpr, "symmetric")))
  expect_true(isTRUE(capv(cpr, "hermitian")))
})

test_that("adjoint(A) %*% A is A^H A: hermitian and PSD by construction", {
  Z <- zmat(4, 3, seed = 3)
  A <- linop(Z)
  g <- adjoint(A) %*% A

  expect_equal(as.matrix(g), Conj(t(Z)) %*% Z, tolerance = 1e-12, ignore_attr = TRUE)
  expect_true(is_herm(Conj(t(Z)) %*% Z))
  expect_true(isTRUE(capv(g, "hermitian")))
  expect_true(isTRUE(capv(g, "positive_semidefinite")))

  ## the evidence is unconditional: it does not depend on anything claimed of A
  ev <- linop:::cape(g, "hermitian")
  expect_equal(ev$source, "construction")
  expect_length(ev$depends_on, 0)

  ## A %*% adjoint(A) likewise
  g2 <- A %*% adjoint(A)
  expect_true(isTRUE(capv(g2, "hermitian")))
})

test_that("scalar multiplication respects the reality of the scalar", {
  Z <- zmat(4, 4, seed = 4)
  H <- linop_dense((Z + Conj(t(Z))) / 2)     # hermitian by computation
  expect_true(isTRUE(capv(H, "hermitian")))

  expect_true(isTRUE(capv(2 * H, "hermitian")))         # real scalar keeps it
  expect_false(isTRUE(capv((1 + 1i) * H, "hermitian"))) # complex scalar does not
  expect_true(isTRUE(capv((1 + 1i) * H, "symmetric")) ||
              is.na(capv((1 + 1i) * H, "symmetric")))
})

test_that("transpose and adjoint move flags differently for complex operators", {
  Z <- zmat(4, 4, seed = 5)
  H <- linop_dense((Z + Conj(t(Z))) / 2)
  S <- linop_dense((Z + t(Z)) / 2)

  ## adjoint preserves hermitian; transpose does not, unless real
  expect_true(isTRUE(capv(adjoint(H), "hermitian")))
  expect_false(isTRUE(capv(t(H), "hermitian")))
  ## transpose preserves symmetric; adjoint does not, unless real
  expect_true(isTRUE(capv(t(S), "symmetric")))
  expect_false(isTRUE(capv(adjoint(S), "symmetric")))
  ## Conj preserves both
  expect_true(isTRUE(capv(Conj(H), "hermitian")))
  expect_true(isTRUE(capv(Conj(S), "symmetric")))
})

test_that("sum requires the property of every term", {
  H1 <- linop_dense(diag(c(1, 2, 3)))
  H2 <- linop_dense(diag(c(4, 5, 6)))
  G  <- linop(rmat(3, 3, seed = 6))
  expect_true(isTRUE(capv(H1 + H2, "symmetric")))
  expect_false(isTRUE(capv(H1 + G, "symmetric")))
})

test_that("brute force: no declared flag contradicts the materialised matrix", {
  set.seed(99)
  for (i in seq_len(300)) {
    cplx <- runif(1) < 0.5
    n <- sample(2:4, 1)
    mk <- function() if (cplx) zmat(n, n) else rmat(n, n)
    Z <- mk()
    base <- list(
      linop(Z),
      linop_dense((Z + Conj(t(Z))) / 2),
      linop_dense((Z + t(Z)) / 2),
      linop_eye(n),
      linop_scaling(if (cplx) complex(real = rnorm(n), imaginary = rnorm(n)) else rnorm(n)))
    A <- base[[sample(seq_along(base), 1)]]
    B <- base[[sample(seq_along(base), 1)]]
    s <- if (cplx && runif(1) < 0.5) complex(real = rnorm(1), imaginary = rnorm(1)) else rnorm(1)
    exprs <- list(
      A + B, s * A, A %*% B, t(A), adjoint(A), Conj(A),
      crossprod(A), tcrossprod(A), adjoint(A) %*% A, A %*% adjoint(A),
      t(A) %*% A, crossprod(A) + s * linop_eye(n))
    for (j in seq_along(exprs)) expect_sound(exprs[[j]], sprintf("iter %d expr %d", i, j))
  }
})
