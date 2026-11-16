## The flag propagation table checked by brute force, including the
## A^H A versus A^T A distinction on complex inputs.
##
## The table follows base R:
## base R's crossprod(A) is A^T A, which is symmetric and generally not
## hermitian.

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

test_that("a Gram product is definite only on the side where full rank means injective", {
  ## A wide operator of full rank has rank m < n, so X^H X is singular and
  ## X X^H definite; a tall one the other way round.
  wide <- linop(rmat(3, 6, seed = 1), properties = c(full_rank = TRUE))
  tall <- linop(rmat(6, 3, seed = 2), properties = c(full_rank = TRUE))
  expect_true(is.na(capv(adjoint(wide) %*% wide, "positive_definite")))
  expect_true(isTRUE(capv(wide %*% adjoint(wide), "positive_definite")))
  expect_true(isTRUE(capv(adjoint(tall) %*% tall, "positive_definite")))
  expect_true(is.na(capv(tall %*% adjoint(tall), "positive_definite")))
  G <- as.matrix(adjoint(wide) %*% wide)
  expect_lt(min(eigen(G, symmetric = TRUE, only.values = TRUE)$values), 1e-10)
})

## ---------------------------------------------------------- structured nodes --

is_psd <- function(M) {
  is_herm(M) && min(eigen((M + Conj(t(M))) / 2, symmetric = TRUE, only.values = TRUE)$values) > -1e-8
}
is_pd <- function(M) {
  is_herm(M) && min(eigen((M + Conj(t(M))) / 2, symmetric = TRUE, only.values = TRUE)$values) > 1e-10
}
is_unitary <- function(M) isTRUE(all.equal(Conj(t(M)) %*% M, diag(ncol(M)), check.attributes = FALSE))
full_rank_of <- function(M) qr(M)$rank == min(dim(M))

## Every capability the lattice claims on a structured node, checked against its
## materialisation, over operands with a mix of established, declared and
## unknown capabilities.
expect_structured_sound <- function(A, label) {
  expect_sound(A, label)
  M <- as.matrix(A)
  if (isTRUE(capv(A, "positive_definite"))) expect_true(is_pd(M), info = paste(label, "pd"))
  if (isTRUE(capv(A, "positive_semidefinite"))) expect_true(is_psd(M), info = paste(label, "psd"))
  if (isTRUE(capv(A, "unitary"))) expect_true(is_unitary(M), info = paste(label, "unitary"))
  if (isTRUE(capv(A, "full_rank"))) expect_true(full_rank_of(M), info = paste(label, "full rank"))
  if (isTRUE(capv(A, "triangular_upper"))) expect_true(all(M[lower.tri(M)] == 0), info = paste(label, "upper"))
}

operand_pool <- function(seed) {
  set.seed(seed)
  H <- crossprod(rmat(3, 3)) + diag(3)
  Z <- zmat(3, 3); Hz <- Conj(t(Z)) %*% Z + diag(3)
  list(
    hpd = linop(H, properties = c(positive_definite = TRUE)),
    hpdz = linop(Hz, properties = c(positive_definite = TRUE)),
    sym_indef = linop(H - 3 * diag(3)),
    general = linop(rmat(3, 3)),
    complex = linop(zmat(3, 3)),
    tall = linop(rmat(4, 2), properties = c(full_rank = TRUE)),
    wide = linop(rmat(2, 4), properties = c(full_rank = TRUE)),
    upper = linop(upper.tri(diag(3), diag = TRUE) * rmat(3, 3),
                  properties = c(triangular_upper = TRUE)),
    perm = linop_perm(c(2, 3, 1)),
    involution = linop_perm(c(3, 2, 1)),
    eye = linop_eye(3),
    diag = linop_scaling(c(1, -2, 3)))
}

test_that("the structured nodes claim nothing their materialisation contradicts", {
  for (seed in 1:5) {
    pool <- operand_pool(seed)
    nm <- names(pool)
    for (i in nm) for (j in nm) {
      A <- pool[[i]]; B <- pool[[j]]
      expect_structured_sound(linop_kron(A, B), paste("kron", i, j, seed))
      expect_structured_sound(linop_blockdiag(A, B), paste("blockdiag", i, j, seed))
      if (A$dim[1L] == B$dim[1L]) expect_structured_sound(cbind(A, B), paste("cbind", i, j, seed))
      if (A$dim[2L] == B$dim[2L]) expect_structured_sound(rbind(A, B), paste("rbind", i, j, seed))
    }
    for (i in nm) {
      A <- pool[[i]]
      if (A$dim[1L] != A$dim[2L]) next
      for (k in 0:3) expect_structured_sound(linop_power(A, k), paste("power", i, k, seed))
      if (i %in% c("hpd", "hpdz", "general", "complex", "upper", "perm", "involution", "eye", "diag")) {
        expect_structured_sound(solve(A), paste("inverse", i, seed))
      }
    }
  }
})

test_that("the structured rules fire where they should", {
  pool <- operand_pool(1)
  expect_true(isTRUE(capv(linop_kron(pool$hpd, pool$hpdz), "positive_definite")))
  expect_true(isTRUE(capv(linop_blockdiag(pool$hpd, pool$eye), "positive_definite")))
  expect_true(isTRUE(capv(linop_power(pool$sym_indef, 2), "positive_semidefinite")))
  expect_true(isTRUE(capv(linop_kron(pool$perm, pool$involution), "unitary")))
  expect_true(isTRUE(capv(pool$involution, "symmetric")))
  expect_false(isTRUE(capv(pool$perm, "symmetric")))
  expect_true(isTRUE(capv(solve(pool$upper), "triangular_upper")))
  ## a tall and a wide factor of full rank give no full rank: 4 x 2 (x) 2 x 4 is
  ## 8 x 8 of rank 4
  expect_true(is.na(capv(linop_kron(pool$tall, pool$wide), "full_rank")))
  expect_false(full_rank_of(as.matrix(linop_kron(pool$tall, pool$wide))))
  expect_true(isTRUE(capv(cbind(pool$wide, linop(rmat(2, 3))), "full_rank")))
  expect_true(isTRUE(capv(rbind(pool$tall, linop(rmat(3, 2))), "full_rank")))
  expect_true(isTRUE(capv(linop_lowrank(rmat(5, 2), d = c(2, 1)), "positive_semidefinite")))
  expect_true(is.na(capv(linop_lowrank(rmat(5, 2), d = c(2, -1)), "positive_semidefinite")))
})

test_that("a bare declaration under a structured node stays visible at the top", {
  ## the laundering case through the new rules
  d1 <- linop(function(X) 2 * X, dim = c(3, 3),
              properties = c(hermitian = TRUE, positive_definite = TRUE))
  strict <- requirement(sources = c("construction", "computation", "theorem"),
                        guarantees = "identity")
  for (A in list(linop_kron(d1, linop_eye(2)), linop_blockdiag(d1, linop_eye(2)),
                 linop_power(d1, 3))) {
    ev <- linop:::cape(A, "positive_definite")
    expect_equal(ev$source, "construction")
    expect_false(evidence_satisfies(ev, strict))
  }
})
