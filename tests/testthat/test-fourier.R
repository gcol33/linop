## Phase 6's first unit: operators on the circle held as banded operators on
## Fourier coefficients. The fixture is closed form twice over: for
## V = a^2 sin^2 x - a cos x the ground state of -u'' + V u is u = exp(a cos x)
## at eigenvalue 0, and the Fourier coefficients of exp(a cos x) are the modified
## Bessel values I_|k|(a).

hill <- function(a) {
  linop_space("fourier", multiplier = c(-a^2 / 4, -a / 2, a^2 / 2, -a / 2, -a^2 / 4))
}

test_that("the ground state is exp(a cos x), value and eigenvector, with a passing certificate", {
  for (a in c(0.5, 2, 5)) {
    fit <- eigs(hill(a), k = 1, which = "smallest_algebraic")
    expect_lt(abs(fit$values), 1e-10, label = sprintf("a = %g", a))
    cert <- fit$certificate
    expect_equal(cert_status(cert, "forward error"), "pass")
    expect_equal(cert_status(cert, "truncation bound"), "pass")
    expect_equal(cert$values$bound >= abs(fit$values), TRUE)
    ## the eigenvector against the Bessel coefficients, up to sign and scale
    n <- fit$n
    u <- Re(fit$vectors[, 1L])
    ## besselI(x, nu): the argument first, the order second
    ref <- besselI(a, abs(seq.int(-n, n)))
    u <- u * sign(sum(u * ref)) * sqrt(sum(ref^2)) / sqrt(sum(u^2))
    expect_lt(max(abs(u - ref)) / max(ref), 1e-9, label = sprintf("a = %g", a))
  }
})

test_that("the truncation term is the exact residual of the extended vector on a wider section", {
  ## the identity behind the certificate: ||(H - q) u~||^2 = ||(H_n - q) u||^2 +
  ## ||r_out||^2, checked by applying a much wider section to the zero-extended u
  a <- 5
  P <- hill(a)
  for (n in c(6L, 10L, 16L)) {
    F <- finite_section(P, n)
    e <- eigen(as.matrix(F), symmetric = TRUE)
    q <- min(e$values); u <- e$vectors[, which.min(e$values)]
    W <- finite_section(P, n + 10L)
    ut <- c(numeric(10), u, numeric(10))
    full <- sqrt(sum((as.matrix(W) %*% ut - q * ut)^2))
    res <- sqrt(sum((as.matrix(F) %*% u - q * u)^2))
    edge <- linop:::fourier_edge(P, matrix(u, ncol = 1L))
    expect_equal(full, sqrt(res^2 + edge^2), tolerance = 1e-10)
  }
})

test_that("the section agrees with its materialisation in every mode, complex multiplier too", {
  ## V = cos x + sin 2x is real and odd in part, so the coefficient matrix is
  ## complex hermitian
  vhat <- c(1i / 2, 1 / 2, 0, 1 / 2, -1i / 2)
  P <- linop_space("fourier", multiplier = vhat)
  expect_equal(P$dtype, "complex")
  expect_true(isTRUE(capv(P, "hermitian")))
  F <- finite_section(P, 7)
  M <- as.matrix(F)
  expect_lt(max(Mod(M - Conj(t(M)))), 1e-15)
  X <- zmat(15, 2, seed = 1)
  for (mode in c("N", "T", "C", "R")) {
    ref <- switch(mode, N = M %*% X, T = t(M) %*% X, C = Conj(t(M)) %*% X, R = Conj(M) %*% X)
    expect_lt(max(Mod(linop:::linop_apply(F, X, mode) - ref)), 1e-13, label = mode)
  }
})

test_that("the free operator has the spectrum k^2, through Lanczos by construction", {
  ## Every nonzero eigenvalue of -d^2/dx^2 on the circle is double (e^{ikx} and
  ## e^{-ikx}), and a Krylov space from one start vector can return one copy:
  ## the same three smallest come back as 0, 1, 1 or 0, 1, 4 depending on k. So
  ## what is asserted is that each value is some k^2 and the smallest is 0.
  fit <- eigs(linop_space("fourier"), k = 3, which = "smallest_algebraic")
  v <- sort(fit$values)
  expect_lt(abs(v[1L]), 1e-10)
  expect_lt(max(abs(v - round(sqrt(pmax(v, 0)))^2)), 1e-9)
  expect_match(fit$method, "lanczos")
})

test_that("a non-self-adjoint operator gets no distance bound", {
  ## d/dx has symbol i k: not real, so not self-adjoint
  D <- linop_space("fourier", derivative = c(0, 1))
  expect_false(isTRUE(capv(D, "hermitian")))
  fit <- eigs(finite_section(D, 8), k = 1, which = "smallest")
  expect_lt(Mod(fit$values), 1e-12)
  expect_equal(cert_status(fit$certificate, "forward error"), "not_checked")
  expect_match(fit$certificate$checks$detail[fit$certificate$checks$check == "forward error"],
               "pseudospectrum")
})

test_that("the width is chosen, and recorded", {
  fit <- eigs(hill(5), k = 1, which = "smallest_algebraic")
  expect_gt(length(fit$widths), 1L)
  expect_match(fit$certificate$dispatch$reason, "width 16 -> ")
  expect_equal(fit$n, fit$widths[length(fit$widths)])
})

test_that("what linop_space() and its section refuse", {
  expect_error(linop_space("legendre"), "kind is")
  expect_error(linop_space("fourier", multiplier = c(1, 2)), "odd number")
  expect_error(finite_section(hill(1), 1), "narrower than the multiplier")
  expect_error(eigs(hill(1), k = 1, v0 = rep(1, 3)), "does not take v0")
  expect_error(hill(1) %*% 1, "no block to act on|finite dimensions")
  expect_error(finite_section(linop(diag(2)), 3), "no truncation rule")
})

test_that("print shows the operator under its section", {
  out <- capture.output(print(finite_section(hill(1), 5)))
  expect_true(any(grepl("section \\|k\\| <= 5", out)))
  expect_true(any(grepl("fourier operator", out)))
})
