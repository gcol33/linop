## Phase 6 step 2 and 3: operators on [-1, 1] in the ultraspherical method.
## Fixtures are closed form: sin(pi x) has Chebyshev coefficients
## 2 (-1)^((k-1)/2) J_k(pi) on odd k, exp(x) has I_0(1) and 2 I_k(1), and the
## Dirichlet Laplacian has eigenvalues (k pi / 2)^2.

cheb_eval <- function(u, x) vapply(x, function(t) sum(u * cos((seq_along(u) - 1) * acos(t))), numeric(1))
xs <- seq(-1, 1, length.out = 41)
dirichlet <- function(l = 0, r = 0) list(at = c(-1, 1), derivative = c(0, 0), value = c(l, r))
sin_pi_coef <- function(K) {
  k <- 0:K
  ifelse(k %% 2 == 1, 2 * (-1)^((k - 1) / 2) * besselJ(pi, k), 0)
}
exp_coef <- function(K) c(besselI(1, 0), 2 * besselI(1, seq_len(K)))

test_that("u'' = f with Dirichlet conditions recovers sin(pi x)", {
  L <- linop_space("chebyshev", derivative = c(0, 0, 1), boundary = dirichlet())
  u <- solve(L, -pi^2 * sin_pi_coef(40))
  expect_lt(max(abs(cheb_eval(u, xs) - sin(pi * xs))), 1e-13)
  ## the coefficients themselves, not only the values
  ref <- sin_pi_coef(length(u) - 1L)
  expect_lt(max(abs(u - ref)), 1e-13)
  cert <- attr(u, "certificate")
  expect_equal(cert_status(cert, "backward error"), "pass")
  expect_equal(cert_status(cert, "truncation residual"), "pass")
  expect_equal(cert_status(cert, "forward error"), "not_checked")
})

test_that("a variable coefficient and inhomogeneous boundary values: u'' + x u = (1 + x) e^x", {
  e <- exp_coef(30)
  f <- c(e, 0) + linop:::cheb_mult(c(0, 1), e)
  L <- linop_space("chebyshev", derivative = c(0, 0, 1), multiplier = c(0, 1),
                   boundary = dirichlet(exp(-1), exp(1)))
  u <- solve(L, f)
  expect_lt(max(abs(cheb_eval(u, xs) - exp(xs))), 1e-13)
})

test_that("Neumann conditions: u'' - u = 0 with u'(+-1) = e^(+-1) is exp(x)", {
  L <- linop_space("chebyshev", derivative = c(-1, 0, 1),
                   boundary = list(at = c(-1, 1), derivative = c(1, 1), value = c(exp(-1), exp(1))))
  u <- solve(L, 0)
  expect_lt(max(abs(cheb_eval(u, xs) - exp(xs))), 1e-13)
})

test_that("the dropped rows are the exact residual on the untruncated system", {
  ## L applied to a finitely supported vector has finite support, so applying it
  ## to the truncated solution padded with zeros gives every row there is
  L <- linop_space("chebyshev", derivative = c(0, 0, 1), multiplier = c(0.5, 0, 1),
                   boundary = dirichlet())
  for (N in c(12L, 20L)) {
    S <- linop:::cheb_system(L, N)
    f <- sin_pi_coef(30)
    F <- linop:::cheb_convert_chain(f, 0L, 2L)
    c <- solve(S, c(0, 0, F[seq_len(N - 2L)]))
    m <- linop:::cheb_measure(L, c, F, c(0, 0))
    big <- linop:::cheb_apply_L(L, c(c, numeric(30)))
    full <- sqrt(sum((linop:::pad_to(big, 80) - linop:::pad_to(F, 80))[-seq_len(N - 2L)]^2))
    expect_equal(m$truncation, full, tolerance = 1e-12)
  }
})

test_that("the Dirichlet Laplacian's smallest eigenvalues are (k pi / 2)^2", {
  E <- linop_space("chebyshev", derivative = c(0, 0, -1), boundary = dirichlet())
  fit <- eigs(E, k = 4, which = "smallest")
  expect_equal(fit$values, ((1:4) * pi / 2)^2, tolerance = 1e-12)
  expect_equal(cert_status(fit$certificate, "truncation residual"), "pass")
  expect_equal(cert_status(fit$certificate, "forward error"), "not_checked")
  ## the eigenvector of (pi/2)^2 is cos(pi x / 2), up to sign and scale
  u <- Re(fit$vectors[, 1L])
  v <- cheb_eval(u, xs); w <- cos(pi * xs / 2)
  v <- v * sum(v * w) / sum(v * v)
  expect_lt(max(abs(v - w)), 1e-10)
})

test_that("the operator records its domain and range spaces", {
  L <- linop_space("chebyshev", derivative = c(0, 0, 1), boundary = dirichlet())
  expect_equal(L$args$domain, "chebyshev")
  expect_equal(L$args$range, "ultraspherical(2)")
  expect_true(any(grepl("chebyshev -> ultraspherical\\(2\\)", capture.output(print(L)))))
})

test_that("what the Chebyshev operator refuses", {
  ## a term without the problem's conditions is allowed; solving it is not
  expect_error(solve(linop_space("chebyshev", derivative = 1, boundary = dirichlet()), 1),
               "order 0")
  expect_error(solve(linop_space("chebyshev", derivative = c(0, 0, 1),
                                 boundary = list(at = 1, derivative = 0, value = 0)), 1),
               "needs 2 boundary conditions; it has 1")
  expect_error(linop_space("chebyshev", derivative = c(0, 0, 1),
                           boundary = list(at = c(0, 1), derivative = c(0, 0), value = c(0, 0))),
               "sit at -1 or 1")
  expect_error(solve(linop_space("chebyshev", derivative = c(0, 0, 1)), 1),
               "needs 2 boundary conditions; it has 0")
  expect_error(linop_space("chebyshev", derivative = c(0, 0, 1), boundary = "x"),
               "boundary is a list")
  E <- linop_space("chebyshev", derivative = c(0, 0, -1), boundary = dirichlet())
  expect_error(eigs(E, k = 2), "smallest magnitude")
  expect_error(E %*% 1, "no block to act on|finite dimensions")
  expect_error(linop_space("legendre"), "kind is")
})
