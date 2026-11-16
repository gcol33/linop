test_that("ODE Gauss nodes recover known low orders and classical moments", {
  g <- linop:::gauss_jacobi("hermite", 3L)
  expect_equal(g$x, c(-sqrt(1.5), 0, sqrt(1.5)), tolerance = 1e-14)
  expect_equal(g$w, c(1/6, 2/3, 1/6), tolerance = 1e-14)
  g <- linop:::gauss_jacobi("laguerre", 2L)
  expect_equal(g$x, c(2-sqrt(2), 2+sqrt(2)), tolerance = 1e-14)
  expect_equal(g$w, c((2+sqrt(2))/4, (2-sqrt(2))/4), tolerance = 1e-14)
  for (kind in c("hermite", "laguerre")) for (n in c(17L, 129L, 513L, 2049L)) {
    g <- linop:::gauss_jacobi(kind, n)
    expect_true(all(diff(g$x) > 0))
    expect_true(all(g$w >= 0))
    for (k in 0:12) {
      truth <- if (kind == "laguerre") factorial(k) else
        if (k %% 2L == 0L) gamma((k+1)/2)/sqrt(pi) else 0
      expect_lt(abs(sum(g$w*g$x^k)-truth)/max(1, abs(truth)), 2e-13)
    }
  }
})

test_that("ODE quadrature recovers analytic Gaussian and Laplace integrals over twenty seeds", {
  set.seed(4519)
  a <- runif(24, .1, 2); b <- runif(24, -3, 3)
  h <- linop:::gauss_jacobi("hermite", 129L)
  l <- linop:::gauss_jacobi("laguerre", 129L)
  for (j in seq_along(a)) {
    truth <- exp(-b[j]^2/(4*(1+a[j])))/sqrt(1+a[j])
    expect_lt(abs(sum(h$w*exp(-a[j]*h$x^2)*cos(b[j]*h$x))-truth), 1e-14)
    expect_lt(abs(sum(l$w*exp(-a[j]*l$x))-1/(1+a[j])), 1e-14)
  }
})
