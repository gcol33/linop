fc_eval <- function(c, x) vapply(x, function(t)
  sum(c * cos((seq_along(c) - 1L) * acos(t))), numeric(1))
fc_bc <- function(value = c(0, 0), derivative = c(0, 0))
  list(at = c(-1, 1), derivative = derivative, value = value)
fc_exp <- function(n) c(besselI(1, 0), 2 * besselI(1, seq_len(n)))
fc_sin <- function(n) {
  k <- 0:n
  ifelse(k %% 2 == 1, 2 * (-1)^((k - 1) / 2) * besselJ(pi, k), 0)
}

test_that("standardChop reproduces the published executable examples", {
  a <- 10^-(1:50); noise <- cos((1:50)^2)
  # Chebfun standardChop.m, examples, lines 30--35, read 2026-11-16.
  expected <- c(18L, 15L, 13L, 50L)
  for (i in seq_along(expected)) {
    s <- c(0, 1e-16, 1e-13, 1e-10)[i]
    expect_identical(linop:::standard_chop(a + s * noise), expected[i])
  }
  expect_identical(linop:::standard_chop(a + 1e-10 * noise, 1e-10), 10L)
  expect_identical(linop:::standard_chop(rep(1, 50)), 50L)
  expect_identical(linop:::standard_chop(numeric(17)), 1L)
  expect_identical(linop:::standard_chop(a[1:16]), 16L)
  expect_identical(linop:::standard_chop(1i * a), 18L)
  expect_identical(linop:::standard_chop(2^500 * a), 18L)
})

test_that("Chebyshev lengths agree with sourced Chebfun examples", {
  # Aurentz and Trefethen, ChebExplain, section 1: degree 14 for scaled exp(x).
  z <- linop:::function_coefficients(exp, "chebyshev", "exp")
  expect_length(z$coefficients, 15L)
  expect_identical(z$grids, c(17L, 33L))
  expect_lt(max(abs(z$coefficients - fc_exp(14))), 3e-16)
  # Aurentz and Trefethen (2017), section 3: length 51 at tolerance 1e-6.
  f <- function(x) 3 * exp(-1 / (x + 1)) - (x + 1)
  z <- linop:::function_coefficients(f, "chebyshev", "paper example", tol = 1e-6)
  expect_length(z$coefficients, 51L)
  x <- seq(-1, 1, length.out = 1001)
  expect_lt(max(abs(fc_eval(z$coefficients, x) - f(x))), 2e-6)
  # The paper, figure 6: cutoff 75 for log(1.1-x).
  z <- linop:::function_coefficients(function(x) log(1.1 - x), "chebyshev", "log")
  expect_length(z$coefficients, 75L)
  expect_lt(max(abs(fc_eval(z$coefficients, x) - log(1.1-x))), 2e-14)
})

test_that("limited smoothness and Runge functions recover independently known values", {
  x <- seq(-1, 1, length.out = 1001)
  fs <- list(function(x) 1 / (1 + 25 * x^2), function(x) abs(x)^3)
  for (i in seq_along(fs)) {
    f <- fs[[i]]
    z <- linop:::function_coefficients(f, "chebyshev", "multiplier", nmax = 65536L)
    h <- (length(z$coefficients) - 1L) %/% 2L
    # The |x|^3 coefficients are 24*(-1)^k/[pi*(4*k^2-9)*(4*k^2-1)].
    # Their absolute tail is bounded by 8/(pi*h^3) for h >= 2.
    bound <- if (i == 1L) 2e-14 else 8 / (pi * h^3)
    expect_lt(max(abs(fc_eval(z$coefficients, x) - f(x))), bound)
    expect_lt(z$tail, 1e-11)
    expect_lte(z$calls, 14L)
  }
})

test_that("function forcing recovers the closed-form Dirichlet and Neumann BVPs", {
  x <- seq(-1, 1, length.out = 41)
  problems <- list(
    list(d = c(0, 0, 1), v = 0, bc = fc_bc(), cf = -pi^2 * fc_sin(40),
         f = function(x) -pi^2 * sin(pi*x), truth = function(x) sin(pi*x)),
    list(d = c(0, 0, 1), v = c(0, 1), bc = fc_bc(exp(c(-1, 1))),
         cf = linop:::add_padded(fc_exp(30), linop:::cheb_mult(c(0, 1), fc_exp(30))),
         f = function(x) (1+x)*exp(x), truth = exp),
    list(d = c(-1, 0, 1), v = 0, bc = fc_bc(exp(c(-1, 1)), c(1, 1)), cf = 0,
         f = function(x) 0*x, truth = exp))
  for (p in problems) {
    A <- linop_space("chebyshev", p$d, p$v, p$bc)
    exact <- solve(A, p$cf); adaptive <- solve(A, p$f)
    expect_lt(max(abs(fc_eval(adaptive, x) - p$truth(x))), 1e-13)
    expect_lt(max(abs(fc_eval(exact, x) - fc_eval(adaptive, x))), 1e-13)
    expect_equal(linop:::cert_status(certificate(adaptive), "backward error"), "pass")
  }
})

test_that("variable terms and boundary functions recover sin and exp", {
  terms <- list(function(x) 0.5+0.2*x, function(x) -x, function(x) -1-x^2/2)
  exact <- linop_space("chebyshev", list(c(.5, .2), c(0, -1), c(-1.25, 0, -.25)), boundary = fc_bc())
  adaptive <- linop_space("chebyshev", terms, boundary = fc_bc())
  f <- function(x) (.5+.2*x)*sin(pi*x) - pi*x*cos(pi*x) + (1+x^2/2)*pi^2*sin(pi*x)
  x <- seq(-1, 1, length.out = 41)
  for (A in list(exact, adaptive)) {
    u <- solve(A, f)
    expect_lt(max(abs(fc_eval(u, x) - sin(pi*x))), 1e-11)
  }
  e <- linop_space("chebyshev", c(0, 0, 1), multiplier = function(x) x,
                    boundary = fc_bc(exp))
  u <- solve(e, function(x) (1+x)*exp(x))
  expect_lt(max(abs(fc_eval(u, x) - exp(x))), 1e-13)
  expect_length(e$args$function_coefficients, 2L)
  expect_true(any(grepl("boundary\\$value", certificate(u)$checks$check)))
})

test_that("Hill functions recover the zero ground state and Bessel eigenvector", {
  for (a in c(.5, 2, 5)) {
    A <- linop_space("fourier", multiplier = function(x) a^2*sin(x)^2-a*cos(x))
    B <- linop_space("fourier", multiplier = c(-a^2/4, -a/2, a^2/2, -a/2, -a^2/4))
    f <- eigs(A, k = 1, which = "smallest_algebraic")
    g <- eigs(B, k = 1, which = "smallest_algebraic")
    expect_lt(abs(f$values), 1e-10)
    expect_equal(f$values, g$values, tolerance = 1e-10)
    ref <- besselI(a, abs(-f$n:f$n))
    v <- f$vectors[, 1L]; v <- v * sum(Conj(v)*ref) / sum(Mod(v)^2)
    expect_lt(max(Mod(v-ref))/max(ref), 1e-9)
    ev <- f$certificate$evidence[["forward error"]]
    expect_false(evidence_satisfies(ev, requirement(guarantees = "deterministic_bound")))
    expect_true(any(vapply(ev$depends_on, function(e) e$guarantee == "estimate", logical(1))))
  }
})

test_that("Gauss transforms recover scaled and shifted oscillators on both domains", {
  for (w in c(.5, 1, 2)) {
    for (b in c(0, 1)) {
      A <- linop_space("hermite", multiplier = function(x) w^2*x^2+b*x)
      B <- linop_space("hermite", multiplier = c(0, b, w^2))
      f <- eigs(A, k = 3, which = "smallest_algebraic")
      g <- eigs(B, k = 3, which = "smallest_algebraic")
      expect_equal(f$values, w*(2*(0:2)+1)-b^2/(4*w^2), tolerance = 1e-10)
      expect_equal(f$values, g$values, tolerance = 1e-10)
    }
  }
  for (b in c(0, .75, 2)) {
    A <- linop_space("laguerre", derivative = 1, multiplier = function(x) b*x)
    B <- linop_space("laguerre", derivative = 1, multiplier = c(0, b))
    f <- eigs(A, k = 3, which = "smallest_algebraic")
    expect_equal(f$values, (2*(0:2)+1)*sqrt(.25+b), tolerance = 1e-10)
    expect_equal(f$values, eigs(B, k = 3, which = "smallest_algebraic")$values, tolerance = 1e-10)
  }
})

test_that("Gauss refinement recovers analytic weighted coefficients", {
  for (kind in c("hermite", "laguerre")) {
    f <- if (kind == "hermite") function(x) exp(-x^2) else function(x) exp(-x)
    z <- linop:::function_coefficients(f, kind, "multiplier", nmax = 257)
    k <- seq_along(z$coefficients) - 1L
    truth <- if (kind == "hermite") ifelse(k %% 2 == 0,
      (-1)^(k/2) * exp(.5*lfactorial(k)-lfactorial(k/2)-(k/2)*log(4)-.5*log(2)), 0)
      else 2^(-k-1)
    expect_lt(max(abs(z$coefficients-truth)), 5e-15)
    expect_gt(length(z$grids), 1L)
    expect_identical(z$calls, length(z$grids)+1L)
    expect_lte(z$calls, 6L)
  }
})

test_that("independent samples reject aliases and random noise", {
  # Aurentz and Trefethen (2017), section 2: T_128 has length 129.
  z <- linop:::function_coefficients(function(x) cos(128*acos(x)), "chebyshev", "alias")
  expect_length(z$coefficients, 129L)
  expect_lt(max(abs(z$coefficients-c(numeric(128), 1))), 2e-14)
  set.seed(632)
  expect_error(linop_space("chebyshev", multiplier = function(x) rnorm(length(x)), nmax = 65),
               "multiplier function unresolved.*nmax = 65")
})

test_that("construction refuses unresolved or malformed functions by argument name", {
  expect_error(linop_space("chebyshev", multiplier = sign, nmax = 33), "multiplier function unresolved.*nmax")
  expect_error(linop_space("fourier", multiplier = function(x) rep(NA_real_, length(x))), "multiplier function.*non-finite")
  expect_error(linop_space("chebyshev", multiplier = function(x) numeric(3)), "multiplier function.*length")
  expect_error(linop_space("chebyshev", list(0, function(x) stop("deliberate"))), "derivative\\[\\[2\\]\\] function failed.*deliberate")
  A <- linop_space("chebyshev", c(0, 0, 1), boundary = fc_bc())
  expect_error(solve(A, sign, control = list(nmax = 33)), "b function unresolved")
  expect_error(solve(A, function(x) x/0), "b function.*non-finite")
  expect_error(linop_space("hermite", multiplier = function(x) rep(Inf, length(x))), "multiplier function.*non-finite")
})

test_that("one detected scalar fallback works and vectorized call counts stay small", {
  calls <- 0L
  f <- function(x) { calls <<- calls+1L; exp(x) }
  A <- linop_space("chebyshev", multiplier = f)
  expect_identical(calls, A$args$function_coefficients[[1L]]$calls)
  expect_identical(calls, 3L)
  vector_errors <- 0L
  scalar <- function(x) {
    if (length(x) > 1L) vector_errors <<- vector_errors+1L
    if (x > 0) exp(x) else exp(x)
  }
  B <- linop_space("chebyshev", multiplier = scalar)
  expect_identical(vector_errors, 1L)
  expect_equal(A$args$terms, B$args$terms, tolerance = 1e-15)
  expect_true(B$args$function_coefficients[[1L]]$scalar)
})

test_that("arithmetic and adjoints preserve approximations and their evidence", {
  for (kind in c("chebyshev", "fourier", "hermite", "laguerre")) {
    f <- switch(kind, chebyshev = function(x) x^2, fourier = cos,
                hermite = function(x) x^2, laguerre = function(x) x)
    v <- switch(kind, chebyshev = c(.5, 0, .5), fourier = c(.5, 0, .5),
                hermite = c(0, 0, 1), laguerre = c(0, 1))
    d <- if (kind == "laguerre") 1 else c(0, 0, -1)
    A <- linop_space(kind, d, f, if (kind == "chebyshev") fc_bc() else NULL)
    B <- linop_space(kind, d, v, if (kind == "chebyshev") fc_bc() else NULL)
    C <- 2*A+B
    expect_identical(C$node, A$node)
    expect_length(C$args$function_coefficients, 1L)
    expect_length(adjoint(A)$args$function_coefficients, 1L)
    if (kind == "chebyshev") expect_equal(C$args$terms, (3*B)$args$terms, tolerance = 1e-14)
    else {
      expect_lt(max(Mod(as.matrix(finite_section(C, 12)) - as.matrix(finite_section(3*B, 12)))), 1e-12)
      cert <- verify(finite_section(C, 12))
      expect_true(any(cert$checks$guarantee == "estimate"))
    }
    expect_output(print(A), "length.*tail/vscale")
    expect_output(explain(A), "length.*tail/vscale")
  }
})

test_that("chop estimates qualify all dependent certificate rows", {
  A <- linop_space("chebyshev", c(0, 0, 1), boundary = fc_bc())
  u <- solve(A, function(x) -pi^2*sin(pi*x))
  cert <- certificate(u)
  i <- grepl("function coefficients", cert$checks$check)
  expect_identical(cert$checks$source[i], "computation")
  expect_identical(cert$checks$guarantee[i], "estimate")
  expect_lt(cert$values$function_coefficients[[1L]]$tail, 1e-14)
  expect_false(evidence_satisfies(cert$evidence[["truncation residual"]],
                                  requirement(guarantees = "deterministic_bound")))
  exact <- certificate(solve(A, -pi^2*fc_sin(40)))
  expect_identical(linop:::function_certificate(exact, list()), exact)
  expect_false(any(grepl("function coefficients", exact$checks$check)))
  expect_true(evidence_satisfies(exact$evidence[["truncation residual"]] %||%
    evidence("theorem", "deterministic_bound", 1), requirement(guarantees = "deterministic_bound")))
})

test_that("exact-coefficient certificates retain their serialized contract", {
  A <- linop_space("chebyshev", c(0, 0, 1), boundary = fc_bc())
  expected <- dget(test_path("fixtures", "function-coefficients-zero.R"))
  actual <- certificate(solve(A, 0))
  expect_identical(serialize(actual, NULL), serialize(expected, NULL))
})

test_that("function-built interval spectra carry evidence in both discretizations", {
  A <- linop_space("chebyshev", derivative = list(0, 0, function(x) -1+0*x), boundary = fc_bc())
  for (method in c("auto", "galerkin")) {
    f <- eigs(A, k = 2, which = "smallest", method = method)
    expect_equal(Re(f$values), ((1:2)*pi/2)^2, tolerance = 1e-12)
    expect_true(any(f$certificate$checks$guarantee == "estimate"))
    row <- if (method == "galerkin") "upper bound" else "truncation residual"
    expect_false(evidence_satisfies(f$certificate$evidence[[row]],
                                    requirement(guarantees = "deterministic_bound")))
  }
})
