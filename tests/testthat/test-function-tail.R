test_that("sampled analytic bounds cover closed-form functions and seeded analytic families", {
  for (kind in c("chebyshev", "fourier")) {
    fs <- analytic_function_fixtures(kind)
    spec <- linop:::coefficient_registry()[[kind]]
    x <- if (kind == "chebyshev") seq(-1, 1, length.out = 4097) else 2*pi*(0:4096)/4097
    for (i in seq_along(fs)) {
      f <- fs[[i]]
      for (tol in c(.Machine$double.eps, 1e-8)) {
        z <- linop:::function_coefficients(f, kind, "multiplier", tol = tol)
        b <- linop:::coefficient_tail_bound(z)
        error <- max(Mod(spec$eval(z$coefficients, x)-f(x)))
        upper <- analytic_error_upper(kind, names(fs)[i], f, z$coefficients)
        expect_true(b$available, info = sprintf("%s fixture %d tol %g: %s", kind, i, tol, b$reason))
        expect_true(b$bound >= error, info = sprintf("%s fixture %d tol %g", kind, i, tol))
        expect_true(b$bound >= upper, info = sprintf("global analytic upper bound: %s fixture %d tol %g", kind, i, tol))
      }
    }
  }
})

test_that("complex refusal and nonanalytic extensions retain honest evidence", {
  for (f in list(function(x) { if (is.complex(x)) stop("real only"); exp(x) }, abs)) {
    A <- linop_space("chebyshev", multiplier = f, tol = 1e-6)
    cert <- linop:::function_certificate(
      linop:::build_certificate(data.frame(check = "data", status = "pass", source = "computation",
        guarantee = "identity", confidence = 1, detail = "data"), "operator"),
      linop:::function_records(A))
    i <- grepl("function tail", cert$checks$check)
    expect_identical(cert$checks$source[i], "computation")
    expect_identical(cert$checks$guarantee[i], "estimate")
    expect_match(cert$checks$detail[i], "complex|analytic|contour samples")
  }
})

test_that("nearby analytic poles receive a quantitative contour estimate", {
  r <- .995
  f <- function(x) (1-r^2)/(1-2*r*cos(x)+r^2)
  attr(f, "rates") <- r
  attr(f, "weights") <- 1
  z <- linop:::function_coefficients(f, "fourier", "multiplier", tol = 1e-8)
  bound <- linop:::coefficient_tail_bound(z)
  expect_true(bound$available)
  expect_gt(bound$samples, 16384)
  expect_gte(bound$bound, analytic_error_upper("fourier", "seed", f, z$coefficients))
})

test_that("Weyl includes multiplier bounds and scaling in the spectral evidence", {
  A <- linop_space("fourier", multiplier = function(x) exp(cos(x)), tol = 1e-7)
  f <- eigs(A, k = 1, which = "smallest_algebraic", tol = 1e-5)
  cert <- certificate(f)
  delta <- cert$values$function_perturbation
  expect_gt(delta, 0)
  expect_true("function eigenvalue perturbation" %in% cert$checks$check)
  expect_false(evidence_satisfies(cert$evidence[["function eigenvalue perturbation"]],
    requirement(guarantees = "deterministic_bound")))
  expect_identical(cert$evidence[["function eigenvalue perturbation"]]$depends_on[[1L]]$source, "theorem")
  B <- linop_space("fourier", multiplier = function(x) exp(cos(x)))
  ref <- eigs(B, k = 1, which = "smallest_algebraic", tol = 1e-11)
  expect_lte(abs(f$values-ref$values), max(cert$values$function_forward_bound))
  g <- eigs(3*A, k = 1, which = "smallest_algebraic", tol = 1e-5)
  expect_equal(certificate(g)$values$function_perturbation, 3*delta, tolerance = 1e-14)
  expect_length(cert$values$function_coefficients, 1L)
})

test_that("request tolerances reduce width within a measured error allocation", {
  A <- linop_space("fourier", multiplier = function(x) exp(cos(x)))
  for (tol in c(1e-6, 1e-8, 1e-10, 1e-12)) {
    B <- linop:::function_request_operator(A, tol)
    expect_lt(B$args$d, A$args$d)
    z <- B$args$function_coefficients[[1L]]
    bound <- linop:::coefficient_tail_bound(z)$bound
    expect_lte(bound, tol*(A$args$vnorm+1)/10)
    fit <- eigs(A, k = 1, which = "smallest_algebraic", tol = tol)
    expect_false(any(certificate(fit)$checks$status == "fail"))
  }
  f <- function(x) { if (is.complex(x)) stop("real only"); exp(x) }
  A <- linop_space("chebyshev", multiplier = f)
  expect_identical(linop:::function_request_operator(A, 1e-6), A)
})

test_that("Fourier iteration uses the function-space certificate target", {
  potential <- function(x) {
    r <- .65; D <- 1-2*r*cos(x)+r^2
    -r*cos(x)/D+3*r^2*sin(x)^2/D^2
  }
  A <- linop_space("fourier", multiplier = potential)
  for (tol in c(1e-6, 1e-8, 1e-10)) {
    fit <- eigs(A, k = 1, which = "smallest_algebraic", tol = tol,
                seed = 237, maxit = 5000, ncv = 32)
    cert <- certificate(fit)
    expect_false(any(cert$checks$status == "fail"))
    expect_lte(Mod(fit$values), max(cert$values$function_forward_bound))
    expect_lt(length(cert$values$function_coefficients[[1L]]$coefficients),
              length(A$args$vhat))
  }
  for (tol in c(1e-12, 1e-14)) expect_identical(linop:::function_request_operator(A, tol), A)
})

test_that("request-sized interval multipliers recover a closed-form boundary value problem", {
  bc <- list(at = c(-1, 1), derivative = c(0, 0), value = c(0, 0))
  A <- linop_space("chebyshev", c(0, 0, -1), exp, bc)
  x <- seq(-1, 1, length.out = 101)
  for (tol in 10^-(seq(6, 14, 2))) {
    u <- solve(A, function(x) (pi^2+exp(x))*sin(pi*x), tol = tol)
    err <- max(Mod(linop:::coefficient_registry()$chebyshev$eval(u, x)-sin(pi*x)))
    expect_lt(err, max(tol, 1e-14))
    expect_false(any(certificate(u)$checks$status == "fail"))
  }
})

test_that("nested transforms recover frequencies and reuse only matching points", {
  for (n in c(32L, 64L)) {
    x <- 2*pi*(0:(n-1L))/n
    c <- linop:::fourier_transform(cos((n/2)*x), NULL)
    expect_equal(c[c(1L, length(c))], c(.5+0i, .5+0i), tolerance = 1e-13)
    expect_equal(sum(Mod(c[-c(1L, length(c))])), 0, tolerance = 1e-13)
  }
  for (kind in c("chebyshev", "fourier")) {
    seen <- numeric()
    f <- function(x) { seen <<- c(seen, x); if (kind == "chebyshev") exp(x) else exp(cos(x)) }
    z <- linop:::function_coefficients(f, kind, "multiplier")
    expect_gt(z$reused, 0)
    expect_identical(length(seen), z$points)
    expect_identical(anyDuplicated(seen), 0L)
  }
})

test_that("additional probes detect a narrow feature missed by the two probes", {
  center <- cos(pi*(5.5+sqrt(2)/10)/(16+sqrt(2)/5))
  f <- function(x) exp(-((x-center)/.001)^2)
  grid <- sin(pi*seq(.5, -.5, length.out = 17))
  expect_lt(max(f(c(grid, -.3579, .1234))), .Machine$double.eps)
  expect_equal(max(f(linop:::coefficient_registry()$chebyshev$probes)), 1)
  expect_error(linop:::function_coefficients(f, "chebyshev", "feature", nmax = 129), "unresolved")
})
