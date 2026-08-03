## eigs() on the operator itself, where there is no matrix and the width of the
## truncation is what the call chooses. Ground truth is S0.6's single-site model:
## for the potential v at the origin the eigenvalue is sign(v) sqrt(v^2 + 4) and
## the eigenvector decays at (-|v| + sqrt(v^2 + 4)) / 2.

single <- function(v) linop_jacobi(diagonal = v)
single_value <- function(v) sign(v) * sqrt(v^2 + 4)
single_rate <- function(v) (-abs(v) + sqrt(v^2 + 4)) / 2
outer_end <- function(v) if (v > 0) "largest_algebraic" else "smallest_algebraic"
cert_status_of <- function(fit, check) linop:::cert_status(fit$certificate, check)

## The level the search accepts at: half the tolerance, in the currency the
## certificate grades in, which is absolute against the largest row sum.
accept_at <- function(v, tol = 1e-10) tol * max(abs(v) + 2, 2) / 2

test_that("the eigenvalue is recovered from an operator with no matrix", {
  for (v in c(0.5, 1, 2, 4, -1, -3)) {
    H <- single(v)
    fit <- eigs(H, k = 1, which = outer_end(v))
    truth <- single_value(v)
    ## the certified bound is a bound on the distance to the true spectrum, so
    ## the closed form has to sit inside it
    expect_lte(abs(fit$values - truth), max(fit$certificate$values$bound))
    expect_lt(abs(fit$values - truth), 1e-12)
    expect_equal(fit$certificate$overall, "pass")
    expect_equal(fit$certificate$subject, "finite section")
    expect_equal(fit$nconv, 1L)
    ## the rate the width was predicted with is the closed-form one
    expect_equal(decay_rate(H, fit$values), single_rate(v), tolerance = 1e-9)
  }
})

test_that("the width is chosen, recorded, and is a section of the operator", {
  H <- single(1)
  fit <- eigs(H, k = 1, which = "largest_algebraic")
  expect_type(fit$n, "integer")
  expect_identical(fit$n, fit$widths[length(fit$widths)])
  expect_true(fit$n > H$args$radius)
  expect_equal(fit$certificate$values$n, fit$n)
  expect_equal(fit$certificate$values$widths, fit$widths)
  expect_match(fit$method, "on a finite section")
  expect_output(print(fit), "section n = \\d+, from \\d+ widths?")
  expect_output(print(fit$certificate), "operator on l\\^2\\(Z\\); width")

  ## and it is exactly what eigs() on that section would have produced: the
  ## search picks the width and changes nothing else
  same <- eigs(finite_section(H, fit$n), k = 1, which = "largest_algebraic")
  expect_identical(fit$values, same$values)
  expect_identical(fit$vectors, same$vectors)
  expect_identical(fit$certificate$checks, same$certificate$checks)
})

test_that("the closed form is what picks the width, not a doubling", {
  ## The prediction rests on eta(n2) / eta(n1) = rho^(n2 - n1), which holds
  ## exactly while eta is above the level the computed eigenvector stores its own
  ## tail at. Measured, not assumed.
  eta <- function(v, n) {
    f <- eigs(finite_section(single(v), n), k = 1, which = outer_end(v))
    max(f$certificate$values$truncation)
  }
  for (v in c(0.3, 0.5)) {
    ratio <- eta(v, 60) / eta(v, 40)
    expect_equal(ratio / single_rate(v)^20, 1, tolerance = 1e-3)
  }

  ## So one measurement plus one division lands inside a factor of five of the
  ## acceptance level, at every rate. A doubling scheme cannot: from 41 it can
  ## only reach 82 or 164, and v = 0.3 needs 150.
  for (v in c(0.2, 0.3, 0.5, 1)) {
    fit <- eigs(single(v), k = 1, which = "largest_algebraic")
    got <- max(fit$certificate$values$truncation)
    expect_lte(got, accept_at(v))
    expect_gt(got, accept_at(v) / 20)
    expect_length(fit$widths, 2L)
  }
  expect_equal(eigs(single(0.3), k = 1, which = "largest_algebraic")$n, 150L)
})

test_that("a width that suffices at once is not widened", {
  ## v = 4 decays at 0.236, so the first section is already past the tolerance
  fit <- eigs(single(4), k = 1, which = "largest_algebraic")
  expect_length(fit$widths, 1L)
  expect_equal(fit$certificate$overall, "pass")
})

test_that("the first section leaves free tail beyond the window, not a fixed width", {
  ## A section that barely contains the window has nowhere for the eigenvector to
  ## have become geometric, so the rate it is predicted with is read off a vector
  ## that is not yet the right one. Measured on a 121-site well: n_start counted
  ## from the window costs one width, counted absolutely it costs three.
  W <- linop_jacobi(diagonal = rep(-1, 121))
  expect_equal(W$args$radius, 61L)
  wide <- eigs(W, k = 1, which = "smallest_algebraic")
  expect_equal(wide$widths[1L], 101L)
  expect_length(wide$widths, 1L)
  expect_equal(wide$certificate$overall, "pass")

  narrow <- eigs(W, k = 1, which = "smallest_algebraic",
                 section = list(n_start = 1))
  expect_equal(narrow$widths[1L], 62L)
  expect_gt(length(narrow$widths), 1L)
  ## both arrive, and they arrive at the same eigenvalue
  expect_equal(narrow$values, wide$values, tolerance = 1e-9)
  expect_equal(narrow$certificate$overall, "pass")
})

test_that("several pairs are found and each predicts its own width", {
  H <- linop_jacobi(diagonal = c(-1.5, -2, -1.5))
  fit <- eigs(H, k = 2, which = "smallest_algebraic")
  expect_equal(fit$nconv, 2L)
  expect_equal(fit$certificate$overall, "pass")
  ## both lie below the band and both decay
  expect_true(all(fit$values < -2))
  expect_true(all(decay_rate(H, fit$values) < 1))
  ## the width taken is the largest of the two predictions, so the slower one
  ## also clears the acceptance level
  expect_lte(max(fit$certificate$values$truncation),
             1e-10 * max(1.5 + 2, 2 + 2, 2) / 2)
})

test_that("a tail that stops falling with n stops the search and is named", {
  ## Past the level at which the finite eigensolve stores the tail at all, eta
  ## measures the eigenvector's own error and no width improves it. The search
  ## has to notice rather than widen to n_max.
  fit <- eigs(single(1), k = 1, which = "largest_algebraic", tol = 1e-18,
              section = list(n_max = 400))
  expect_lt(fit$n, 400L)
  expect_match(fit$certificate$dispatch$reason, "widening stopped shrinking the tail")
  ## the tail did reach the floor. The row is qualified rather than failed,
  ## because at a tolerance of 1e-18 the arithmetic floor is what covers it, and
  ## it says so rather than reporting a truncation that is not there
  expect_lt(max(fit$certificate$values$truncation), 1e-15)
  expect_equal(cert_status_of(fit, "truncation bound"), "qualified")
  detail <- fit$certificate$checks$detail[
    fit$certificate$checks$check == "truncation bound"]
  expect_match(detail, "the tail is at the arithmetic floor")
})

test_that("a width budget caps the search and the certificate reports it", {
  fit <- eigs(single(0.3), k = 1, which = "largest_algebraic",
              section = list(n_max = 60))
  expect_equal(fit$n, 60L)
  expect_match(fit$certificate$dispatch$reason, "n_max = 60")
  expect_equal(fit$certificate$overall, "fail")
  expect_equal(cert_status_of(fit, "truncation bound"), "fail")
  ## the eigensolve itself was fine; the section is what was short, and the two
  ## rows say so separately
  expect_equal(cert_status_of(fit, "finite residual"), "pass")
  expect_equal(cert_status_of(fit, "isolation"), "pass")
})

test_that("a value inside the band has no rate, and the search says so", {
  ## V = 0 is the free operator: its spectrum is the band and nothing else, so
  ## there is no eigenvalue at any width. With no rate to predict with the width
  ## doubles, and the certificate refuses the value rather than reporting one.
  fit <- eigs(single(0), k = 1, which = "largest_algebraic",
              section = list(n_max = 200))
  expect_equal(fit$widths, c(41L, 82L, 164L, 200L))
  expect_equal(fit$certificate$overall, "fail")
  expect_equal(cert_status_of(fit, "isolation"), "fail")
  expect_equal(cert_status_of(fit, "enclosure"), "not_checked")
  expect_true(is.na(decay_rate(single(0), fit$values)))
  expect_equal(fit$nconv, 0L)
})

test_that("shift-invert runs through the chosen section", {
  H <- single(1)
  fit <- eigs(H, k = 1, sigma = 2.3)
  expect_lt(abs(fit$values - sqrt(5)), 1e-12)
  expect_equal(fit$certificate$overall, "pass")
  expect_match(fit$method, "shift-invert on a finite section")
})

test_that("the refusals name what to do instead", {
  H <- single(1)
  expect_error(eigs(H, k = 1, v0 = rep(1, 83)), "does not take v0")
  expect_error(eigs(H, k = 1, v0 = rep(1, 83)), "finite_section\\(A, n\\)")
  expect_error(eigs(H, k = 1, section = list(n = 40)), "it does not take: n")
  expect_error(eigs(H, k = 1, section = list(n_max = 5)),
               "n_max = 5 is below the first width")
  expect_error(eigs(H, k = 1, section = list(n_start = 0)), "at least 1")
  expect_error(eigs(H, k = 1, section = 40), "section must be a list")

  ## and the knob has no meaning on an operator whose truncation is already fixed
  expect_error(eigs(finite_section(H, 20), k = 1, section = list(n_start = 5)),
               "section applies to an operator on a sequence space")
  expect_error(eigs(finite_section(H, 20), k = 1, section = list(n_start = 5)),
               "41 x 41")

  ## an infinite operator whose node type has no route to a matrix is refused as
  ## it always was
  expect_error(eigs(H %*% H, k = 1), "needs an operator with finite dimensions")
})

test_that("svds still refuses an operator on a sequence space", {
  ## The hook is eigs()'s, and a singular value of a self-adjoint operator on
  ## l^2(Z) is not a second thing this layer answers.
  expect_error(svds(single(1), k = 1), "needs an operator with finite dimensions")
  expect_error(solve(single(1), rep(1, 5)), "needs an operator with finite dimensions")
})
