## Gate 1: linsolve conformance is a separate suite, and running the linop suite
## against a linsolve errors rather than failing purity. Section 4.

mk_exact_solve <- function(M) {
  linop:::new_linsolve(function(R) solve(M, R), dim = dim(M),
                       fidelity = "exact", determinacy = "fixed", method = "lu")
}

## A warm-started CG stands in for the history-dependent case: same input,
## different output, by design.
mk_history_solve <- function(n) {
  last <- NULL
  linop:::new_linsolve(
    function(R) { out <- if (is.null(last)) R else (R + last) / 2; last <<- out; out },
    dim = c(n, n), fidelity = "variable_inexact", determinacy = "history_dependent",
    method = "warm-started cg")
}

test_that("running the operator suite against a linsolve is an error", {
  S <- mk_exact_solve(diag(c(2, 3, 4)))
  expect_error(linop:::verify.linop(S), "not a linop")
  expect_error(linop:::verify.linop(S), "conformance suite does not apply")
})

test_that("verify() dispatches a linsolve to its own suite", {
  S <- mk_exact_solve(diag(c(2, 3, 4)))
  cert <- verify(S)
  expect_s3_class(cert, "linop_certificate")
  expect_equal(cert$subject, "linsolve")
  expect_false(cert$overall == "fail")
})

test_that("a history-dependent solve is not judged by purity", {
  S <- mk_history_solve(4)
  cert <- verify(S)
  ## it must not fail: it never claimed repeatability
  expect_false(cert$overall == "fail")
  expect_equal(cert$checks$status[cert$checks$check == "repeatability"], "not_checked")
  expect_equal(cert$checks$status[cert$checks$check == "additivity"], "not_checked")
})

test_that("a solve that claims fixed determinacy is held to it", {
  bad <- linop:::new_linsolve(function(R) R * runif(1), dim = c(3, 3),
                              fidelity = "linear_approximation", determinacy = "fixed")
  cert <- verify(bad)
  expect_equal(cert$checks$status[cert$checks$check == "repeatability"], "fail")
})

test_that("contract acceptance keeps history-dependent solves away from ordinary Arnoldi", {
  exact <- mk_exact_solve(diag(c(2, 3, 4)))
  hist <- mk_history_solve(4)
  accepted <- c("exact", "fixed_linear")

  expect_silent(linop:::check_solve_contract(exact, accepted, "shift-invert Lanczos"))
  expect_error(linop:::check_solve_contract(hist, accepted, "shift-invert Lanczos"),
               "does not accept this solve object")
})

test_that("variable_inexact is accepted only with an error bound", {
  no_bound <- linop:::new_linsolve(function(R) R, dim = c(3, 3),
                                   fidelity = "variable_inexact",
                                   determinacy = "input_dependent", error_bound = FALSE)
  with_bound <- linop:::new_linsolve(function(R) R, dim = c(3, 3),
                                     fidelity = "variable_inexact",
                                     determinacy = "input_dependent", error_bound = TRUE)
  accepted <- c("exact", "fixed_linear", "variable_inexact")
  expect_error(linop:::check_solve_contract(no_bound, accepted, "inexact Krylov"))
  expect_silent(linop:::check_solve_contract(with_bound, accepted, "inexact Krylov"))
})

test_that("linsolve rejects unknown contract values", {
  expect_error(linop:::new_linsolve(function(R) R, c(2, 2), "approximate", "fixed"),
               "fidelity must be one of")
  expect_error(linop:::new_linsolve(function(R) R, c(2, 2), "exact", "sometimes"),
               "determinacy must be one of")
})

## ------------------------------------------------------------ preconditioner --

test_that("each solver refuses a preconditioner lacking a property it requires", {
  flexible <- preconditioner(function(R) R, fixed = FALSE,
                             hermitian = TRUE, positive_definite = TRUE)
  indefinite <- preconditioner(function(R) R, fixed = TRUE,
                               hermitian = TRUE, positive_definite = FALSE)
  unknown_pd <- preconditioner(function(R) R, fixed = TRUE,
                               hermitian = TRUE, positive_definite = NA)
  good <- preconditioner(function(R) R, fixed = TRUE,
                         hermitian = TRUE, positive_definite = TRUE)

  ## one row per line of the section 4.3 table
  for (m in c("cg", "minres")) {
    expect_error(linop:::check_preconditioner(flexible, m), "fixed = TRUE")
    expect_error(linop:::check_preconditioner(indefinite, m), "positive_definite = TRUE")
    expect_error(linop:::check_preconditioner(unknown_pd, m), "positive_definite = TRUE")
    expect_silent(linop:::check_preconditioner(good, m))
  }
  for (m in c("gmres", "bicgstab", "lsqr", "lsmr")) {
    expect_error(linop:::check_preconditioner(flexible, m), "fixed = TRUE")
    expect_silent(linop:::check_preconditioner(indefinite, m))
  }
  ## FGMRES is the one that accepts a flexible preconditioner
  expect_silent(linop:::check_preconditioner(flexible, "fgmres"))
})

test_that("an unknown property is refused, not assumed", {
  ## NA must never be read as TRUE
  p <- preconditioner(function(R) R, fixed = TRUE, hermitian = TRUE,
                      positive_definite = NA)
  expect_error(linop:::check_preconditioner(p, "cg"), "positive_definite")
  q <- preconditioner(function(R) R, fixed = TRUE, hermitian = NA,
                      positive_definite = TRUE)
  expect_error(linop:::check_preconditioner(q, "cg"), "hermitian")
})

test_that("the error names the flag and suggests the alternative", {
  flexible <- preconditioner(function(R) R, fixed = FALSE)
  expect_error(linop:::check_preconditioner(flexible, "cg"), "fgmres")
})

test_that("both preconditioner directions produce the same internal form", {
  M <- diag(c(2, 4, 8))
  Pm <- as_preconditioner(linop(M), solver = function(op, R) solve(as.matrix(op), R))
  Pi <- as_preconditioner_inverse(linop(solve(M)))
  r <- matrix(c(1, 1, 1), 3, 1)
  expect_equal(Pm$apply_inverse(r), Pi$apply_inverse(r), tolerance = 1e-12,
               ignore_attr = TRUE)
})

test_that("as_preconditioner refuses to invent a solver", {
  expect_error(as_preconditioner(linop(diag(3))), "needs a solver")
})
