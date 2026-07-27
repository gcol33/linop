## Gate 1: every conditional construction rule records depends_on, and
## evidence_satisfies() recursing through a composite of user_declaration inputs
## fails a requirement that the inputs would have failed directly.

## CG's default requirement for positive_definite: construction, a trusted
## adapter contract, or a theorem. Not ten positive Rayleigh quotients.
cg_requirement <- function() {
  requirement(sources = c("construction", "adapter_contract", "theorem"),
              guarantees = "identity", min_confidence = 1)
}

test_that("evidence rejects unknown fields", {
  expect_error(evidence("guesswork", "identity"), "source must be one of")
  expect_error(evidence("probe", "vibes"), "guarantee must be one of")
  expect_error(requirement(sources = "telepathy"), "unknown evidence source")
})

test_that("a capability with a known value needs evidence", {
  expect_error(capability(TRUE), "needs evidence")
  expect_silent(capability(NA))
  expect_null(capability(NA)$evidence)
})

test_that("NA is never read as FALSE", {
  A <- linop(rmat(4, 4, seed = 1))
  expect_true(is.na(capv(A, "positive_definite")))
  expect_false(isTRUE(capv(A, "positive_definite")))
  ## and the three-valued conjunction keeps it
  expect_true(is.na(linop:::tv_and(TRUE, NA)))
  expect_false(linop:::tv_and(FALSE, NA))
})

test_that("a user declaration does not satisfy CG's requirement", {
  n <- 4
  A <- linop(function(X) X, adjoint = function(X) X, dim = c(n, n),
             properties = c(hermitian = TRUE, positive_definite = TRUE))
  ev <- linop:::cape(A, "positive_definite")
  expect_equal(ev$source, "user_declaration")
  expect_false(evidence_satisfies(ev, cg_requirement()))
})

test_that("construction evidence from crossprod satisfies CG's requirement", {
  A <- linop(rmat(5, 4, seed = 2))
  g <- adjoint(A) %*% A
  ev <- linop:::cape(g, "positive_semidefinite")
  expect_equal(ev$source, "construction")
  expect_true(evidence_satisfies(ev, cg_requirement()))
})

test_that("depends_on is recorded for conditional rules and empty for unconditional ones", {
  ## unconditional: A^H A is hermitian whatever is known about A
  A <- linop(zmat(4, 3, seed = 3))
  ev_uncond <- linop:::cape(adjoint(A) %*% A, "hermitian")
  expect_length(ev_uncond$depends_on, 0)

  ## conditional: A + B is hermitian only if both are
  n <- 4
  H1 <- linop(function(X) X, adjoint = function(X) X, dim = c(n, n),
              properties = c(hermitian = TRUE))
  H2 <- linop(function(X) 2 * X, adjoint = function(X) 2 * X, dim = c(n, n),
              properties = c(hermitian = TRUE))
  ev_cond <- linop:::cape(H1 + H2, "hermitian")
  expect_equal(ev_cond$source, "construction")
  expect_gt(length(ev_cond$depends_on), 0)
})

test_that("the laundering case: construction over declarations still fails", {
  ## This is the case section 5.3 exists to prevent. Both operands carry
  ## unconditional construction evidence at the top, but the composite is built
  ## by a conditional rule, so the declarations underneath must remain visible.
  n <- 4
  mk_declared <- function(s) {
    linop(function(X) s * X, adjoint = function(X) s * X, dim = c(n, n),
          properties = c(hermitian = TRUE, positive_definite = TRUE))
  }
  A <- mk_declared(1); B <- mk_declared(2)

  ## direct: fails, as it must
  expect_false(evidence_satisfies(linop:::cape(A, "positive_definite"), cg_requirement()))

  ## composite by the conditional sum rule: must also fail
  S <- A + B
  ev <- linop:::cape(S, "positive_definite")
  expect_equal(ev$source, "construction")
  expect_false(evidence_satisfies(ev, cg_requirement()),
               info = "a sum of declared-PD operators laundered its evidence")

  ## and one level deeper
  S2 <- 2 * (A + B)
  ev2 <- linop:::cape(S2, "positive_definite")
  expect_false(evidence_satisfies(ev2, cg_requirement()))
})

test_that("evidence_satisfies recurses to arbitrary depth", {
  deep <- evidence("user_declaration", "identity", 1)
  for (i in 1:5) deep <- evidence("construction", "identity", 1, depends_on = list(deep))
  expect_false(evidence_satisfies(deep, cg_requirement()))

  clean <- evidence("theorem", "identity", 1)
  for (i in 1:5) clean <- evidence("construction", "identity", 1, depends_on = list(clean))
  expect_true(evidence_satisfies(clean, cg_requirement()))
})

test_that("confidence is enforced", {
  ev <- evidence("computation", "probabilistic_bound", 0.99)
  expect_false(evidence_satisfies(ev, requirement(min_confidence = 1)))
  expect_true(evidence_satisfies(ev, requirement(min_confidence = 0.95)))
  expect_false(evidence_satisfies(evidence("probe", "heuristic", NA),
                                  requirement(min_confidence = 0.5)))
})

test_that("verify adds a probe record beside a declaration and does not replace it", {
  n <- 4
  A <- linop(function(X) X, adjoint = function(X) X, dim = c(n, n),
             properties = c(hermitian = TRUE))
  before <- linop:::cape(A, "hermitian")$source
  cert <- verify(A)
  expect_equal(before, "user_declaration")
  ## the operator is immutable, so the declaration is untouched
  expect_equal(linop:::cape(A, "hermitian")$source, "user_declaration")
  ## and the probe result is reported separately, as heuristic evidence
  row <- cert$checks[cert$checks$check == "declared capabilities", ]
  expect_equal(row$source, "probe")
  expect_equal(row$guarantee, "heuristic")
})
