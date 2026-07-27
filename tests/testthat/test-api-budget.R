## Gate 1: the exported name set is asserted against the budget. Adding an export
## fails this test until the budget below is edited, which is the only mechanism
## that reliably stops a public surface from growing by accretion.
##
## Section 1.1 budgets the v0.1 surface. Phase 1 delivers the core, so solve(),
## eigs() and svds() are absent here and arrive with Phase 2. Their names are
## listed as `planned` so the intent stays visible and the test still fails if
## one appears early without a deliberate edit.

BUDGET <- list(
  common = c("linop", "adjoint", "verify"),

  ## constructors for the leaves a user builds directly
  constructors = c("linop_eye", "linop_scaling"),

  advanced = c("preconditioner", "as_preconditioner", "as_preconditioner_inverse",
               "collapse", "explain", "as_sparse"),

  developer = c("linop_register_node", "linop_nodes",
                "evidence", "capability", "requirement", "evidence_satisfies",
                "cap", "dtype", "is_linop", "linop_cost",
                "provenance", "set_provenance", "strip_provenance",
                "provenance_lift", "provenance_refine",
                "provenance_original_residual", "provenance_summary")
)

PLANNED_PHASE2 <- c("solve", "eigs", "svds", "solver", "linsolve",
                    "spectral_approximation", "plan_eigs", "linop_register_backend")

test_that("the exported surface matches the budget exactly", {
  actual <- sort(getNamespaceExports("linop"))
  actual <- actual[!grepl("^\\.", actual)]
  expected <- sort(unlist(BUDGET, use.names = FALSE))

  extra <- setdiff(actual, expected)
  missing <- setdiff(expected, actual)

  expect_equal(extra, character(0),
               info = paste0("exports not in the budget: ", paste(extra, collapse = ", "),
                             "\n  Add them to BUDGET in this file, deliberately."))
  expect_equal(missing, character(0),
               info = paste0("budgeted but not exported: ", paste(missing, collapse = ", ")))
})

test_that("Phase 2 names are not exported yet", {
  actual <- getNamespaceExports("linop")
  early <- intersect(PLANNED_PHASE2, actual)
  expect_equal(early, character(0),
               info = paste0("Phase 2 names exported during Phase 1: ",
                             paste(early, collapse = ", ")))
})

test_that("the budget stays small", {
  expect_lte(length(BUDGET$common), 6L)
  ## the whole surface, so a slow drift upward is visible in the diff
  expect_lte(length(unlist(BUDGET, use.names = FALSE)), 32L)
})

test_that("R's own generics work without being exported", {
  A <- linop(rmat(4, 3, seed = 1))
  ## these come from S3 method registration, not from the export list
  expect_s3_class(t(A), "linop")
  expect_s3_class(Conj(A), "linop")
  expect_s3_class(crossprod(A), "linop")
  expect_s3_class(tcrossprod(A), "linop")
  expect_s3_class(A + A, "linop")
  expect_equal(dim(A), c(4L, 3L))
  ## nrow and ncol are not generic in base R and need no methods (S0.1 section 5)
  expect_equal(nrow(A), 4L)
  expect_equal(ncol(A), 3L)
})
