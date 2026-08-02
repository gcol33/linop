## Gate 1: the exported name set is asserted against the budget. Adding an export
## fails this test until the budget below is edited, which is the only mechanism
## that reliably stops a public surface from growing by accretion.
##
## Section 1.1 budgets the v0.1 surface. eigs() and svds() have now arrived and
## cost one export each, which is the whole of what Phase 2 spends: they are new
## names rather than base generics, so unlike solve() there is nothing to register
## a method on.
##
## solve() cost no export at all. It is a base generic, so an S3 method reaches it
## the way t(), crossprod() and %*% already do, and section 1.1's "keeps the verb
## budget" turns out to be literal for it. The last test in this file is where
## that is asserted, beside the other generics that work the same way.
##
## print() and summary() on the two spectral results are the same mechanism and
## cost nothing either.
##
## Eleven names came OUT in the one-package redesign, which is the first time this
## budget has moved down. dev_notes/one-package-and-the-abstractions-that-were-holding-a-boundary.md
## records the reassessment; the short version is that every one of them existed to
## hold a package boundary that no longer exists:
##
##   provenance, set_provenance, strip_provenance and the four generics
##     -- an opaque envelope buys core not needing to know a payload's structure,
##        which is worth nothing once core owns the structure. A truncation is a
##        node with a child instead, which is structural rather than annotated.
##   build_certificate, cert_rows
##     -- internal for all of Phase 2, public for one day, for one consumer in
##        another package. The public certificate surface is a reader, not a builder.
##   linop_register_node, linop_nodes
##     -- how a package outside this one would add a composition type. There is no
##        package outside this one, and a user bringing their own operator writes
##        linop.<class>(), which needs no registry.

BUDGET <- list(
  common = c("linop", "adjoint", "verify", "eigs", "svds"),

  ## constructors for the leaves a user builds directly
  constructors = c("linop_eye", "linop_scaling"),

  ## An operator on a sequence space, the truncation that makes it computable,
  ## and the rate the closed forms in its certificate are written in. Three names
  ## for a whole layer, because the truncation returns an ordinary linop and
  ## everything downstream of it is the surface that already existed.
  sequence_space = c("linop_jacobi", "finite_section", "decay_rate"),

  advanced = c("preconditioner", "as_preconditioner", "as_preconditioner_inverse",
               "collapse", "explain", "as_sparse"),

  developer = c("evidence", "capability", "requirement", "evidence_satisfies",
                "cap", "dtype", "is_linop", "linop_cost")
)

PLANNED_PHASE2 <- c("solver", "linsolve",
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

test_that("the names Phase 2 deliberately keeps private are not exported", {
  actual <- getNamespaceExports("linop")
  early <- intersect(PLANNED_PHASE2, actual)
  expect_equal(early, character(0),
               info = paste0("names section 1.1 keeps private were exported: ",
                             paste(early, collapse = ", ")))
})

test_that("the budget stays small", {
  expect_lte(length(BUDGET$common), 6L)
  ## The whole surface, so a slow drift upward is visible in the diff. The cap
  ## came down from 32 with the eleven names the one-package redesign removed.
  ## Three of the four the Hilbert layer adds are in; certificate() is the last.
  expect_lte(length(unlist(BUDGET, use.names = FALSE)), 25L)
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

  ## and solve(), which is the whole of Phase 2's verb surface
  S <- linop(diag(c(2, 3, 4)), properties = c(hermitian = TRUE,
                                              positive_definite = TRUE))
  x <- solve(S, c(2, 3, 4))
  expect_equal(as.numeric(x), c(1, 1, 1))
  ## registered as a method, not exported as a name
  expect_false("solve" %in% getNamespaceExports("linop"))
  expect_false(is.null(utils::getS3method("solve", "linop", optional = TRUE)))
})
