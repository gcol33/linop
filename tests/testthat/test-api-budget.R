## The exported name set is asserted against the budget. Adding an export
## fails this test until the budget below is edited, which is the only mechanism
## that reliably stops a public surface from growing by accretion.
##
## The budget covers the v0.1 surface. eigs() and svds() cost one export each:
## they are new
## names rather than base generics, so unlike solve() there is nothing to register
## a method on.
##
## solve() cost no export at all. It is a base generic, so an S3 method reaches it
## the way t(), crossprod() and %*% already do, so keeping the verb
## budget is literal for it. The last test in this file is where
## that is asserted, beside the other generics that work the same way.
##
## print() and summary() on the two spectral results are the same mechanism and
## cost nothing either.
##
## Eleven names came OUT in the one-package redesign, which is the first time this
## budget has moved down. Every one of them existed to
## hold a package boundary that no longer exists:
##
##   provenance, set_provenance, strip_provenance and the four generics
##     -- an opaque envelope buys core not needing to know a payload's structure,
##        which is worth nothing once core owns the structure. A truncation is a
##        node with a child instead, which is structural rather than annotated.
##   build_certificate, cert_rows
##     -- internal before, public for one day, for one consumer in
##        another package. The public certificate surface is a reader, not a builder.
##   linop_register_node, linop_nodes
##     -- how a package outside this one would add a composition type. There is no
##        package outside this one, and a user bringing their own operator writes
##        linop.<class>(), which needs no registry.

BUDGET <- list(
  ## certificate() is the reader that replaced the two builders. Every verb here
  ## returns a certificate somewhere in its result, and one accessor over the
  ## four shapes is what a public certificate surface is.
  common = c("linop", "adjoint", "verify", "eigs", "svds", "certificate"),

  ## constructors for the leaves a user builds directly
  constructors = c("linop_eye", "linop_scaling"),

  ## An operator on a sequence space, the truncation that makes it computable,
  ## and the rate the closed forms in its certificate are written in. Three names
  ## for a whole layer, because the truncation returns an ordinary linop and
  ## everything downstream of it is the surface that already existed.
  sequence_space = c("linop_jacobi", "finite_section", "decay_rate"),

  ## The structured node types, decided 2026-11-15 with the cap raised
  ## from 25 to 30.
  ## Five constructors, because base R has no generic to hang them on; the
  ## stacks come through cbind() and rbind() methods and the exact inverse
  ## through solve(a), at no export cost.
  structured = c("linop_kron", "linop_blockdiag", "linop_lowrank", "linop_perm",
                 "linop_power"),

  ## The one front door to operators on function spaces, decided
  ## 2026-11-15 with the cap raised from 30 to 31: every space is a kind of it,
  ## so later spaces add arguments rather than names.
  function_space = "linop_space",

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

test_that("the names deliberately kept private are not exported", {
  actual <- getNamespaceExports("linop")
  early <- intersect(PLANNED_PHASE2, actual)
  expect_equal(early, character(0),
               info = paste0("names kept private were exported: ",
                             paste(early, collapse = ", ")))
})

test_that("the budget stays small", {
  expect_lte(length(BUDGET$common), 6L)
  ## The whole surface, so a slow drift upward is visible in the diff. The cap
  ## came down from 32 with the eleven names the one-package redesign removed,
  ## then went to 30 for the five structured constructors and to 31 for the one
  ## function-space front door, and all 31 are spent.
  ## The next name is a decision with a note behind it, in either direction.
  expect_lte(length(unlist(BUDGET, use.names = FALSE)), 31L)
  expect_equal(length(unlist(BUDGET, use.names = FALSE)), 31L)
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
  ## nrow and ncol are not generic in base R and need no methods
  expect_equal(nrow(A), 4L)
  expect_equal(ncol(A), 3L)

  ## and solve(), which is the whole of the solver verb surface
  S <- linop(diag(c(2, 3, 4)), properties = c(hermitian = TRUE,
                                              positive_definite = TRUE))
  x <- solve(S, c(2, 3, 4))
  expect_equal(as.numeric(x), c(1, 1, 1))
  ## registered as a method, not exported as a name
  expect_false("solve" %in% getNamespaceExports("linop"))
  expect_false(is.null(utils::getS3method("solve", "linop", optional = TRUE)))

  ## the stacks and the exact inverse, which are node types without names
  expect_s3_class(cbind(A, A), "linop")
  expect_s3_class(rbind(A, A), "linop")
  expect_s3_class(solve(S), "linop")
  expect_false(any(c("cbind", "rbind") %in% getNamespaceExports("linop")))
})
