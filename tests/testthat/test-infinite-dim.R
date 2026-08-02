## Step 2 of the one-package redesign: a dimension may be Inf, because an operator
## on a sequence space has an infinite index set. Inf is the extent, not a flag, so
## the shape arithmetic that already exists has to keep working on it, and every
## numeric path has to refuse it by name rather than fail inside an allocation.
##
## The operator used here is a placeholder for the Jacobi operator step 3 brings:
## it has structure and no apply, which is exactly the case the contract must hold
## for.

inf_op <- function(dtype = "double", caps = c(hermitian = TRUE, symmetric = TRUE)) {
  linop:::new_linop("fun", c(Inf, Inf), dtype,
                    do.call(linop:::new_caps, linop:::normalise_properties(caps)),
                    args = list(apply = function(X) stopf("unreachable"),
                                adjoint = NULL))
}

test_that("a dimension may be Inf, and stays double so it can be", {
  A <- inf_op()
  expect_equal(dim(A), c(Inf, Inf))
  expect_type(dim(A), "double")
  expect_false(linop:::is_finite_dim(A))

  ## and a finite dimension is still integer, so nothing that compared against
  ## integers before now compares against doubles
  B <- linop(rmat(4, 3, seed = 1))
  expect_type(dim(B), "integer")
  expect_true(linop:::is_finite_dim(B))
})

test_that("the dimension contract refuses what is not a dimension", {
  mk <- function(d) linop:::new_linop("fun", d, "double", args = list(apply = identity))
  expect_error(mk(c(-1, 2)), "non-negative")
  expect_error(mk(c(2.5, 2)), "whole number")
  expect_error(mk(c(NA, 2)), "whole number or Inf")
  expect_error(mk(c(2, 2, 2)), "whole number or Inf")
  expect_error(mk("4"), "whole number or Inf")
  ## -Inf is not a size
  expect_error(mk(c(-Inf, 2)), "non-negative")
})

test_that("length() is NA rather than a number no integer can hold", {
  expect_true(is.na(length(inf_op())))
  expect_type(length(inf_op()), "integer")
  expect_equal(length(linop(rmat(4, 3, seed = 2))), 12L)
})

test_that("every numeric path refuses an infinite operator by name", {
  A <- inf_op()
  expect_error(A %*% matrix(1, 3, 1), "apply\\(\\) needs an operator with finite")
  expect_error(as.matrix(A), "as.matrix\\(\\) needs an operator with finite")
  expect_error(as.array(A), "as.matrix\\(\\) needs an operator with finite")
  expect_error(collapse(A), "collapse\\(\\) needs an operator with finite")
  expect_error(verify(A), "verify\\(\\) needs an operator with finite")
  expect_error(eigs(A, k = 1), "eigs\\(\\) needs an operator with finite")
  expect_error(svds(A, k = 1), "svds\\(\\) needs an operator with finite")
  expect_error(linop:::norm2(A), "norm2\\(\\) needs an operator with finite")
  expect_error(solve(A, matrix(1, 3, 1)), "solve\\(\\) needs an operator with finite")
})

test_that("the refusal names the route out, since there is one", {
  ## The message has to say what to do, because an operator on a sequence space is
  ## not a mistake: it is the thing the caller meant to build.
  err <- tryCatch(as.matrix(inf_op()), error = conditionMessage)
  expect_match(err, "Inf x Inf")
  expect_match(err, "no matrix and no block to act on")
  expect_match(err, "finite_section\\(\\)")
})

test_that("the algebra composes infinite operators without touching a block", {
  A <- inf_op()
  B <- inf_op()

  expect_equal(dim(t(A)), c(Inf, Inf))
  expect_equal(dim(adjoint(A)), c(Inf, Inf))
  expect_equal(dim(Conj(A)), c(Inf, Inf))
  expect_equal(dim(2 * A), c(Inf, Inf))
  expect_equal(dim(A + B), c(Inf, Inf))
  expect_equal(dim(A %*% B), c(Inf, Inf))

  ## and the composite is still an operator with no block, refused the same way
  expect_error(as.matrix(A %*% B), "needs an operator with finite")
})

test_that("squareness and capability closure work on Inf without a special case", {
  ## Inf == Inf, so the hermitian squareness check passes on the same line that
  ## checks it for a finite operator.
  expect_silent(inf_op())
  expect_true(isTRUE(cap(inf_op(), "hermitian")$value))
  ## a rectangular infinite operator cannot be declared hermitian, same rule
  expect_error(
    linop:::new_linop("fun", c(Inf, 4), "double",
                      do.call(linop:::new_caps,
                              linop:::normalise_properties(c(hermitian = TRUE))),
                      args = list(apply = identity)),
    "must be square")
})

test_that("a mixed shape is a dimension too, and composes by the usual rule", {
  ## An operator from a finite space into a sequence space. Nothing here needs it
  ## yet; the contract admits it because the conformability rule is the same one.
  E <- linop:::new_linop("fun", c(Inf, 4), "double", args = list(apply = identity))
  F <- linop:::new_linop("fun", c(4, Inf), "double", args = list(apply = identity))
  expect_equal(dim(E %*% F), c(Inf, Inf))
  expect_equal(dim(F %*% E), c(4, 4))
  ## and the finite composite is still refused, because its factors have no block
  expect_error(as.matrix(F %*% E), "needs an operator with finite")
  expect_error(E %*% E, "non-conformable product: Inf x 4 times Inf x 4")
})

test_that("printing an infinite operator says Inf rather than erroring", {
  ## sprintf("%d", Inf) is an error, so every shape message goes through fmt_dim().
  A <- inf_op()
  expect_output(print(A), "<linop> Inf x Inf, double")
  expect_output(print(summary(A)), "<linop> Inf x Inf, double, root node 'fun'")
  expect_output(explain(A), "operator Inf x Inf")
  expect_output(print(A %*% A), "Inf x Inf")
})

test_that("a non-conformable sum of infinite operators reports both shapes", {
  A <- inf_op()
  E <- linop:::new_linop("fun", c(Inf, 4), "double", args = list(apply = identity))
  expect_error(A + E, "non-conformable sum: Inf x Inf and Inf x 4")
})

test_that("cost is Inf, which is what an ordering heuristic should say", {
  ## linop_cost() orders expressions. An operator with no block cannot be cheaper
  ## than one that has one, and Inf sorts that way without a special case.
  expect_equal(linop_cost(inf_op()), Inf)
  expect_equal(linop_cost(2 * inf_op()), Inf)
})
