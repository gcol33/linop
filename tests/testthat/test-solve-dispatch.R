## Section 1.1: one verb over the seven methods. What auto is allowed to choose
## and on what evidence, what naming a method means instead, where the verb
## refuses, and what the result carries back.

cert_status <- linop:::cert_status

chosen <- function(x) attr(x, "certificate")$dispatch$chosen
why <- function(x) attr(x, "certificate")$dispatch$reason

## An operator whose definiteness rests on something other than the caller
## saying so. No constructor in the package produces one today, which is the
## subject of its own test below; this reaches the branch so the branch is
## covered rather than assumed.
spd_by_construction <- function(M) {
  A <- linop(M)
  A$caps$hermitian <- capability(TRUE, evidence("construction", "identity", 1))
  A$caps$positive_definite <- capability(TRUE, evidence("construction", "identity", 1))
  A
}

## -------------------------------------------------------------- what auto picks

test_that("auto picks cg for a hermitian positive definite operator, and says why", {
  n <- 30
  M <- spd_prescribed(n, seq(1, 6, length.out = n))
  A <- spd_by_construction(M)
  set.seed(1)
  b <- matrix(stats::rnorm(n), n, 1)

  x <- solve(A, b, tol = 1e-11)
  expect_equal(chosen(x), "cg")
  expect_match(why(x), "hermitian positive definite")
  expect_equal(as.matrix(x), solve(M, b), tolerance = 1e-9, ignore_attr = TRUE)
})

test_that("auto picks minres where definiteness is not established", {
  ## A dense leaf computes symmetry from the entries, which is an identity on
  ## data the operator holds, and computes no definiteness, which would cost a
  ## factorization. So the capability that is established selects the method.
  n <- 30
  M <- spd_prescribed(n, seq(1, 6, length.out = n))
  A <- linop(M)
  expect_true(isTRUE(capv(A, "hermitian")))
  expect_equal(linop:::cape(A, "hermitian")$source, "computation")
  expect_true(is.na(capv(A, "positive_definite")))
  set.seed(2)
  b <- matrix(stats::rnorm(n), n, 1)

  x <- solve(A, b, tol = 1e-11, maxit = 2000L)
  expect_equal(chosen(x), "minres")
  expect_match(why(x), "definiteness is not established")
  expect_equal(as.matrix(x), solve(M, b), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("no constructor in the package reaches the cg branch of auto today", {
  ## Measured, and recorded here rather than left to be rediscovered. The cg
  ## branch requires positive_definite behind construction, adapter_contract,
  ## theorem or computation evidence. linop_scaling() and linop_eye() establish
  ## it that way and are diagonal, so the direct route takes them first; every
  ## other route to positive_definite in the package is a bare user declaration,
  ## which auto does not act on. The branch is correct and its supply is empty.
  ##
  ## When a constructor starts proving definiteness, this test is the one that
  ## fails, and it should be replaced rather than adjusted.
  n <- 12
  Mspd <- spd_prescribed(n, seq(1, 4, length.out = n))
  ops <- list(
    eye = linop_eye(n),
    scaling = linop_scaling(seq(1, 4, length.out = n)),
    dense_bare = linop(Mspd),
    dense_declared = linop(Mspd, properties = c(hermitian = TRUE,
                                                positive_definite = TRUE)),
    fun_declared = laplacian_1d(n),
    gram = adjoint(linop(rmat(2 * n, n, seed = 1))) %*% linop(rmat(2 * n, n, seed = 1)))
  picked <- vapply(ops, function(A) linop:::auto_method(A)$method, character(1))
  expect_false(any(picked == "cg"))
  expect_equal(unname(picked[c("eye", "scaling")]), c("direct", "direct"))
  expect_equal(unname(picked[c("dense_bare", "gram")]), c("minres", "minres"))
  expect_equal(unname(picked[c("dense_declared", "fun_declared")]), c("gmres", "gmres"))
})

test_that("auto falls back to gmres, which is what makes it total", {
  n <- 40
  A <- convdiff_1d(n, 0.3)
  expect_true(is.na(capv(A, "hermitian")))
  set.seed(3)
  b <- matrix(stats::rnorm(n), n, 1)

  x <- solve(A, b, tol = 1e-11, maxit = 2000L, control = list(restart = n))
  expect_equal(chosen(x), "gmres")
  expect_match(why(x), "requires none")
  expect_equal(as.matrix(x), convdiff_1d_solve(n, 0.3, b),
               tolerance = 1e-9, ignore_attr = TRUE)
})

test_that("auto takes the direct route on a diagonal operator, at no iterations", {
  d <- c(2, -3, 0.5, 4, 10)
  A <- linop_scaling(d)
  set.seed(4)
  b <- matrix(stats::rnorm(5), 5, 1)

  fit <- solve(A, b, details = TRUE)
  expect_equal(fit$method, "direct")
  expect_equal(fit$iterations, 0L)
  expect_match(fit$reason, "diagonal and of full rank")
  expect_equal(fit$x, b / d, tolerance = 1e-14)
  expect_equal(cert_status(fit$certificate, "residual"), "pass")

  ## and the identity leaf reaches it too, since it declares both
  I <- linop_eye(6)
  expect_equal(chosen(solve(I, rep(3, 6))), "direct")
})

test_that("auto never acts on a capability whose evidence it does not accept", {
  ## The laundering case of section 5.3 reaching dispatch. A sum of two declared
  ## positive definite operators is positive definite, and the argument for it
  ## bottoms out in user declarations, which is below what auto requires. The
  ## same operator solved by cg on request runs, because naming the method is the
  ## caller asserting their own declaration.
  n <- 20
  M <- spd_prescribed(n, seq(1, 4, length.out = n))
  P <- linop(M, properties = c(hermitian = TRUE, positive_definite = TRUE))
  S <- P + P
  expect_true(isTRUE(capv(S, "positive_definite")))

  set.seed(5)
  b <- matrix(stats::rnorm(n), n, 1)
  x <- solve(S, b, tol = 1e-11, maxit = 500L)
  expect_false(chosen(x) == "cg")

  named <- solve(S, b, method = "cg", tol = 1e-11, maxit = 500L)
  expect_equal(chosen(named), "cg")
  expect_equal(as.numeric(named), as.numeric(solve(2 * M, b)), tolerance = 1e-8)
})

test_that("the same evidence source can pass or fail, and depends_on is why", {
  ## Both operators below carry positive_definite behind a construction source,
  ## and only one of them is a proof. The sum's construction evidence depends on
  ## two user declarations, so evidence_satisfies() recurses and refuses it; the
  ## other depends on nothing, so it holds unconditionally. Evidence is three
  ## fields and a dependency chain, not a rank, and this is where dispatch reads
  ## the chain.
  n <- 20
  M <- spd_prescribed(n, seq(1, 4, length.out = n))
  bare <- spd_by_construction(M)
  P <- linop(M, properties = c(hermitian = TRUE, positive_definite = TRUE))
  laundered <- P + P

  expect_equal(linop:::cape(bare, "positive_definite")$source, "construction")
  expect_equal(linop:::cape(laundered, "positive_definite")$source, "construction")
  expect_true(isTRUE(capv(laundered, "positive_definite")))

  expect_true(evidence_satisfies(linop:::cape(bare, "positive_definite"),
                                 linop:::CG_PD_REQUIREMENT))
  expect_false(evidence_satisfies(linop:::cape(laundered, "positive_definite"),
                                  linop:::CG_PD_REQUIREMENT))

  set.seed(5)
  b <- matrix(stats::rnorm(n), n, 1)
  expect_equal(chosen(solve(bare, b, tol = 1e-11)), "cg")
  expect_false(chosen(solve(laundered, b, tol = 1e-11, maxit = 500L)) == "cg")
})

## ------------------------------------------------------- naming a method

test_that("every method of section 7.1 is reachable by name", {
  n <- 30
  spd <- as_spd_linop(kms_matrix(n, 0.4))
  set.seed(6)
  b <- matrix(stats::rnorm(n), n, 1)
  truth <- kms_inverse(n, 0.4) %*% b

  for (m in c("cg", "minres", "gmres", "fgmres", "bicgstab")) {
    x <- solve(spd, b, method = m, tol = 1e-11, maxit = 2000L)
    expect_equal(chosen(x), m, info = m)
    expect_match(why(x), "named by the caller", info = m)
    expect_equal(as.matrix(x), truth, tolerance = 1e-8, ignore_attr = TRUE,
                 info = m)
  }
  ## the least-squares pair on the same square system
  for (m in c("lsqr", "lsmr")) {
    x <- solve(spd, b, method = m, tol = 1e-11, maxit = 4000L)
    expect_equal(chosen(x), m, info = m)
    expect_equal(as.matrix(x), truth, tolerance = 1e-7, ignore_attr = TRUE,
                 info = m)
  }
})

test_that("naming a method is not filtered by evidence, and auto is", {
  ## The distinction the plan asks not to be collapsed. This operator declares
  ## positive definiteness and nothing establishes it beyond the declaration.
  n <- 20
  M <- spd_prescribed(n, seq(1, 3, length.out = n))
  A <- linop(M, properties = c(hermitian = TRUE, positive_definite = TRUE))
  set.seed(7)
  b <- matrix(stats::rnorm(n), n, 1)

  expect_equal(chosen(solve(A, b, method = "cg", tol = 1e-11)), "cg")
  expect_false(chosen(solve(A, b, tol = 1e-11)) == "cg")
})

test_that("a named method still refuses an operator that contradicts it", {
  ## Naming the method waives the evidence minimum, not the requirement.
  n <- 20
  A <- convdiff_1d(n, 0.3)
  expect_error(solve(A, rep(1, n), method = "cg"), "hermitian")
  expect_error(solve(A, rep(1, n), method = "minres"), "hermitian")
})

## -------------------------------------------------------------- the shapes

test_that("a rectangular operator is a different request, and the message says so", {
  f <- lsq_prescribed(30, 10, seq(1, 4, length.out = 10), seed = 8)
  A <- linop(f$A)
  set.seed(9)
  b <- matrix(stats::rnorm(30), 30, 1)

  expect_error(solve(A, b), "least-squares")
  expect_error(solve(A, b), "lsqr")
  expect_error(solve(A, b), "lsmr")

  for (m in c("lsqr", "lsmr")) {
    x <- solve(A, b, method = m, tol = 1e-12, maxit = 2000L)
    expect_equal(as.matrix(x), lsq_prescribed_solve(f, b), tolerance = 1e-8,
                 ignore_attr = TRUE, info = m)
  }
  ## and the square methods refuse it by name
  expect_error(solve(A, b, method = "cg"), "square")
  expect_error(solve(A, b, method = "gmres"), "square")
})

test_that("solve() without a right-hand side names what it will not form", {
  A <- laplacian_1d(5)
  expect_error(solve(A), "inverse")
})

## -------------------------------------------------------- what comes back

test_that("the default return is a plain vector or matrix with the certificate attached", {
  n <- 20
  A <- laplacian_1d(n)
  set.seed(10)
  b <- stats::rnorm(n)

  x <- solve(A, b, tol = 1e-11)
  expect_null(dim(x))
  expect_type(x, "double")
  expect_s3_class(attr(x, "certificate"), "linop_certificate")

  B <- matrix(stats::rnorm(n * 3), n, 3)
  X <- solve(A, B, tol = 1e-11)
  expect_equal(dim(X), c(n, 3L))
  expect_s3_class(attr(X, "certificate"), "linop_certificate")
})

test_that("indexing and coercion drop the certificate, and arithmetic keeps it", {
  ## Plan section 1.1 states the tradeoff the other way round, that arithmetic on
  ## the result drops the attribute. Measured, R's arithmetic operators copy
  ## attributes from their longer operand, so the certificate survives every one
  ## of them and it is subsetting and coercion that drop it. Both halves are
  ## asserted so the direction stays recorded: a certificate can outlive the
  ## value it describes, which is a sharper cost than the one the plan names.
  n <- 20
  A <- linop(spd_prescribed(n, seq(1, 4, length.out = n)))
  set.seed(30)
  b <- stats::rnorm(n)
  x <- solve(A, b, tol = 1e-11, maxit = 500L)

  for (e in list(x + 0, 2 * x, x + x, -x)) {
    expect_s3_class(attr(e, "certificate"), "linop_certificate")
  }
  for (e in list(x[1:3], as.numeric(x), c(x), sum(x))) {
    expect_null(attr(e, "certificate"))
  }

  X <- solve(A, cbind(b, b), tol = 1e-11, maxit = 500L)
  expect_s3_class(attr(X + 0, "certificate"), "linop_certificate")
  expect_null(attr(X[, 1], "certificate"))
})

test_that("details = TRUE returns the full object and the plain one is its solution", {
  n <- 20
  A <- linop(spd_prescribed(n, seq(1, 4, length.out = n)))
  set.seed(11)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- solve(A, b, tol = 1e-11, maxit = 500L, details = TRUE)
  plain <- solve(A, b, tol = 1e-11, maxit = 500L)
  expect_equal(fit$x, plain, ignore_attr = TRUE)
  expect_true(fit$converged)
  expect_equal(fit$requested, "auto")
  expect_equal(fit$method, "minres")
  expect_true(fit$iterations > 0L)
})

test_that("the certificate prints which route ran and why", {
  A <- linop(spd_prescribed(15, seq(1, 4, length.out = 15)))
  cert <- attr(solve(A, rep(1, 15), tol = 1e-10, maxit = 300L), "certificate")
  out <- utils::capture.output(print(cert))
  expect_match(out[1L], "solved by minres")
  expect_match(out[1L], "definiteness is not established")
})

test_that("a solve that does not converge comes back as a certificate, not an error", {
  n <- 60
  A <- convdiff_1d(n, 0.5)
  set.seed(12)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- solve(A, b, tol = 1e-14, maxit = 2L, details = TRUE)
  expect_false(fit$converged)
  expect_equal(cert_status(fit$certificate, "convergence"), "fail")
})

## ------------------------------------------------------------ control

test_that("control reaches the method and a name it does not take is an error", {
  n <- 30
  A <- convdiff_1d(n, 0.3)
  set.seed(13)
  b <- matrix(stats::rnorm(n), n, 1)

  ## restart is a gmres knob and changes the answer's cost
  short <- solve(A, b, method = "gmres", tol = 1e-10, maxit = 300L,
                 control = list(restart = 5L), details = TRUE)
  long <- solve(A, b, method = "gmres", tol = 1e-10, maxit = 300L,
                control = list(restart = n), details = TRUE)
  expect_gt(short$restarts, long$restarts)

  expect_error(solve(A, b, method = "cg", control = list(restart = 5L)),
               "does not take: restart")
  expect_error(solve(A, b, method = "gmres", control = list(nonsense = 1)),
               "does not take: nonsense")
  ## and control may not set what solve() already has an argument for
  expect_error(solve(A, b, method = "gmres", control = list(tol = 1e-4)),
               "does not take: tol")
})

test_that("a preconditioner reaches the method through the same verb", {
  n <- 40
  A <- laplacian_1d(n)
  d <- diag(as.matrix(A))
  P <- preconditioner(function(R) R / d, fixed = TRUE, hermitian = TRUE,
                      positive_definite = TRUE, side = "left")
  set.seed(14)
  b <- matrix(stats::rnorm(n), n, 1)

  fit <- solve(A, b, method = "cg", tol = 1e-11, preconditioner = P, details = TRUE)
  expect_equal(fit$method, "cg")
  expect_true(fit$converged)
  ## and the refusals of section 4.3 survive the extra layer
  flex <- preconditioner(function(R) R, fixed = FALSE, side = "right")
  expect_error(solve(A, b, preconditioner = flex), "fixed = TRUE")
})

test_that("the direct route refuses what it has nothing to do with", {
  A <- linop_scaling(c(1, 2, 3, 4))
  P <- preconditioner(function(R) R, fixed = TRUE, hermitian = TRUE,
                      positive_definite = TRUE)
  expect_error(solve(A, rep(1, 4), preconditioner = P), "takes no preconditioner")
  expect_error(solve(A, rep(1, 4), x0 = rep(0, 4)), "no starting iterate")
  expect_error(solve(A, rep(1, 4), control = list(restart = 3L)),
               "does not take: restart")
  ## naming an iterative method is how the preconditioner becomes the point
  expect_silent(solve(A, rep(1, 4), method = "cg", preconditioner = P))
})

test_that("a diagonal declaration the operator contradicts is caught, not divided by", {
  ## A zero on the diagonal of an operator declaring full_rank. It is the
  ## operator's own action that gives the declaration away, which is where the
  ## check has to be.
  ##
  ## Reached directly, because auto will not take the direct route on a bare
  ## declaration and no constructor produces a diagonal operator that lies about
  ## its rank. The guard is on the route rather than on the way in, so it holds
  ## whenever a constructor starts establishing diagonality some other way.
  A <- linop(function(X) X * c(1, 0, 3), dim = c(3L, 3L),
             properties = c(diagonal = TRUE, full_rank = TRUE))
  expect_error(linop:::diag_solve(A, rep(1, 3), tol = 1e-8, floor_const = 4,
                                  norm_control = list()),
               "contradicted by the operator's own action")
  ## and auto does not take that route at all on a bare declaration
  expect_false(linop:::auto_method(A)$method == "direct")
})

## ---------------------------------------------------------------- arguments

test_that("bad arguments are refused by name", {
  A <- laplacian_1d(10)
  b <- rep(1, 10)
  expect_error(solve(A, b, method = "nope"), "unknown method")
  expect_error(solve(A, b, method = c("cg", "gmres")), "single string")
  expect_error(solve(A, b, control = "restart"), "control must be a list")
  expect_error(solve(A, rep(1, 9)), "non-conformable")
})
