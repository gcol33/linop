## The truncation, and the certificate that carries a statement about an operator
## no computation touched. Ground truth is S0.6's single-site model throughout,
## and the bound is checked against the true error rather than against itself.

sect <- function(v, n, ...) finite_section(linop_jacobi(diagonal = v, ...), n = n)
top <- function(F, k = 1, ...) eigs(F, k = k, which = "largest_algebraic", ...)

test_that("a section is a finite operator holding an infinite one", {
  H <- linop_jacobi(diagonal = 1)
  F <- finite_section(H, n = 6)
  expect_equal(dim(F), c(13L, 13L))
  expect_type(dim(F), "integer")
  expect_true(linop:::is_finite_dim(F))
  ## the child is the operator itself, so the relationship is structural
  expect_identical(linop:::node_children(F)[[1L]], H)
  expect_false(linop:::is_finite_dim(linop:::node_children(F)[[1L]]))
  expect_output(print(F), "section \\|j\\| <= 6")
  expect_output(print(F), "jacobi on l\\^2\\(Z\\), sigma_ess \\[-2, 2\\]")
  expect_output(ex <- explain(F), "operator 13 x 13")
  expect_equal(nrow(ex), 2L)
  expect_equal(ex$nrow, c(13, Inf))
})

test_that("the section is the tridiagonal it says it is", {
  F <- finite_section(linop_jacobi(diagonal = c(1.5, -0.5, 2), offdiagonal = c(2, 3, 4, 5)),
                      n = 5)
  M <- as.matrix(F)
  expect_true(isTRUE(all.equal(M, t(M))))
  expect_equal(diag(M), c(0, 0, 0, 0, 1.5, -0.5, 2, 0, 0, 0, 0))
  expect_equal(M[cbind(1:10, 2:11)], c(1, 1, 1, 2, 3, 4, 5, 1, 1, 1))
  ## every entry off the three diagonals is zero
  expect_equal(sum(abs(M[abs(row(M) - col(M)) > 1])), 0)
})

test_that("the apply is the matrix, in every mode and on a complex block", {
  F <- sect(1, n = 7)
  M <- as.matrix(F)
  X <- rmat(15, 3, seed = 4)
  for (mode in c("N", "T", "C", "R")) {
    expect_equal(linop:::linop_apply(F, X, mode), M %*% X, ignore_attr = TRUE)
  }
  Z <- zmat(15, 2, seed = 5)
  expect_equal(F %*% Z, M %*% Z, ignore_attr = TRUE)
  expect_equal(verify(F)$overall, "pass")
})

test_that("a section that cuts through the window is refused, with the reason", {
  H <- linop_jacobi(diagonal = c(1, 2, 3))
  expect_equal(H$args$radius, 2L)
  expect_error(finite_section(H, n = 2), "cuts through the window")
  expect_error(finite_section(H, n = 2), "free recurrence")
  expect_silent(finite_section(H, n = 3))
  ## and a section of a section is an ordinary submatrix, which says nothing
  ## further about the operator underneath
  expect_error(finite_section(finite_section(H, n = 5), n = 4),
               "Truncating a truncation")
  expect_error(finite_section(linop(diag(3)), n = 4), "linop_jacobi\\(\\)")
})

test_that("the bound holds against the true error over the whole S0.6 ladder", {
  for (v in c(0.5, 1, 2, 4)) {
    truth <- sqrt(v^2 + 4)
    for (n in c(5, 10, 20, 40, 80)) {
      fit <- top(sect(v, n))
      cv <- fit$certificate$values
      expect_lte(abs(fit$values - truth), cv$bound)
      ## and the bound is the two terms in quadrature plus the floor, exactly
      expect_equal(cv$bound, sqrt(cv$residual^2 + cv$truncation^2) + cv$floor)
    }
  }
})

test_that("the truncation term decays at exactly the closed-form rate", {
  ## eta / rho^n is constant in n while the tail is resolved, which is what says
  ## the rate is rho with a computable constant rather than merely O(rho^n).
  v <- 1
  rho <- (-abs(v) + sqrt(v^2 + 4)) / 2
  ratio <- vapply(c(10, 20, 30, 40), function(n) {
    top(sect(v, n))$certificate$values$truncation / rho^n
  }, numeric(1))
  expect_equal(max(ratio) / min(ratio), 1, tolerance = 1e-3)
  expect_equal(ratio[[1L]], 0.5847, tolerance = 1e-3)
})

test_that("the truncation term stops measuring truncation at the floor", {
  ## Past the point where the tail falls below what the finite eigensolve stores,
  ## eta plateaus: the analytic value at n = 80 is 1.1e-17 and the measured one is
  ## four orders above it. The bound still holds, because the floor covers exactly
  ## that gap, and the certificate says which regime it is in.
  v <- 1
  rho <- (-abs(v) + sqrt(v^2 + 4)) / 2
  fit <- top(sect(v, 80))
  cv <- fit$certificate$values
  expect_gt(cv$truncation / rho^80, 1e3)
  expect_lte(abs(fit$values - sqrt(5)), cv$bound)
  ## and the row says so where the tail has reached the floor
  fit40 <- top(sect(4, 40))
  cv40 <- fit40$certificate$values
  expect_lt(cv40$truncation, cv40$floor)
  expect_match(fit40$certificate$checks$detail[fit40$certificate$checks$check ==
                                                 "truncation bound"],
               "no longer measures truncation")
})

test_that("the residual identity the certificate rests on is exact", {
  ## ||(H - q) u~||^2 = ||(H_n - q) u||^2 + b^2 (u_n^2 + u_{-n}^2). Checked by
  ## embedding the computed pair in a wider section, which agrees with H on every
  ## index the extended vector can reach.
  n <- 12
  F <- sect(1, n)
  fit <- top(F)
  u <- as.numeric(fit$vectors)
  cv <- fit$certificate$values

  N <- 40
  G <- sect(1, N)
  ub <- numeric(2 * N + 1)
  ub[(N - n + 1):(N + n + 1)] <- u
  wide <- sqrt(sum((as.numeric(G %*% ub) - fit$values * ub)^2)) / sqrt(sum(ub^2))

  expect_equal(wide, sqrt(cv$residual^2 + cv$truncation^2), tolerance = 1e-12)
})

test_that("the eigenvector decays at the predicted rate well inside the section", {
  ## Near the cut the section imposes u_{n+1} = 0, so the tail there is the
  ## geometric solution plus its reflection and the ratio is not rho. Inside, it
  ## is, and the deviation falls with the section width.
  H <- linop_jacobi(diagonal = 1)
  dev <- vapply(c(20, 40), function(n) {
    fit <- top(finite_section(H, n))
    u <- as.numeric(fit$vectors)
    j <- (n + 1) + seq_len(n %/% 2)
    max(abs(u[j + 1] / u[j] - decay_rate(H, fit$values)))
  }, numeric(1))
  expect_lt(dev[[1L]], 1e-4)
  expect_lt(dev[[2L]], 1e-8)
})

test_that("Temple's bracket contains the truth and beats the linear bound", {
  for (n in c(5, 10, 20, 40)) {
    fit <- top(sect(1, n))
    cv <- fit$certificate$values
    half <- cv$bound^2 / cv$gap
    ## the variational end is exact in exact arithmetic and carries the floor in
    ## this one, which is the S0.6 correction reaching a third row
    expect_lte(sqrt(5), fit$values + half)
    expect_gte(sqrt(5), fit$values - cv$floor)
    ## quadratic against linear: four orders at n = 20
    expect_lt(half, cv$bound)
  }
  expect_match(top(sect(1, 20))$certificate$checks$detail[6L],
               "assuming no spectrum between the band edge")
})

test_that("the bracket covers the value adjacent to the band and no other", {
  ## A well deep enough to carry several eigenvalues above the band. Only the
  ## lowest of them has nothing that could be between it and the edge; for any
  ## other the near end of Temple's interval is a neighbouring eigenvalue that no
  ## residual locates.
  W <- linop_jacobi(diagonal = rep(2, 5))
  one <- top(finite_section(W, n = 60), k = 1)
  expect_equal(linop:::cert_status(one$certificate, "enclosure"), "pass")

  three <- top(finite_section(W, n = 60), k = 3)
  expect_equal(sum(three$certificate$values$gap > 0), 3L)
  expect_match(three$certificate$checks$detail[three$certificate$checks$check ==
                                                 "enclosure"],
               "for the 1 of 3 values adjacent to the band")
  ## the forward-error bound covers all three: it needs no gap at all
  expect_equal(linop:::cert_status(three$certificate, "forward error"), "pass")
  expect_true(all(three$certificate$values$bound < 1e-8))
})

test_that("a discretisation of continuous spectrum is refused, at every width", {
  ## V = 0 has no eigenvalues at all. Every finite-section value lies strictly
  ## inside the band and approaches the edge from below, and the eigensolver is
  ## not what notices: it converges cleanly and returns one.
  H0 <- linop_jacobi(diagonal = 0)
  for (n in c(10, 40, 100)) {
    fit <- top(finite_section(H0, n = n))
    expect_lt(fit$values, 2)
    expect_equal(linop:::cert_status(fit$certificate, "isolation"), "fail")
    expect_equal(fit$certificate$overall, "fail")
    expect_match(fit$certificate$checks$detail[fit$certificate$checks$check ==
                                                 "isolation"],
                 "have not separated from the continuum")
    ## and with nothing isolated there is no eigenvalue to bracket either
    expect_equal(linop:::cert_status(fit$certificate, "enclosure"), "not_checked")
  }
})

test_that("the certificate is the fifth shape and rolls up over it", {
  fit <- top(sect(1, 50))
  cert <- fit$certificate
  expect_s3_class(cert, "linop_certificate")
  expect_equal(cert$subject, "finite section")
  expect_equal(cert$checks$check,
               c("arithmetic floor", "finite residual", "orthogonality",
                 "truncation bound", "isolation", "enclosure", "convergence",
                 "forward error"))
  expect_equal(cert$overall, "pass")
  ## no residual against a right-hand side and no backward error, because a
  ## finite-section eigenvalue answers no equation with data in it
  expect_false(any(c("residual", "backward error") %in% cert$checks$check))
  ## the one line without a deterministic bound is the one with an unverified
  ## hypothesis under it
  expect_equal(cert$without_deterministic_bound, "enclosure")
})

test_that("the forward-error bound rests on nothing declared", {
  ## eigs() on an ordinary operator records the hermitian capability under the
  ## Weyl bound, so a bound resting on a bare declaration fails a requirement the
  ## declaration would have failed directly. Here the class is self-adjoint by
  ## construction and the residual identity is arithmetic, so the row depends on
  ## nothing: the first unconditional deterministic bound in the package.
  cert <- top(sect(1, 50))$certificate
  ev <- cert$evidence[["forward error"]]
  expect_equal(ev$source, "theorem")
  expect_equal(ev$guarantee, "deterministic_bound")
  expect_length(ev$depends_on, 0L)
  expect_true(evidence_satisfies(ev, requirement(guarantees = "deterministic_bound")))

  declared <- linop(as.matrix(sect(1, 20)), properties = c(hermitian = TRUE),
                    check = FALSE)
  dev <- eigs(declared, k = 1, which = "largest_algebraic")$certificate$evidence[["forward error"]]
  expect_false(evidence_satisfies(dev, requirement(sources = "theorem")))
})

test_that("the convergence row names the knob that is short", {
  ## An eigensolve that met its target inside a section too narrow to carry the
  ## tolerance is not a failed eigensolve, and saying "stopped early" would send
  ## the reader to maxit.
  fit <- top(sect(1, 10))
  cert <- fit$certificate
  expect_equal(linop:::cert_status(cert, "finite residual"), "pass")
  expect_equal(linop:::cert_status(cert, "convergence"), "fail")
  expect_match(cert$checks$detail[cert$checks$check == "convergence"],
               "the eigensolve met its target and the section did not")
})

test_that("a section gets a wider default subspace than any other operator", {
  ## Measured, not assumed: at v = 0.3 the ordinary default of 21 reaches 4.2e-7
  ## at n = 40 and 9.8e-3 at n = 80, a wider section giving the worse answer,
  ## because enlarging n densifies the cluster the value has to be told apart
  ## from. dev_notes/spikes/section-ncv-probe.R.
  expect_equal(linop:::default_ncv(sect(1, 8), 1L), 40L)
  expect_equal(linop:::default_ncv(linop(diag(50)), 1L), 21L)
  ## and it is a floor rather than a formula, so a large request still widens it
  expect_equal(linop:::default_ncv(sect(1, 200), 60L), 121L)

  narrow <- top(sect(0.3, 80), ncv = 21L)
  wide <- top(sect(0.3, 80))
  expect_gt(abs(narrow$values - sqrt(0.09 + 4)), 1e-4)
  expect_lt(abs(wide$values - sqrt(0.09 + 4)), 1e-10)
  ## the narrow run is wrong and says so; it is never a silent wrong answer
  expect_equal(narrow$nconv, 0L)
  expect_equal(narrow$certificate$overall, "fail")
})

test_that("a section is an ordinary operator everywhere else", {
  F <- sect(1, 8)
  b <- as.numeric(as.matrix(F) %*% rep(1, 17))
  x <- solve(F, b, method = "minres")
  expect_equal(as.numeric(x), rep(1, 17), tolerance = 1e-8)
  expect_equal(dim(t(F) %*% F), c(17L, 17L))
  expect_equal(linop_cost(F), 85)
  sv <- svds(F, k = 2)
  expect_equal(sv$certificate$subject, "svd")
})

test_that("a three-site potential recovers against its own reference", {
  H <- linop_jacobi(diagonal = c(1.5, -0.5, 2))
  ref <- top(finite_section(H, n = 200))$values
  for (n in c(5, 10, 20, 40)) {
    fit <- top(finite_section(H, n = n))
    expect_lte(abs(fit$values - ref), fit$certificate$values$bound)
  }
  expect_equal(top(finite_section(H, n = 20))$values, 2.848291705548,
               tolerance = 1e-11)
})
