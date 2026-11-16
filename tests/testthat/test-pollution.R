## The pollution test. A Jacobi operator with
## period-2 couplings (1, 0.4) has two bands, [-1.4, -0.6] and [0.6, 1.4], and no
## eigenvalue at all; every finite section of it puts a point at 0 in the gap,
## with an eigenvector sitting at the cut. A defect at the origin creates a
## genuine eigenvalue in the gap. The certificate's exact truncation term is what
## tells the two apart: the invented point comes back labelled, not refused.

ssh <- function(defect = 0) {
  linop_jacobi(diagonal = defect, offdiagonal = c(0.4, 1), offdiagonal_limit = c(1, 0.4))
}

test_that("periodic limits give bands with a gap, from the Bloch matrices", {
  H <- ssh()
  B <- linop:::jacobi_bands(H)
  expect_equal(B, rbind(c(-1.4, -0.6), c(0.6, 1.4)), tolerance = 1e-12)
  expect_true(any(grepl("\\[-1.4, -0.6\\] u \\[0.6, 1.4\\]", capture.output(print(H)))))
  ## inside a band nothing decays; in the gap the transfer matrix gives a rate
  expect_true(is.na(decay_rate(H, 1)))
  expect_true(decay_rate(H, 0) > 0 && decay_rate(H, 0) < 1)
})

test_that("the periodic decay rate is what a dense eigenvector's tail does", {
  Hd <- ssh(0.3)
  M <- as.matrix(finite_section(Hd, 120))
  e <- eigen(M, symmetric = TRUE)
  i <- which.min(abs(e$values + 0.584))
  u <- abs(e$vectors[, i])
  rho <- decay_rate(Hd, e$values[i])
  ## two periods apart, far out on the right, the tail falls by rho^4; the
  ## section at n = 120 perturbs the eigenpair by about rho^120 = 2e-6, which
  ## is the agreement to expect (measured 4e-6 relative)
  j <- 121L + 60L
  expect_equal(u[j + 4L] / u[j], rho^4, tolerance = 1e-4)
})

test_that("the point a section invents in the gap certifies as fail, at every width", {
  H <- ssh()
  for (n in c(20L, 41L)) {
    fit <- eigs(finite_section(H, n), k = 1, sigma = 0.1)
    cert <- fit$certificate
    expect_lt(abs(fit$values), 1e-8)
    expect_equal(cert_status(cert, "truncation bound"), "fail")
    expect_equal(cert_status(cert, "isolation"), "fail")
    expect_equal(cert_status(cert, "forward error"), "fail")
    ## the tail is pinned to the cut: the same truncation term at both widths
    expect_equal(max(cert$values$truncation), 0.916515138991168, tolerance = 1e-8)
  }
})

test_that("a genuine eigenvalue in the gap certifies as pass and matches a dense section", {
  Hd <- ssh(0.3)
  fit <- eigs(Hd, k = 1, sigma = -0.58)
  M <- as.matrix(finite_section(Hd, 400))
  e <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  ref <- e[which.min(abs(e + 0.584))]
  expect_equal(fit$values, ref, tolerance = 1e-10)
  cert <- fit$certificate
  for (row in c("truncation bound", "isolation", "forward error")) {
    expect_equal(cert_status(cert, row), "pass", info = row)
  }
  ## a value in an interior gap gets no bracket: the bracket's free end is the
  ## variational principle, which holds only beyond every band
  expect_equal(cert_status(cert, "enclosure"), "not_checked")
})

test_that("eigs() on the operator with no eigenvalue stops widening and says why", {
  fit <- eigs(ssh(), k = 1, sigma = 0.1)
  expect_equal(fit$certificate$overall, "fail")
  expect_lte(length(fit$widths), 3L)
  expect_match(fit$certificate$dispatch$reason, "widening stopped shrinking the tail")
})
