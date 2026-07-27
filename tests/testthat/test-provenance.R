## Section 5.11. Core carries provenance and never inspects it.

test_that("provenance is NULL by default and round-trips", {
  A <- linop(rmat(4, 3, seed = 1))
  expect_null(provenance(A))

  B <- set_provenance(A, "linop.hilbert", list(scheme = "finite_section", n = 128))
  expect_equal(provenance(B)$provider, "linop.hilbert")
  expect_equal(provenance(B)$payload$n, 128)
  expect_null(provenance(strip_provenance(B)))
})

test_that("core never inspects the payload", {
  ## an payload that would break anything trying to read it
  hostile <- structure(list(), class = "no_methods_at_all")
  A <- set_provenance(linop(rmat(3, 3, seed = 2)), "somebody", hostile)
  expect_silent(dim(A))
  expect_silent(as.matrix(A))
  expect_silent(t(A))
  expect_silent(A %*% rnorm(3))
  cert <- verify(A)
  expect_false(cert$overall == "fail")
})

test_that("provenance propagates through the algebra without being depended on", {
  A <- set_provenance(linop(rmat(4, 4, seed = 3)), "prov", list(tag = "x"))
  B <- linop(rmat(4, 4, seed = 4))

  expect_equal(provenance(t(A))$provider, "prov")
  expect_equal(provenance(adjoint(A))$provider, "prov")
  expect_equal(provenance(2 * A)$provider, "prov")
  expect_equal(provenance(A + B)$provider, "prov")
  expect_equal(provenance(A %*% B)$provider, "prov")
  ## an expression with no provenance anywhere stays NULL
  expect_null(provenance(B %*% B))
})

test_that("the provider generics error informatively when unregistered", {
  p <- list(provider = "nobody", payload = NULL)
  expect_error(provenance_lift(p, 1), "no provenance_lift\\(\\) method")
  expect_error(provenance_refine(p, 2), "no provenance_refine\\(\\) method")
  expect_error(provenance_original_residual(p, NULL), "no provenance_original_residual")
  expect_match(provenance_summary(p), "no summary method registered")
})

test_that("a provider can register methods and core will route to them", {
  registerS3method("provenance_summary", "demo_prov",
                   function(p, ...) sprintf("finite section at n = %d", p$payload$n),
                   envir = globalenv())
  p <- structure(list(provider = "demo", payload = list(n = 256)),
                 class = "demo_prov")
  expect_equal(provenance_summary(p), "finite section at n = 256")
})

test_that("set_provenance validates its arguments", {
  A <- linop(rmat(3, 3, seed = 5))
  expect_error(set_provenance(A, 42, list()), "single string")
  expect_error(set_provenance(matrix(1), "p", list()), "expects a linop")
})
