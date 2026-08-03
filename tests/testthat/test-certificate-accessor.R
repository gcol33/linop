## certificate(): one reader over every place a verb has room to put one.

spd <- function() linop(diag(c(2, 3, 4)),
                        properties = c(hermitian = TRUE, positive_definite = TRUE))

test_that("it reads the certificate off every result that carries one", {
  A <- spd()
  b <- c(2, 3, 4)

  ## a solve puts it in an attribute, so the result stays what a matrix solve
  ## returns
  x <- solve(A, b)
  expect_s3_class(certificate(x), "linop_certificate")
  expect_identical(certificate(x), attr(x, "certificate"))
  expect_equal(certificate(x)$subject, "solve")

  ## details = TRUE returns the solve itself, where it is a field
  full <- solve(A, b, details = TRUE)
  expect_identical(certificate(full), full$certificate)

  ## the spectral verbs carry it in a field of a classed result
  e <- eigs(A, k = 1)
  expect_identical(certificate(e), e$certificate)
  expect_equal(certificate(e)$subject, "eigen")
  s <- svds(A, k = 1)
  expect_identical(certificate(s), s$certificate)

  ## verify() returns one directly, so the accessor is the identity on it
  v <- verify(A)
  expect_identical(certificate(v), v)

  ## and the fifth shape, off a section
  H <- linop_jacobi(diagonal = 1)
  f <- eigs(finite_section(H, 30), k = 1, which = "largest_algebraic")
  expect_equal(certificate(f)$subject, "finite section")
})

test_that("the attribute survives arithmetic and not coercion, and the error says which", {
  A <- spd()
  x <- solve(A, c(2, 3, 4))

  ## what plan section 1.1 has backwards: arithmetic keeps it, so a certificate
  ## can outlive the value it describes
  expect_s3_class(certificate(x + 0), "linop_certificate")
  expect_s3_class(certificate(2 * x), "linop_certificate")
  expect_identical(certificate(2 * x), certificate(x))

  ## coercion, indexing and reduction drop it
  expect_error(certificate(as.numeric(x)), "no certificate here")
  expect_error(certificate(x[1]), "no certificate here")
  expect_error(certificate(sum(x)), "no certificate here")
  expect_error(certificate(as.numeric(x)), "arithmetic on the solution")
  expect_error(certificate(as.numeric(x)), "details = TRUE")
})

test_that("an operator has no certificate, and the error names verify()", {
  A <- spd()
  expect_error(certificate(A), "a certificate belongs to a result")
  expect_error(certificate(A), "verify\\(A\\)")
  expect_error(certificate(linop_jacobi(diagonal = 1)), "belongs to a result")
  ## and nothing at all is the same refusal as a solution that lost its attribute
  expect_error(certificate(42), "no certificate here")
  expect_error(certificate(list(a = 1)), "no certificate here")
})

test_that("it is one exported name and the methods are registered, not exported", {
  exports <- getNamespaceExports("linop")
  expect_true("certificate" %in% exports)
  for (cls in c("linop", "linop_eigen", "linop_svd", "linop_certificate", "default")) {
    expect_false(paste0("certificate.", cls) %in% exports)
    expect_false(is.null(utils::getS3method("certificate", cls, optional = TRUE)))
  }
})
