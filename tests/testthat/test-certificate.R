## Section 6. build_certificate() and cert_rows() are the provider surface: a
## package outside this one computes its own quantities and reports them in this
## object rather than in a second one of its own. So the contract these assert is
## the one an outside caller sees -- the roll-up rule, what the vocabulary refuses,
## and that a row's evidence survives into evidence_satisfies().

test_that("overall is a roll-up of the four statuses and nothing else", {
  mk <- function(...) {
    statuses <- c(...)
    r <- cert_rows()
    for (i in seq_along(statuses)) r$add(paste0("check ", i), statuses[i])
    build_certificate(r$collect(), subject = "test")
  }
  expect_equal(mk("pass", "pass")$overall, "pass")
  expect_equal(mk("pass", "qualified")$overall, "qualified")
  expect_equal(mk("pass", "not_checked")$overall, "qualified")
  expect_equal(mk("pass", "fail")$overall, "fail")
  ## a failure outranks a qualification, whichever order they arrive in
  expect_equal(mk("qualified", "fail")$overall, "fail")
  expect_equal(mk("fail", "qualified")$overall, "fail")
  expect_equal(mk("not_checked", "fail")$overall, "fail")
})

test_that("without_deterministic_bound names the rows with no bound under them", {
  r <- cert_rows()
  r$add("exact", "pass", source = "computation", guarantee = "identity")
  r$add("bounded", "pass", source = "theorem", guarantee = "deterministic_bound")
  r$add("estimated", "pass", source = "computation", guarantee = "estimate")
  r$add("probabilistic", "pass", source = "probe", guarantee = "probabilistic_bound",
        confidence = 0.99)
  cert <- build_certificate(r$collect(), subject = "test")
  expect_equal(cert$without_deterministic_bound, c("estimated", "probabilistic"))
})

test_that("a not_checked row is listed even when its guarantee is deterministic", {
  ## The latent bug the eigenpair certificate exposed, from the other side: a row
  ## that was never made carries no bound whatever its guarantee field says, and a
  ## deterministic_bound row that WAS made must not be listed beside it.
  r <- cert_rows()
  r$add("never made", "not_checked", source = "computation", guarantee = "identity")
  r$add("weyl", "pass", source = "theorem", guarantee = "deterministic_bound")
  cert <- build_certificate(r$collect(), subject = "test")
  expect_equal(cert$without_deterministic_bound, "never made")
})

test_that("the status vocabulary is closed, at both ends", {
  r <- cert_rows()
  expect_error(r$add("x", "warn"), "status must be one of")
  expect_error(r$add("x", "ok"), "status must be one of")
  expect_error(r$add("x", c("pass", "fail")), "status must be one of")

  ## and again where a frame arrives from somewhere other than cert_rows()
  df <- data.frame(check = "x", status = "warn", source = "computation",
                   guarantee = "identity", confidence = 1, detail = "",
                   stringsAsFactors = FALSE)
  expect_error(build_certificate(df, subject = "test"), "unknown status")
  expect_error(build_certificate(df, subject = "test"), "counted as a pass")
})

test_that("the evidence vocabulary is the same one evidence() uses", {
  r <- cert_rows()
  expect_error(r$add("x", "pass", source = "vibes"), "source must be one of")
  expect_error(r$add("x", "pass", guarantee = "probably"), "guarantee must be one of")
  expect_error(r$add("x", "pass", evidence = list(source = "theorem")),
               "must come from evidence")
})

test_that("build_certificate refuses a frame that is not a certificate table", {
  r <- cert_rows()
  r$add("x", "pass")
  df <- r$collect()
  expect_error(build_certificate(df[, setdiff(names(df), "detail")], subject = "t"),
               "missing the column detail")
  expect_error(build_certificate(df[, c("check", "status")], subject = "t"),
               "missing the columns")
  expect_error(build_certificate(df[0, ], subject = "t"), "at least one row")
  expect_error(build_certificate(list(check = "x"), subject = "t"), "data frame")
  expect_error(build_certificate(df, subject = c("a", "b")), "single string")
})

test_that("a row's evidence fills the flat fields and survives into the object", {
  ev <- evidence("theorem", "deterministic_bound", 1,
                 depends_on = list(evidence("user_declaration", "identity", 1)))
  r <- cert_rows()
  r$add("forward error", "pass", "Weyl", evidence = ev)
  df <- r$collect()
  ## the table shows the flat fields, read off the object rather than typed twice
  expect_equal(df$source, "theorem")
  expect_equal(df$guarantee, "deterministic_bound")
  expect_equal(df$confidence, 1)

  cert <- build_certificate(df, subject = "test", evidence = r$collect_evidence())
  expect_s3_class(cert$evidence[["forward error"]], "linop_evidence")

  ## and the laundering case reaches a provider's certificate the way it reaches a
  ## propagated capability: the bound is a theorem, what it rests on is not.
  req <- requirement(sources = c("construction", "computation", "theorem"))
  expect_false(evidence_satisfies(cert$evidence[["forward error"]], req))
  expect_true(evidence_satisfies(ev$depends_on[[1L]],
                                 requirement(sources = "user_declaration")))
})

test_that("collect_evidence is NULL when no row carried an object", {
  r <- cert_rows()
  r$add("x", "pass")
  expect_null(r$collect_evidence())
  cert <- build_certificate(r$collect(), subject = "test",
                            evidence = r$collect_evidence())
  expect_null(cert$evidence)
})

test_that("rows come back in the order they were added", {
  r <- cert_rows()
  for (nm in c("third", "first", "second")) r$add(nm, "pass")
  expect_equal(r$collect()$check, c("third", "first", "second"))
})

test_that("a provider builds and prints a certificate of its own shape", {
  ## The three-part finite-section certificate of S0.6, built with exported names
  ## only. Its rows are not core's rows: there is no residual line and no
  ## backward-error line, because a finite-section eigenvalue is not the answer to
  ## an equation with a right-hand side.
  eta <- 3.863955e-05
  gap <- 0.236067977
  floor_v <- 4 * 3 * .Machine$double.eps
  herm <- evidence("user_declaration", "identity", 1)

  r <- cert_rows()
  r$add("arithmetic floor", "pass",
        sprintf("||H|| <= 2 + max|V| = 3; c = 4, floor %.3e", floor_v),
        source = "construction", guarantee = "identity")
  r$add("truncation bound", "pass",
        sprintf("dist(theta, spec H) <= %.3e", eta + floor_v),
        evidence = evidence("theorem", "deterministic_bound", 1,
                            depends_on = list(herm)))
  r$add("isolation gap", "pass", sprintf("q - a = %.3e against the band edge 2", gap),
        source = "computation", guarantee = "identity")
  cert <- build_certificate(r$collect(), subject = "finite section",
                            values = list(eta = eta, gap = gap),
                            evidence = r$collect_evidence())

  expect_s3_class(cert, "linop_certificate")
  expect_equal(cert$overall, "pass")
  expect_equal(cert$subject, "finite section")
  expect_equal(cert$checks$check,
               c("arithmetic floor", "truncation bound", "isolation gap"))
  expect_equal(cert$without_deterministic_bound, character(0))
  expect_equal(cert$values$eta, eta)
  ## the printed form needs neither a dispatch record nor probes
  expect_output(print(cert), "finite section|truncation bound")
  expect_output(print(cert), "overall\\s+pass")
})

test_that("a refused finite section reports as a failure with its reason", {
  ## V = 0: the top Ritz value has not separated from the essential spectrum, so
  ## q - a > 0 fails and no eigenvalue may be claimed.
  r <- cert_rows()
  r$add("arithmetic floor", "pass", "||H|| <= 2", source = "construction")
  r$add("truncation bound", "not_checked",
        "no eigenvalue to bound: the Ritz value is inside the essential spectrum",
        source = "computation", guarantee = "identity")
  r$add("isolation gap", "fail", "q - a = -2.036e-02, the Ritz value is inside [-2, 2]",
        source = "computation", guarantee = "identity")
  cert <- build_certificate(r$collect(), subject = "finite section")

  expect_equal(cert$overall, "fail")
  expect_equal(cert$without_deterministic_bound, "truncation bound")
  expect_output(print(cert), "failures:")
  expect_output(print(cert), "inside \\[-2, 2\\]")
})
