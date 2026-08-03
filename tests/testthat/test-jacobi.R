## The first operator whose index set is infinite, and the closed forms S0.6
## derived for it. Ground truth throughout is the single-site model
##
##     V = v delta_0,   lambda = sign(v) sqrt(v^2 + 4),
##     rho = (-|v| + sqrt(v^2 + 4)) / 2
##
## which is exact, so these are parameter-recovery tests rather than agreement
## against another implementation.

single_site <- function(v) list(op = linop_jacobi(diagonal = v),
                                lambda = sign(v) * sqrt(v^2 + 4),
                                rho = (-abs(v) + sqrt(v^2 + 4)) / 2)

test_that("a Jacobi operator lives on a sequence space and says so", {
  H <- linop_jacobi(diagonal = 1)
  expect_equal(dim(H), c(Inf, Inf))
  expect_false(linop:::is_finite_dim(H))
  expect_true(is.na(length(H)))
  expect_equal(linop_cost(H), Inf)
})

test_that("a constructor states construction, not a declaration", {
  ## The hilbert spike found the matrix-free route reporting user_declaration for
  ## a property the constructor knows the way linop_eye() knows the identity is
  ## unitary. Inside one package there is a constructor to say it, and the
  ## forward-error bound the certificate carries rests on this evidence.
  H <- linop_jacobi(diagonal = c(1, -2, 3))
  for (nm in c("hermitian", "symmetric", "real")) {
    expect_true(isTRUE(cap(H, nm)$value))
    expect_equal(cap(H, nm)$evidence$source, "construction")
    expect_equal(cap(H, nm)$evidence$guarantee, "identity")
    expect_length(cap(H, nm)$evidence$depends_on, 0L)
  }
})

test_that("nothing numeric runs on it, and every refusal names the way out", {
  H <- linop_jacobi(diagonal = 1)
  expect_error(as.matrix(H), "finite_section\\(\\)")
  expect_error(svds(H, k = 1), "svds\\(\\) needs an operator with finite")
  expect_error(verify(H), "verify\\(\\) needs an operator with finite")
  expect_error(solve(H, matrix(1, 3, 1)), "solve\\(\\) needs an operator with finite")
  expect_error(H %*% matrix(1, 3, 1), "apply\\(\\) needs an operator with finite")
  ## and reaching the handler past the gate gets the same sentence
  expect_error(linop:::jacobi_apply(H, matrix(1, 3, 1), "N"),
               "no block to act on")
  ## eigs() is the one verb that is not refused, because the node type registers
  ## a way to make a matrix; test-sequence-eigs.R is that route
  expect_s3_class(eigs(H, k = 1, which = "largest_algebraic"), "linop_eigen")
})

test_that("the essential spectrum is read off the limiting coefficients", {
  ## Weyl: a finitely supported perturbation is finite rank, hence compact, so the
  ## essential spectrum is the free one and the symbol gives it exactly.
  expect_equal(linop:::jacobi_band(linop_jacobi(diagonal = 1)), c(-2, 2))
  expect_equal(linop:::jacobi_band(linop_jacobi(diagonal = 1, offdiagonal = 3)),
               c(-6, 6))
  H <- linop_jacobi(diagonal = 1, offdiagonal = 3,
                    diagonal_limit = 5, offdiagonal_limit = 2)
  expect_equal(linop:::jacobi_band(H), c(1, 9))
  ## the sign of the coupling does not move the band, since the symbol is
  ## a + 2 b cos(theta) over a full period
  expect_equal(linop:::jacobi_band(linop_jacobi(diagonal = 1, offdiagonal = -1)),
               c(-2, 2))
})

test_that("the decay rate is the closed form, over the whole model", {
  for (v in c(0.5, 1, 2, 4, -1, -3)) {
    m <- single_site(v)
    expect_equal(decay_rate(m$op, m$lambda), m$rho, tolerance = 1e-14)
  }
  ## v = 1 is the golden-ratio conjugate, which is the case worth naming
  expect_equal(decay_rate(linop_jacobi(diagonal = 1), sqrt(5)),
               (sqrt(5) - 1) / 2, tolerance = 1e-15)
})

test_that("the decay rate is NA where nothing decays", {
  H <- linop_jacobi(diagonal = 1)
  ## inside the essential spectrum both roots are on the unit circle
  expect_true(is.na(decay_rate(H, 0)))
  expect_true(is.na(decay_rate(H, 1.9)))
  ## and at the edge itself they coincide at 1
  expect_true(is.na(decay_rate(H, 2)))
  expect_true(is.na(decay_rate(H, -2)))
  expect_equal(is.na(decay_rate(H, c(0, 3, 1, -3))), c(TRUE, FALSE, TRUE, FALSE))
})

test_that("the decay rate does not cancel for a well separated eigenvalue", {
  ## Written as |mu|/2 - sqrt((mu/2)^2 - 1) this is the difference of two nearly
  ## equal numbers and reaches exactly zero well before the true value does. The
  ## two roots are reciprocal, so rho + 1/rho = |mu| is the identity to check.
  H <- linop_jacobi(diagonal = 1)
  for (lambda in c(10, 1e4, 1e8, 1e12)) {
    rho <- decay_rate(H, lambda)
    expect_gt(rho, 0)
    expect_equal(rho + 1 / rho, lambda, tolerance = 1e-12)
  }
  naive <- function(mu) abs(mu) / 2 - sqrt((mu / 2)^2 - 1)
  expect_equal(naive(1e12), 0)
  expect_gt(decay_rate(H, 1e12), 0)
})

test_that("the coefficients sit where the window says they do", {
  H <- linop_jacobi(diagonal = c(7, 8, 9), from = -1)
  expect_equal(linop:::jacobi_a(H, -3:3), c(0, 0, 7, 8, 9, 0, 0))
  ## the default centres the window on the origin
  expect_equal(linop_jacobi(diagonal = c(7, 8, 9))$args$first, -1L)
  expect_equal(linop_jacobi(diagonal = 5)$args$first, 0L)

  ## b_j couples j and j+1, so the couplings touching a window of three diagonal
  ## entries are the one on each edge and the two between them
  G <- linop_jacobi(diagonal = c(7, 8, 9), offdiagonal = c(2, 3, 4, 5), from = -1)
  expect_equal(linop:::jacobi_b(G, -3:3), c(1, 2, 3, 4, 5, 1, 1))
})

test_that("the norm bound is the largest row sum, closed form", {
  expect_equal(linop:::jacobi_norm_bound(linop_jacobi(diagonal = 1)), 3)
  expect_equal(linop:::jacobi_norm_bound(linop_jacobi(diagonal = 4)), 6)
  ## and it is an upper bound on the section's spectral radius, which is what the
  ## arithmetic floor needs of it
  for (v in c(0.5, 4)) {
    M <- as.matrix(finite_section(linop_jacobi(diagonal = v), n = 30))
    expect_lte(max(abs(eigen(M, symmetric = TRUE, only.values = TRUE)$values)),
               linop:::jacobi_norm_bound(linop_jacobi(diagonal = v)))
  }
})

test_that("the constructor refuses what the class is not", {
  expect_error(linop_jacobi(diagonal = 1i), "must be real")
  expect_error(linop_jacobi(diagonal = 1, offdiagonal = 1i), "must be real")
  expect_error(linop_jacobi(diagonal = numeric(0)), "at least one entry")
  expect_error(linop_jacobi(diagonal = NA_real_), "at least one entry")
  expect_error(linop_jacobi(diagonal = c(1, 2), offdiagonal = c(1, 2)),
               "offdiagonal is one number, or 3 of them")
  expect_error(linop_jacobi(diagonal = 1, offdiagonal_limit = 0),
               "essential spectrum is a set of points")
  expect_error(linop_jacobi(diagonal = 1, diagonal_limit = Inf), "finite real")
})

test_that("a constant off-diagonal sets its own limit and a windowed one does not", {
  expect_equal(linop_jacobi(diagonal = 1, offdiagonal = 3)$args$b_inf, 3)
  expect_equal(linop_jacobi(diagonal = 1, offdiagonal = c(3, 3))$args$b_inf, 1)
  expect_equal(linop_jacobi(diagonal = 1, offdiagonal = c(3, 3),
                            offdiagonal_limit = 2)$args$b_inf, 2)
})

test_that("decay_rate refuses an operator that carries no limiting coefficients", {
  expect_error(decay_rate(linop(diag(3)), 1), "linop_jacobi\\(\\)")
  expect_error(decay_rate(3, 1), "linop_jacobi\\(\\)")
  ## but it reads through a section, since the section holds the operator
  H <- linop_jacobi(diagonal = 1)
  expect_equal(decay_rate(finite_section(H, n = 8), sqrt(5)), decay_rate(H, sqrt(5)))
})

test_that("the general free part is the single-site model shifted and scaled", {
  ## H = c I + s H0 for the standard H0 with potential v, so every eigenvalue is
  ## c + s lambda_0 and the section recovers it. This is the whole content of the
  ## limits being arguments rather than fixed at 0 and 1.
  cc <- 1.5; s <- 0.75; v <- 2
  H <- linop_jacobi(diagonal = cc + s * v, offdiagonal = s,
                    diagonal_limit = cc, offdiagonal_limit = s)
  expect_equal(linop:::jacobi_band(H), cc + s * c(-2, 2))
  fit <- eigs(finite_section(H, n = 40), k = 1, which = "largest_algebraic")
  expect_equal(fit$values, cc + s * sqrt(v^2 + 4), tolerance = 1e-12)
  ## and the decay rate is the same as the unscaled model's, since it depends on
  ## the eigenvalue only through (lambda - a_inf) / b_inf
  expect_equal(decay_rate(H, cc + s * sqrt(v^2 + 4)),
               decay_rate(linop_jacobi(diagonal = v), sqrt(v^2 + 4)),
               tolerance = 1e-14)
})
