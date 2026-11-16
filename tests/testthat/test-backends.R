## The backend registry and the PRIMME engine behind method = "primme". The
## other suites resolve auto to the reference (helper-backend.R); this one names
## its engines, and clears that option where it tests auto's choice.

with_priority <- function(priority, code) {
  old <- options(linop.backend_priority = priority)
  on.exit(options(old), add = TRUE)
  force(code)
}

with_default_priority <- function(code) with_priority(NULL, code)

with_backend <- function(name, family, code, ...) {
  key <- linop:::register_backend(name, family, ...)
  on.exit(rm(list = key, envir = linop:::.backends), add = TRUE)
  force(code)
}

triplet_backward <- function(A, d, u, v, anorm) {
  ru <- linop:::linop_apply(A, v, "N") - u * rep(d, each = nrow(u))
  rv <- linop:::linop_apply(A, u, "C") - v * rep(d, each = nrow(v))
  max(sqrt(colSums(Mod(ru)^2) + colSums(Mod(rv)^2)) / sqrt(2)) / anorm
}

rq_on <- function(A, X) {
  AX <- linop:::linop_apply(A, X, "N")
  Re(colSums(Conj(X) * AX)) / Re(colSums(Conj(X) * X))
}

expect_eigen_recovery <- function(A, fit, truth, tolerance = 1e-8) {
  expect_equal(fit$nconv, length(truth))
  expect_lt(max(abs(fit$values - truth)), tolerance)
  residual <- linop:::linop_apply(A, fit$vectors, "N") -
    sweep(fit$vectors, 2L, truth, "*")
  expect_lt(max(sqrt(colSums(Mod(residual)^2))), tolerance)
}

## ---------------------------------------------------------------- registry --

test_that("every engine is listed with its verdict, and the one that ran says so", {
  fit <- eigs(laplacian_1d(30), k = 2, method = "lanczos")
  d <- fit$certificate$dispatch
  expect_equal(d$requested, "lanczos")
  expect_equal(d$chosen, "lanczos")
  expect_setequal(d$candidates$backend, c("lanczos", "arnoldi", "primme", "rspectra"))
  expect_equal(d$candidates$status[d$candidates$backend == "lanczos"], "selected")
  expect_match(d$reason, "requested by name")

  s <- svds(linop(rmat(12, 8, seed = 1L)), k = 2, method = "golub-kahan")
  expect_setequal(s$certificate$dispatch$candidates$backend, c("golub-kahan", "primme", "rspectra"))
  expect_equal(s$method, "golub-kahan")
})

test_that("an engine whose package is absent is rejected by name, and auto passes over it", {
  A <- laplacian_1d(30)
  with_backend("ghost", "eigs", priority = 1000,
               provides_label = "nothing", package = "linop.ghost.absent",
               info_fn = "ghost_engine", run = function(...) stop("never run"), {
    fit <- with_default_priority(eigs(A, k = 2))
    cand <- fit$certificate$dispatch$candidates
    expect_equal(cand$backend[1L], "ghost")
    expect_equal(cand$status[1L], "rejected")
    expect_equal(cand$reason[1L], "package linop.ghost.absent is not installed")
    expect_false(fit$certificate$dispatch$chosen == "ghost")
    expect_match(fit$certificate$dispatch$reason, "ghost package linop.ghost.absent is not installed")
    expect_error(eigs(A, k = 2, method = "ghost"),
                 "method = 'ghost'\\) cannot run: package linop.ghost.absent is not installed")
  })
})

test_that("supports() returns a reason, and the reason is what the rejection reports", {
  A <- laplacian_1d(30)
  with_backend("picky", "eigs", priority = 1000, provides_label = "nothing",
               supports = function(A, request) "it only takes k = 1",
               run = function(...) stop("never run"), {
    fit <- with_default_priority(eigs(A, k = 2))
    cand <- fit$certificate$dispatch$candidates
    expect_equal(cand$reason[cand$backend == "picky"], "it only takes k = 1")
    expect_error(eigs(A, k = 2, method = "picky"), "it only takes k = 1")
  })
})

test_that("a method that names no engine lists the ones there are", {
  expect_error(eigs(laplacian_1d(10), k = 2, method = "davidson"),
               "method is one of 'auto', 'arnoldi', 'lanczos', 'primme', 'rspectra'")
})

## -------------------------------------------------------- with the engine --

test_that("an engine speaking another contract version is rejected, not attempted", {
  skip_if_not_installed("linop.primme")
  A <- laplacian_1d(30)
  with_backend("primme-old", "eigs", priority = 1000,
               provides_label = "nothing", package = "linop.primme",
               info_fn = "primme_engine", contract_version = "0.9",
               run = function(...) stop("never run"), {
    fit <- with_default_priority(eigs(A, k = 2))
    cand <- fit$certificate$dispatch$candidates
    expect_equal(cand$reason[cand$backend == "primme-old"],
                 "linop.primme speaks contract 1.1; this linop speaks 0.9")
  })
})

test_that("auto takes primme when it is installed, and the priority option overrides it", {
  skip_if_not_installed("linop.primme")
  ## a dense leaf computes its symmetry, which auto acts on
  A <- linop(as.matrix(laplacian_1d(40)))
  with_default_priority({
    expect_equal(eigs(A, k = 2)$method, "primme")
    expect_equal(svds(A, k = 2)$method, "primme")
  })
  with_priority(c(primme = -1), {
    expect_equal(eigs(A, k = 2)$method, "lanczos")
    expect_equal(svds(A, k = 2)$method, "golub-kahan")
  })
})

test_that("primme refuses a tolerance below machine epsilon, and auto takes the reference", {
  skip_if_not_installed("linop.primme")
  A <- linop(as.matrix(laplacian_1d(30)))
  for (tol in c(0, 1e-17)) {
    fit <- with_default_priority(eigs(A, k = 2, tol = tol, maxit = 40L))
    expect_equal(fit$method, "lanczos")
    cand <- fit$certificate$dispatch$candidates
    expect_match(cand$reason[cand$backend == "primme"], "below machine epsilon")
    expect_error(eigs(A, k = 2, tol = tol, method = "primme"), "below machine epsilon")
  }
})

test_that("two installed engines report candidates, refusals and priority overrides", {
  skip_if_not_installed("linop.primme")
  skip_if_not_installed("RSpectra")
  A <- linop(as.matrix(laplacian_1d(40)))
  check_dispatch <- function(fit, requested, chosen, priorities, rejected = character()) {
    dispatch <- fit$certificate$dispatch
    candidates <- dispatch$candidates
    expect_equal(dispatch$requested, requested)
    expect_equal(dispatch$chosen, chosen)
    expect_equal(candidates$backend, names(priorities))
    expect_equal(candidates$priority, unname(priorities))
    expected <- ifelse(candidates$backend %in% names(rejected), "rejected", "eligible")
    expected[candidates$backend == chosen] <- "selected"
    expect_equal(candidates$status, expected)
    reasons <- rep("eligible", nrow(candidates))
    reasons[match(names(rejected), candidates$backend)] <- unname(rejected)
    expect_equal(candidates$reason, reasons)
    expect_match(dispatch$reason, if (requested == "auto") "highest priority" else "requested by name")
  }
  default <- c(primme = 100, lanczos = 0, arnoldi = -10, rspectra = -20)
  with_default_priority({
    check_dispatch(eigs(A, k = 2), "auto", "primme", default)
    for (method in c("primme", "rspectra")) {
      fit <- eigs(A, k = 2, method = method)
      check_dispatch(fit, method, method, default)
      expect_eigen_recovery(A, fit, rev(laplacian_1d_eigenvalues(40))[1:2])
    }
    rejected <- c(primme = "it refuses tol = 0, below machine epsilon")
    check_dispatch(eigs(A, k = 2, tol = 0, maxit = 80L), "auto", "lanczos",
                   default, rejected)
    expect_error(eigs(A, k = 2, tol = 0, method = "primme"), rejected[[1L]], fixed = TRUE)
  })
  preferred <- c(rspectra = 200, primme = 100, lanczos = 0, arnoldi = -10)
  with_priority(c(rspectra = 200), {
    check_dispatch(eigs(A, k = 2), "auto", "rspectra", preferred)
    check_dispatch(eigs(A, k = 2, method = "primme"), "primme", "primme", preferred)
    rejected <- c(primme = "with a block of k = 2 it needs ncv >= 5")
    check_dispatch(eigs(A, k = 2, ncv = 4L, maxit = 3000L), "auto", "rspectra",
                   preferred, rejected)
    Z <- linop(hpd_prescribed(40, seq(5, 1, length.out = 40), seed = 8L))
    rejected <- c(rspectra = "does not take a complex operator")
    check_dispatch(eigs(Z, k = 2), "auto", "primme", preferred, rejected)
    expect_error(eigs(Z, k = 2, method = "rspectra"), rejected[[1L]], fixed = TRUE)
  })
  with_priority(c(primme = -100, rspectra = -200), {
    priorities <- c(lanczos = 0, arnoldi = -10, primme = -100, rspectra = -200)
    check_dispatch(eigs(A, k = 2), "auto", "lanczos", priorities)
    check_dispatch(eigs(A, k = 2, method = "rspectra"), "rspectra", "rspectra", priorities)
  })
})

test_that("svds with a starting vector is the reference's, since primme needs both sides", {
  skip_if_not_installed("linop.primme")
  A <- linop(rmat(12, 8, seed = 2L))
  fit <- with_default_priority(svds(A, k = 2, v0 = rep(1, 8)))
  expect_equal(fit$method, "golub-kahan")
  expect_error(svds(A, k = 2, v0 = rep(1, 8), method = "primme"), "only the right one")
})

## --------------------------------------------- agreement with the reference --

test_that("PRIMME's Arnoldi preset and the reference agree at every step count tried", {
  ## Both keep a fully orthogonalised basis of the Krylov space of the same
  ## starting vector, so the Ritz values after m applies are the same numbers.
  ## Measured at 1 to 40 steps on
  ## four fixtures and five seeds: worst 4.7e-16 relative, no drift with m.
  ## Like BiCGSTAB and unlike LSQR, LSMR or MINRES, there is no step count at
  ## which this stops holding, because there is no orthogonality to lose.
  skip_if_not_installed("linop.primme")
  fixtures <- list(
    as_spd_linop(spd_prescribed(60, exp(seq(log(1e6), 0, length.out = 60)), seed = 3L)),
    linop(hpd_prescribed(60, exp(seq(log(1e4), 0, length.out = 60)), seed = 4L),
          properties = c(hermitian = TRUE)))
  for (A in fixtures) {
    v <- linop:::eigen_start_vector(60L, A$dtype, NULL, 7L)
    for (m in c(1L, 2L, 4L, 8L, 16L, 32L)) {
      ref <- eigs(A, k = 1L, which = "largest_algebraic", method = "lanczos",
                  maxit = m, ncv = m + 1L, tol = 0, v0 = v)
      out <- linop.primme::primme_eigsh(
        function(X) linop:::linop_apply(A, X, "N"), 60L, 1L, which = "LA",
        tol = .Machine$double.eps, maxit = m, ncv = 80L, min_restart = 1L,
        block = 1L, v0 = v, complex = A$dtype == "complex", floor = 0,
        preset = "PRIMME_Arnoldi")
      expect_equal(out$stats$numMatvecs, m)
      expect_lt(abs(rq_on(A, out$vectors) - ref$values) / abs(ref$values), 1e-14)
    }
  }
})

test_that("both engines reach the closed-form spectrum at every end, real and complex", {
  skip_if_not_installed("linop.primme")
  n <- 80
  A <- laplacian_1d(n)
  ev <- laplacian_1d_eigenvalues(n)
  want <- list(largest_algebraic = rev(ev)[1:3], smallest_algebraic = ev[1:3],
               largest = rev(ev)[1:3], smallest = ev[1:3])
  for (w in names(want)) {
    for (m in c("lanczos", "primme")) {
      fit <- eigs(A, k = 3, which = w, method = m, maxit = 3000L)
      expect_equal(fit$nconv, 3L, info = paste(w, m))
      expect_lt(max(abs(fit$values - want[[w]])), 1e-9, label = paste(w, m))
      expect_eigen_recovery(A, fit, want[[w]])
      expect_equal(cert_status(fit$certificate, "backward error"), "pass")
      expect_equal(cert_status(fit$certificate, "forward error"), "pass")
      ## The residual row holds the pair to tol relative to |value|, which only
      ## primme's stopping test asks for. The reference locks on the backward
      ## error alone, and at the small end of this spectrum (1.5e-3 against
      ## ||A|| near 4) it locks at 7.7e-10 against 1e-10.
      if (m == "primme") {
        expect_equal(cert_status(fit$certificate, "residual"), "pass")
      }
    }
  }

  lam <- c(8, 7.5, 7, seq(3, 1, length.out = 37))
  Z <- linop(hpd_prescribed(40, lam, seed = 9L), properties = c(hermitian = TRUE))
  fp <- eigs(Z, k = 3, which = "largest_algebraic", method = "primme")
  fl <- eigs(Z, k = 3, which = "largest_algebraic", method = "lanczos")
  expect_true(is.complex(fp$vectors))
  expect_lt(max(abs(fp$values - lam[1:3])), 1e-10)
  expect_lt(max(abs(fp$values - fl$values)), 1e-10)
  expect_eigen_recovery(Z, fp, lam[1:3])
})

test_that("shift-invert runs through primme on the same transformed operator", {
  skip_if_not_installed("linop.primme")
  n <- 60
  A <- laplacian_1d(n)
  ev <- laplacian_1d_eigenvalues(n)
  sigma <- 1
  near <- ev[order(abs(ev - sigma))][1:3]
  fit <- eigs(A, k = 3, sigma = sigma, method = "primme")
  expect_equal(fit$method, "primme shift-invert")
  expect_eigen_recovery(A, fit, near)
  ## PRIMME's stopping test runs on (A - sigma I)^-1, which bounds the residual
  ## on A by tol ||A - sigma I||, not by tol |value|; measured 1.13e-10 against
  ## 1e-10 here. The backward error and the forward bound hold.
  expect_equal(cert_status(fit$certificate, "backward error"), "pass")
  expect_equal(cert_status(fit$certificate, "forward error"), "pass")
})

test_that("PRIMME signed targets and shift-invert recover paired Laplacian eigenvectors", {
  skip_if_not_installed("linop.primme")
  n <- 40L
  values <- laplacian_1d_eigenvalues(n) - 2.63
  targets <- list(largest_algebraic = rev(values)[1:3], smallest_algebraic = values[1:3],
                  largest = values[order(abs(values), decreasing = TRUE)][1:3],
                  smallest = values[order(abs(values))][1:3])
  for (dtype in c("double", "complex")) {
    phase <- if (dtype == "complex") exp(1i * seq_len(n) / 7) else rep(1, n)
    apply <- function(X) {
      Y <- Conj(phase) * X
      phase * (laplacian_1d_apply(Y) - 2.63 * Y)
    }
    A <- linop(apply, adjoint = apply, dim = c(n, n), dtype = dtype,
               properties = c(hermitian = TRUE))
    for (which in names(targets)) {
      fit <- eigs(A, k = 3, which = which, method = "primme", maxit = 3000L)
      expect_equal(is.complex(fit$vectors), dtype == "complex")
      expect_eigen_recovery(A, fit, targets[[which]])
    }
    sigma <- 0.3
    fit <- eigs(A, k = 3, sigma = sigma, method = "primme", maxit = 3000L)
    expect_equal(fit$method, "primme shift-invert")
    expect_eigen_recovery(A, fit, values[order(abs(values - sigma))][1:3])
  }
})

test_that("PRIMME warm starts recover closed-form eigenpairs in both dtypes", {
  skip_if_not_installed("linop.primme")
  n <- 40L
  values <- laplacian_1d_eigenvalues(n)
  vectors <- laplacian_1d_eigenvectors(n)
  for (dtype in c("double", "complex")) {
    phase <- if (dtype == "complex") exp(1i * seq_len(n) / 7) else rep(1, n)
    apply <- function(X) phase * laplacian_1d_apply(Conj(phase) * X)
    A <- linop(apply, adjoint = apply, dim = c(n, n), dtype = dtype,
               properties = c(hermitian = TRUE))
    v0 <- phase * (vectors[, n] + 0.05 * rowSums(vectors[, 1:5]))
    fit <- eigs(A, k = 3, which = "largest_algebraic", method = "primme",
                v0 = v0, maxit = 3000L)
    expect_eigen_recovery(A, fit, rev(values)[1:3])
  }
})

test_that("primme's smallest singular value near zero meets the residual row", {
  ## sigma_min = 1.67e-8 against ||A|| = 2: the residual row's floor is set by
  ## the lower bound on ||A||, which has to include PRIMME's own estimate for
  ## the converged triplet to clear it.
  skip_if_not_installed("linop.primme")
  n <- 60
  M <- diag(n)[c(2:n, 1), ] - diag(n)
  M[1, 1] <- M[1, 1] + 1e-6
  truth <- min(svd(M)$d)
  fit <- svds(linop(M), k = 1, which = "smallest", method = "primme")
  expect_equal(fit$nconv, 1L)
  expect_lt(fit$iterations, 2000L)
  expect_equal(fit$d, truth, tolerance = 1e-6)
  expect_lte(sqrt(sum((M %*% fit$v - fit$d * fit$u)^2)), 1e-14)
  expect_lte(sqrt(sum((crossprod(M, fit$u) - fit$d * fit$v)^2)), 1e-14)
})

test_that("singular triplets through primme, real and complex", {
  skip_if_not_installed("linop.primme")
  sigma <- c(20, 9, 5, seq(2, 0.4, length.out = 17))
  for (dt in c("double", "complex")) {
    f <- lsq_prescribed(30, 20, sigma, seed = 5L, dtype = dt)
    for (M in list(f$A, Conj(t(f$A)))) {
      A <- linop(M)
      for (w in c("largest", "smallest")) {
        fit <- svds(A, k = 3, which = w, method = "primme", maxit = 3000L)
        tru <- if (w == "largest") sigma[1:3] else rev(sigma)[1:3]
        expect_equal(fit$nconv, 3L)
        expect_lt(max(abs(fit$d - tru)), 1e-9)
        expect_lte(triplet_backward(A, tru, fit$u, fit$v, max(sigma)), 1e-9)
        expect_true(cert_status(fit$certificate, "residual") %in% c("pass", "qualified"))
        expect_equal(is.complex(fit$u), dt == "complex")
      }
    }
  }
})

test_that("a repeated eigenvalue comes back with every copy, through a block of k", {
  ## An 8-fold eigenvalue at order 30: a block of one returned 6 or 7 copies on
  ## every one of ten seeds, and a
  ## block of k all eight on every seed.
  skip_if_not_installed("linop.primme")
  n <- 30
  lambda <- c(rep(10, 8), seq(2, 0.5, length.out = n - 8L))
  for (dtype in c("double", "complex")) {
    for (seed in 1:3) {
      M <- if (dtype == "complex") hpd_prescribed(n, lambda, seed) else
        spd_prescribed(n, lambda, seed)
      A <- linop(M)
      fit <- eigs(A, k = 8, which = "largest_algebraic", method = "primme")
      expect_equal(fit$values, rep(10, 8), tolerance = 1e-9)
      expect_lt(max(Mod(Conj(t(fit$vectors)) %*% fit$vectors - diag(8))), 1e-10)
      expect_equal(fit$engine$params$maxBlockSize, 8L)
      expect_eigen_recovery(A, fit, rep(10, 8))
    }
  }
})

test_that("a basis too small for a block of k is primme's to refuse, with the size it needs", {
  skip_if_not_installed("linop.primme")
  A <- laplacian_1d(200)
  expect_error(eigs(A, k = 4, ncv = 8L, method = "primme"), "needs ncv >= 9")
  expect_error(eigs(A, k = 1, ncv = 3L, which = "largest", method = "primme"),
               "needs ncv >= 4")
  expect_equal(eigs(A, k = 1, ncv = 3L, which = "largest_algebraic",
                    method = "primme", maxit = 3000L)$method, "primme")
  fit <- with_default_priority(eigs(linop(as.matrix(A)), k = 4, ncv = 8L))
  expect_equal(fit$method, "lanczos")
})

## ------------------------------------------------- the certificate survives --

test_that("primme's values are Rayleigh quotients measured on A, under a shift too", {
  ## Under shift-invert the engine iterates on (A - sigma I)^-1, so its own value
  ## for a pair is sigma + 1/theta and carries the inner solve's error. What is
  ## reported is x^H A x / x^H x on A, which is exactly what is checked here.
  skip_if_not_installed("linop.primme")
  A <- laplacian_1d(50)
  for (fit in list(eigs(A, k = 3, which = "largest_algebraic", method = "primme"),
                   eigs(A, k = 2, sigma = 0.5, method = "primme",
                        inner = list(tol = 1e-6)))) {
    expect_equal(fit$values, rq_on(A, fit$vectors), tolerance = 1e-15)
  }
  ## and svds reports Re(u^H A v) of unit vectors
  B <- linop(rmat(15, 10, seed = 3L))
  s <- svds(B, k = 2, method = "primme")
  expect_equal(s$d, Re(colSums(Conj(s$u) * (as.matrix(B) %*% s$v))), tolerance = 1e-15)
})

test_that("the forward bound rests on the hermitian declaration whichever engine ran", {
  skip_if_not_installed("linop.primme")
  M <- spd_prescribed(20, seq(10, 1, length.out = 20), seed = 5L)
  strict <- requirement(sources = c("theorem", "construction", "computation"),
                        guarantees = c("identity", "deterministic_bound"))
  computed <- eigs(linop(M), k = 2, which = "largest_algebraic", method = "primme")
  expect_true(evidence_satisfies(computed$certificate$evidence[["forward error"]], strict))
  declared <- eigs(linop(function(X) M %*% X, dim = c(20, 20),
                         properties = c(hermitian = TRUE)),
                   k = 2, which = "largest_algebraic", method = "primme")
  ev <- declared$certificate$evidence[["forward error"]]
  expect_false(evidence_satisfies(ev, strict))
  expect_true(evidence_satisfies(ev, requirement(sources = c("theorem", "user_declaration"))))
  row <- declared$certificate$checks
  expect_equal(row$guarantee[row$check == "forward error"], "deterministic_bound")
})

test_that("primme's convergence test is the certificate's residual row too", {
  ## PRIMME's own test is the backward error alone, ||r|| <= tol ||A||, and a
  ## pair stopping exactly there fails the residual row whenever |value| is
  ## below ||A||: on this well, values near -0.9 against ||A|| near 4 read
  ## 4.0e-10 against 1e-10. The engine's test holds the pair to both.
  skip_if_not_installed("linop.primme")
  n <- 2000
  v <- numeric(n); v[(n %/% 2 - 14):(n %/% 2 + 15)] <- 1
  ap <- function(X) laplacian_1d_apply(X) - v * X
  A <- linop(ap, adjoint = ap, dim = c(n, n), properties = c(hermitian = TRUE))
  fit <- eigs(A, k = 4, which = "smallest_algebraic", tol = 1e-10, method = "primme",
              maxit = 5000L)
  expect_equal(cert_status(fit$certificate, "residual"), "pass")
  expect_equal(fit$nconv, 4L)
})

test_that("an exhausted budget comes back as a fail certificate, even with no pair", {
  skip_if_not_installed("linop.primme")
  fit <- eigs(laplacian_1d(400), k = 4, which = "smallest_algebraic",
              method = "primme", maxit = 30L)
  expect_equal(cert_status(fit$certificate, "convergence"), "fail")
  expect_match(fit$certificate$checks$detail[fit$certificate$checks$check == "convergence"],
               "budget exhausted")

  ## PRIMME's singular value solver hands back only what it converged, so a
  ## budget that converges nothing returns nothing, and says so.
  s <- svds(diff_1d(80), k = 4, maxit = 6L, method = "primme")
  expect_length(s$d, 0L)
  expect_equal(s$nconv, 0L)
  expect_equal(s$certificate$overall, "fail")
  expect_equal(cert_status(s$certificate, "residual"), "not_checked")
  expect_match(s$certificate$checks$detail[s$certificate$checks$check == "convergence"],
               "0 of 4 requested pairs")
})

test_that("the priority option is left as the suite set it", {
  expect_identical(getOption("linop.backend_priority"), c(primme = -Inf))
})

test_that("an error inside the operator surfaces through the engine as an R error", {
  skip_if_not_installed("linop.primme")
  bad <- linop(function(X) stop("the operator failed"), dim = c(30, 30),
               properties = c(hermitian = TRUE))
  expect_error(eigs(bad, k = 2, method = "primme"), "the operator failed")
})

test_that("auto acts on a hermitian declaration only at the evidence minimum", {
  ## A bare user declaration is honoured by name and passed over by auto, the
  ## split solve() makes between method = "minres" and auto's MINRES.
  A <- laplacian_1d(30)
  fit <- eigs(A, k = 2, which = "largest_algebraic")
  expect_equal(fit$method, "arnoldi")
  cand <- fit$certificate$dispatch$candidates
  expect_match(cand$reason[cand$backend == "lanczos"],
               "auto acts on hermitian only at computed, constructed or proven evidence")
  expect_false(grepl("\n", cand$reason[cand$backend == "lanczos"]))
  expect_equal(eigs(A, k = 2, which = "largest_algebraic", method = "lanczos")$method,
               "lanczos")
  ## the certificate still reads the declaration, recorded under the forward row
  expect_equal(cert_status(fit$certificate, "forward error"), "pass")
  ## a dense leaf computes its symmetry, and auto takes Lanczos
  expect_equal(eigs(linop(as.matrix(A)), k = 2, which = "largest_algebraic")$method,
               "lanczos")
})
