## What each method costs and what it comes back with, on operators whose answer
## is known before the solve runs.
##
## Every method on a row's operator gets the same tolerance, the same budget, the
## same right-hand side and the same zero start. A comparison in which two
## methods ran to different maxit is a measurement of the budget, and the MINRES
## suite already shipped one draft that made exactly that mistake.
##
## The roster per fixture is not a ranking. It is what the operator supplies:
## convdiff_1d has no adjoint, so five of the seven cannot run on it at all, and
## CG is absent wherever the operator is indefinite because running it there
## contradicts a declaration rather than converging slowly.
##
## Each fixture also carries a ||A|| row. The certificate's arithmetic floor
## rests on that estimate, every solve row above it paid for one, and the counted
## operators are callback leaves, so norm2() takes the power-iteration route on
## all of them -- the route a matrix-free caller gets. Subtracting the ||A|| row
## from a solve row leaves what the iteration itself spent.
##
## A right-hand side never draws from the seed its fixture drew from. Both
## generators here build their operator out of a QR factorisation of a random
## matrix, and the first n draws from a seed are the first column of that matrix,
## which the factorisation puts exactly in the span of its own first column: with
## a shared seed b is a multiple of the first eigenvector of spd_prescribed, and
## it is exactly in the range of lsq_prescribed, so the incompatible fixture
## becomes compatible and CG converges in one step. The rectangular case measures
## and reports how far outside the range its b actually falls.

SOLVER_TOL <- 1e-8

solver_cases <- function() {
  cases <- list()
  add <- function(...) cases[[length(cases) + 1L]] <<- list(...)

  ## Matrix-free, hermitian positive definite, spectrum in closed form.
  add(label = "laplacian_1d(200)",
      methods = c("cg", "minres", "gmres", "bicgstab", "lsqr", "lsmr"),
      truth = "closed-form eigendecomposition",
      build = function(cn) {
        n <- 200L
        lam <- laplacian_1d_eigenvalues(n)
        set.seed(101)
        b <- stats::rnorm(n)
        list(A = counted_linop(laplacian_1d_apply, laplacian_1d_apply, c(n, n), cn,
                               properties = c(hermitian = TRUE,
                                              positive_definite = TRUE)),
             b = b, truth = shifted_laplacian_solve(n, 0, b),
             kappa = max(lam) / min(lam))
      })

  ## The same operator shifted into the spectrum: hermitian, indefinite, still
  ## closed form in both its eigenvalues and its eigenvectors. CG is not on the
  ## roster because p^H A p <= 0 here is a contradicted declaration.
  add(label = "shifted_laplacian_1d(200, 0.9)",
      methods = c("minres", "gmres", "bicgstab", "lsqr", "lsmr"),
      truth = "closed-form eigendecomposition",
      build = function(cn) {
        n <- 200L
        sigma <- 0.9
        lam <- laplacian_1d_eigenvalues(n) - sigma
        apply_fn <- function(X) laplacian_1d_apply(X) - sigma * X
        set.seed(102)
        b <- stats::rnorm(n)
        list(A = counted_linop(apply_fn, apply_fn, c(n, n), cn,
                               properties = c(hermitian = TRUE)),
             b = b, truth = shifted_laplacian_solve(n, sigma, b),
             kappa = max(abs(lam)) / min(abs(lam)))
      })

  ## Nonsymmetric and matrix free, supplying the forward action and nothing else.
  ## n is 80 and mu is 0.1 because the closed form stops being ground truth long
  ## before the operator stops being well conditioned: the helper's table puts
  ## rho^n = 3e3 here, inside the documented budget of about 1e4.
  add(label = "convdiff_1d(80, 0.1)",
      methods = c("gmres", "bicgstab"),
      truth = "closed-form similarity",
      build = function(cn) {
        n <- 80L
        mu <- 0.1
        M <- convdiff_1d_apply(diag(1, n), mu)
        set.seed(103)
        b <- stats::rnorm(n)
        list(A = counted_linop(function(X) convdiff_1d_apply(X, mu), NULL,
                               c(n, n), cn),
             b = b, truth = convdiff_1d_solve(n, mu, b),
             kappa = {
               s <- svd(M)$d
               max(s) / min(s)
             })
      })

  ## Stored, hermitian positive definite, well conditioned, with a tridiagonal
  ## inverse in closed form.
  add(label = "kms(300, 0.7)",
      methods = c("cg", "minres", "gmres", "bicgstab", "lsqr", "lsmr"),
      truth = "closed-form inverse",
      build = function(cn) {
        n <- 300L
        rho <- 0.7
        M <- kms_matrix(n, rho)
        set.seed(104)
        b <- stats::rnorm(n)
        list(A = counted_dense(M, cn,
                               properties = c(hermitian = TRUE,
                                              positive_definite = TRUE)),
             b = b, truth = kms_inverse(n, rho) %*% b,
             kappa = {
               s <- svd(M)$d
               max(s) / min(s)
             })
      })

  ## Stored, with the condition number dialled to 1e4 exactly.
  add(label = "spd_prescribed(300, kappa 1e4)",
      methods = c("cg", "minres", "gmres", "bicgstab", "lsqr", "lsmr"),
      truth = "dense LU",
      build = function(cn) {
        n <- 300L
        M <- spd_prescribed(n, 10^seq(-4, 0, length.out = n), seed = 5L)
        set.seed(105)
        b <- stats::rnorm(n)
        list(A = counted_dense(M, cn,
                               properties = c(hermitian = TRUE,
                                              positive_definite = TRUE)),
             b = b, truth = solve(M, b), kappa = 1e4)
      })

  ## Complex hermitian positive definite. Nothing about the roster changes.
  add(label = "hpd_prescribed(200, kappa 1e3)",
      methods = c("cg", "minres", "gmres", "bicgstab", "lsqr", "lsmr"),
      truth = "dense LU",
      build = function(cn) {
        n <- 200L
        M <- hpd_prescribed(n, 10^seq(-3, 0, length.out = n), seed = 6L)
        set.seed(106)
        b <- complex(real = stats::rnorm(n), imaginary = stats::rnorm(n))
        list(A = counted_dense(M, cn,
                               properties = c(hermitian = TRUE,
                                              positive_definite = TRUE)),
             b = b, truth = solve(M, b), kappa = 1e3)
      })

  ## Rectangular, with a right-hand side outside the range of A, so the residual
  ## does not go to zero and the certificate reports Stewart's reading instead.
  add(label = "lsq_prescribed(600 x 200, kappa 1e3)",
      methods = c("lsqr", "lsmr"),
      truth = "closed-form pseudoinverse",
      build = function(cn) {
        m <- 600L
        n <- 200L
        f <- lsq_prescribed(m, n, 10^seq(0, -3, length.out = n), seed = 7L)
        set.seed(107)
        b <- stats::rnorm(m)
        outside <- vnorm(lsq_prescribed_residual(f, b)) / vnorm(b)
        list(A = counted_dense(f$A, cn), b = b,
             truth = lsq_prescribed_solve(f, b), kappa = 1e3,
             note = sprintf("||(I - U U^H) b|| / ||b|| = %.3f", outside))
      })

  ## Rectangular, matrix free and rank deficient: the constant vector is in the
  ## nullspace, so what is checked is which of the infinitely many solutions
  ## comes back. kappa is over the nonzero singular values.
  add(label = "diff_1d(400)",
      methods = c("lsqr", "lsmr"),
      truth = "closed-form minimum-norm solution",
      build = function(cn) {
        n <- 400L
        s <- diff_1d_singular_values(n)
        set.seed(108)
        b <- stats::rnorm(n - 1L)
        list(A = counted_linop(function(X) X[-1L, , drop = FALSE] - X[-n, , drop = FALSE],
                               function(X) rbind(0, X) - rbind(X, 0),
                               c(n - 1L, n), cn),
             b = b, truth = diff_1d_min_norm_solve(n, b),
             kappa = max(s) / min(s[s > 0]))
      })

  cases
}

bench_solvers <- function() {
  out <- bench_rows()

  for (case in solver_cases()) {
    cat(sprintf("  %s\n", case$label))
    cn <- new_counter()
    fx <- case$build(cn)
    A <- fx$A
    maxit <- 10L * as.integer(A$dim[2L])
    case_note <- fx$note %||% ""
    if (nzchar(case_note)) cat(sprintf("    %s\n", case_note))

    for (m in case$methods) {
      run <- tryCatch(
        bench_run(solve(A, fx$b, method = m, tol = SOLVER_TOL, maxit = maxit,
                        details = TRUE),
                  counter = cn),
        error = function(e) structure(list(message = conditionMessage(e)),
                                      class = "bench_error"))
      if (inherits(run, "bench_error")) {
        out$add(fixture = case$label, rows = A$dim[1L], cols = A$dim[2L],
                kappa = fx$kappa, method = m, iterations = NA_integer_,
                applies = NA_integer_, column_applies = NA_integer_,
                seconds = NA_real_, residual = NA_real_,
                backward_error = NA_real_, reading = NA_character_,
                residual_status = NA_character_, backward_status = NA_character_,
                overall = "refused", forward_error = NA_real_,
                truth = case$truth, note = run$message)
        next
      }
      fit <- run$value
      cert <- fit$certificate
      out$add(fixture = case$label, rows = A$dim[1L], cols = A$dim[2L],
              kappa = fx$kappa, method = m, iterations = fit$iterations,
              applies = run$applies, column_applies = run$cols,
              seconds = run$seconds,
              residual = max(cert$values$residual),
              backward_error = max(cert$values$backward_error),
              reading = cert_reading(cert),
              residual_status = cert_line(cert, "residual"),
              backward_status = cert_line(cert, "backward error"),
              overall = cert$overall,
              forward_error = rel_error(fit$x, fx$truth),
              truth = case$truth, note = case_note)
      cat(sprintf("    %-9s %5d it  %6d applies  %7.3f s  backward %.2e  %s\n",
                  m, fit$iterations, run$applies, run$seconds,
                  max(cert$values$backward_error), cert$overall))
    }

    ## What the certificate's floor rests on, priced on its own.
    nrun <- bench_run(norm2(A), counter = cn)
    out$add(fixture = case$label, rows = A$dim[1L], cols = A$dim[2L],
            kappa = fx$kappa, method = "norm2()", iterations = NA_integer_,
            applies = nrun$applies, column_applies = nrun$cols,
            seconds = nrun$seconds, residual = NA_real_,
            backward_error = NA_real_, reading = NA_character_,
            residual_status = NA_character_, backward_status = NA_character_,
            overall = NA_character_, forward_error = NA_real_,
            truth = NA_character_,
            note = sprintf("||A|| ~ %.4g by %s, included in every solve row above",
                           nrun$value$value, nrun$value$method))
  }

  out$collect()
}
