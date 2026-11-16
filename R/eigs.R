## Section 7.2, the eigensolver verb, and the reference hermitian engine.
## Lanczos with full reorthogonalisation, locking, and explicit restarts, for an
## operator established hermitian; Arnoldi with Krylov-Schur restarting
## (R/eigs-arnoldi.R) for every other square operator.
##
## What each requires. Lanczos needs a hermitian operator and a forward action;
## Arnoldi needs the forward action alone. Every apply in both is mode "N", so
## neither asks for an adjoint, and for Lanczos the reason is structural -- a
## hermitian operator IS its own adjoint, so there is no second callback to ask
## for.
##
## The hermitian requirement is a capability in the backend registry, so
## method = "auto" reads it the way solve()'s auto reads its minima: an operator
## that establishes it gets Lanczos (or PRIMME), and one that does not gets
## Arnoldi, with the dispatch table naming the reason. A named method is the
## caller's choice and is refused with the same reason. There is still no
## evidence minimum on the declaration: a bare user_declaration of hermitian
## selects Lanczos, and is carried into the certificate instead, where the
## forward-error row records it in depends_on, so a bound resting on it fails any
## requirement the declaration would have failed directly.
##
## The Rayleigh quotient is what gets reported, not the Ritz value. For a plain
## Lanczos run the two agree to rounding, since a Ritz value of a hermitian
## operator is the Rayleigh quotient of its own Ritz vector. Under shift-invert
## they do not: sigma + 1/theta inherits whatever the inner solve got wrong, and
## x^H A x / (x^H x) is the minimiser of ||A x - mu x|| over mu, so it is the value
## for which the certificate's exhibited perturbation is smallest. Measuring the
## eigenvalue on A rather than reading it off the transformed problem is the same
## discipline as the outer loop measuring b - A x, and it is what makes an inexact
## shift-invert reportable rather than merely fast.
##
## Arithmetic. alpha is a Rayleigh quotient and beta is a norm, so the projected
## tridiagonal is real and symmetric whether the operator is real or complex, and
## its eigenvectors are real. Only the basis carries complex storage.

## The wanted eigenvalues under shift-invert are the ones nearest sigma, which is
## the largest end of the transformed spectrum whatever the caller asked about the
## original one. Selecting by anything else would need a second transformation and
## would not mean what it said.
SIGMA_SELECT <- "LM"

## The GMRES cycle length of the inner solve under a shift, for an operator not
## established hermitian. See shift_invert_operator().
SHIFT_INVERT_RESTART <- 100L

#' Eigenvalues and eigenvectors of an operator
#'
#' A few eigenpairs of a square operator, computed from its action alone. An
#' operator established hermitian gets Lanczos with full reorthogonalisation;
#' any other gets Arnoldi with Krylov-Schur restarting. The result records what
#' converged, what residual each pair reached, and what argument establishes the
#' bound on the eigenvalues.
#'
#' The eigenvalue reported for each pair is the Rayleigh quotient
#' `x^H A x / (x^H x)`, measured on `A` itself. It is the value that minimises
#' `||A x - mu x||`, so no other choice of `mu` certifies better.
#'
#' For an operator that is not established hermitian the values are complex
#' wherever the spectrum is, and a real operator's non-real eigenvalues come in
#' conjugate pairs that are computed together. The eigenvectors of a non-normal
#' operator are not orthogonal, and a small residual places a value in the
#' operator's pseudospectrum rather than near an eigenvalue: the certificate's
#' `forward error` row says so instead of claiming a bound. An operator with no
#' declaration at all runs too, through Arnoldi; declaring `hermitian = TRUE`
#' buys the real recurrence and Weyl's bound.
#'
#' `target identity` is `not_checked`, always. A small residual proves the pair is
#' an approximate eigenpair; that it is the largest, smallest or nearest one
#' requires an inertia count, an enclosure or a separation bound, and matrix-free
#' there is generally none. See `vignette("certificates")`.
#'
#' @param A A square `linop`. No adjoint is required.
#' @param k How many eigenpairs.
#' @param which Which end of the spectrum: `"largest"` or `"smallest"` by
#'   magnitude, `"largest_algebraic"` or `"smallest_algebraic"` by value (by real
#'   part for a complex spectrum, also spelled `"largest_real"` and
#'   `"smallest_real"`), and `"largest_imaginary"` or `"smallest_imaginary"` by
#'   the magnitude of the imaginary part. RSpectra's `"LM"`, `"SM"`, `"LA"`,
#'   `"SA"`, `"LR"`, `"SR"`, `"LI"` and `"SI"` are accepted for the same orders.
#'   A hermitian spectrum is real, so the two imaginary orders are refused for an
#'   operator established hermitian.
#' @param sigma A shift. Given one, the eigenvalues nearest `sigma` are sought,
#'   by running the same recurrence on `(A - sigma I)^-1` through the package's
#'   own solvers: MINRES for a hermitian operator and a real shift, GMRES
#'   otherwise. A complex shift is accepted for any operator and runs Arnoldi.
#'   `which` does not apply and passing both is an error.
#' @param B For the generalized problem `A x = lambda B x`: a square `linop` or
#'   matrix of the same shape as `A`, or `NULL` for the standard problem. With
#'   `A` established hermitian and `B` established hermitian positive definite,
#'   Lanczos runs on `B^-1 A` in the `B` inner product, the returned vectors are
#'   `B`-orthonormal and the `forward error` row bounds each value; otherwise
#'   Arnoldi runs on `B^-1 A`. A shift runs on `(A - sigma B)^-1 B`. Each step
#'   solves with `B`, or with `A - sigma B`, through the package's own solvers,
#'   and `inner` carries the arguments. Values are measured as
#'   `x^H A x / x^H B x` and the backward error perturbs both operators,
#'   `||A x - mu B x|| / ((||A|| + |mu| ||B||) ||x||)`.
#' @param tol Convergence tolerance. Every method holds each pair, measured on
#'   `A`, to the backward error `||A x - mu x|| / (||A|| ||x||)` and to the
#'   relative residual `||A x - mu x|| / (|mu| ||x||)` wherever `|mu|` is above
#'   the arithmetic floor, each plus that floor: the certificate's
#'   `backward error` and `residual` rows.
#' @param maxit Budget of applies across all restarts, of the operator the
#'   iteration runs on: `A`, or under `sigma` the inner solve. `NULL` takes
#'   `20 n` for an `n x n` operator, twenty times what forming it column by
#'   column would cost. Hard ends of a spectrum need applies that grow with `n`
#'   at a fixed `ncv` (seven times `n` for the largest eigenvalue of the 1-D
#'   Laplacian at `n = 2000`), so the default grows with it. A run that stops
#'   improving ends before the budget does, and a run that exhausts it returns a
#'   `fail` certificate.
#' @param ncv Subspace size per restart. Storage is `ncv` vectors. `NULL`
#'   takes the engine's default: `max(2k + 1, k + 20)` for `"lanczos"` and
#'   `"arnoldi"`, and PRIMME's own basis size for `"primme"`.
#' @param v0 Starting vector, or `NULL` for a seeded random one.
#' @param seed Seed for the starting vector. The caller's random stream is
#'   restored afterwards.
#' @param method `"auto"`, `"lanczos"`, `"arnoldi"`, `"primme"` or `"rspectra"`.
#'   `"lanczos"` is the reference implementation for a hermitian operator and
#'   `"arnoldi"` for any square one. `"primme"` is the PRIMME library, available
#'   for a hermitian operator when the `linop.primme` package is installed, which
#'   restarts implicitly and so keeps a basis of fixed small size. `"rspectra"` is
#'   the Spectra library through the RSpectra package, for a real operator only,
#'   and stops on its own convergence rule; `"auto"` never takes it. `"auto"` takes PRIMME
#'   when it is available and the operator is established hermitian, Lanczos when
#'   only the second holds, and Arnoldi otherwise. Established means computed,
#'   constructed or proven: a bare `hermitian = TRUE` declaration on a callback
#'   is acted on only by naming `"lanczos"` or `"primme"`, as `solve()` acts on
#'   one only for a named method. The certificate's `dispatch`
#'   records every engine considered and why. Whichever runs, the reported
#'   values are Rayleigh quotients measured on `A` and the certificate is the
#'   same. For a Jacobi operator with periodic limits, `"second_order"`
#'   solves the quadratic relative-spectrum pencil and `"dissipative"` solves
#'   the section with an imaginary projection barrier. Both require a real
#'   `sigma` strictly inside a gap between essential bands. They use dense
#'   base R eigenproblems and the section width budget; `maxit`, `ncv`, `v0`,
#'   `B`, `inner` and `norm_control` are refused. Their extraction is
#'   deterministic and independent of `seed`. Values are real parts of the
#'   computed roots. Each returned `enclosures` row is validated on the original
#'   operator. An interval wholly inside the gap contains a discrete eigenvalue.
#'   The certificate reports its precision separately from its isolation.
#' @param barrier For `method = "dissipative"`, a list with positive `gamma`
#'   (default `0.01`) and projection half-width `n`. The projection is onto
#'   sites `abs(j) <= n`, must cover the perturbation and must end strictly
#'   inside the section. With `n = NULL`, it covers half the free tail and
#'   grows with the section. `Im(root)/gamma > 1/2` is the barrier separation
#'   diagnostic, recorded as an estimate. Spectral inclusion is independently
#'   validated on the original operator. A finite calculation establishes
#'   existence in its isolated intervals; completeness of the gap stays
#'   unchecked. For these two methods `floor_const` adds to a mandatory
#'   operation-count roundoff bound.
#' @param inner Arguments for the inner solve when `sigma` or `B` is given,
#'   such as `method` or `tol`.
#' @param section Arguments for the truncation, when `A` is an operator on a
#'   sequence space and has no matrix: `n_start` for how much free tail the first
#'   section leaves beyond the window, and `n_max` for the widest section the
#'   search may reach. The width itself is not among them, because
#'   `eigs(finite_section(A, n))` is already how a caller fixes it.
#' @param floor_const `c` in the arithmetic floor of the certificate.
#' @param norm_control Arguments for the `||A||` estimate.
#' @return An object of class `linop_eigen`, with `values`, `vectors` and the
#'   `certificate`. A run that exhausts its budget returns the best pairs it
#'   reached and a `fail` certificate rather than an error, and it may hold fewer
#'   than `k` of them: nothing is padded. On a sequence space it also carries the
#'   width `n` that was chosen and the `widths` that were tried. The gap methods
#'   also return `roots`, `enclosures` (columns `lower`, `upper`) and
#'   `all_roots` from the extraction. For `"second_order"`, `all_roots` holds
#'   the upper-half-plane representatives. The exact second-order enclosure
#'   radius is `abs(Im(root))`; the returned radius is recomputed from its
#'   vector to validate floating-point roots, including near-real pairs.
#'   `certificate(fit)$values$relative_enclosures` additionally holds the
#'   second-order theorem's intervals enlarged by the checked pencil defect.
#' @examples
#' n <- 40
#' L <- linop(function(X) rbind(X[-1, , drop = FALSE], 0) +
#'                        rbind(0, X[-n, , drop = FALSE]) - 2 * X,
#'            dim = c(n, n), properties = c(hermitian = TRUE))
#' fit <- eigs(L, k = 3, which = "largest_algebraic")
#' fit$values
#'
#' ## a non-normal operator, convection-diffusion: its spectrum is real and
#' ## known, and its eigenvectors are far from orthogonal
#' C <- linop(function(X) 1.5 * rbind(X[-1, , drop = FALSE], 0) +
#'                        0.5 * rbind(0, X[-n, , drop = FALSE]) - 2 * X,
#'            dim = c(n, n))
#' fit <- eigs(C, k = 2, which = "largest_real")
#' fit$values - (-2 + 2 * sqrt(0.75) * cos(pi * (1:2) / (n + 1)))
#'
#' ## on an operator with no matrix, where the width is chosen rather than given
#' H <- linop_jacobi(diagonal = 1)
#' fit <- eigs(H, k = 1, which = "largest_algebraic")
#' fit$n
#' fit$values - sqrt(5)
#' @export
eigs <- function(A, k, which = "largest", sigma = NULL, B = NULL,
                 tol = 1e-10, maxit = NULL, ncv = NULL, v0 = NULL, seed = 1L,
                 method = "auto", inner = list(), section = list(),
                 floor_const = SOLVE_FLOOR_CONST, norm_control = list(), barrier = list()) {
  if (!is_linop(A)) stopf("eigs() expects a linop")
  if (!is.list(section)) stopf("section must be a list")
  own <- if (is_scalar_string(method)) get_node(A$node)$methods[[method]]
  if (!is.null(own)) {
    return(own(A, k, section, list(
      which = which, sigma = sigma, B = B, tol = tol, maxit = maxit, ncv = ncv,
      v0 = v0, seed = seed, method = method, inner = inner, barrier = barrier,
      floor_const = floor_const, norm_control = norm_control)))
  }
  owners <- if (is_scalar_string(method)) node_method_owners(method) else character()
  if (length(owners)) {
    stopf("method = '%s' is defined on the %s node type%s, and this operator is a '%s' node",
          method, paste(sprintf("'%s'", owners), collapse = " and "),
          if (length(owners) > 1L) "s" else "", A$node)
  }
  if (length(barrier)) stopf("barrier applies to method = 'dissipative'")
  if (!is.null(B)) {
    if (is.matrix(B) || inherits(B, "Matrix")) B <- linop(B)
    if (!is_linop(B)) stopf("B must be a linop or a matrix")
    if (!is_finite_dim(A)) {
      stopf(paste0("eigs() takes B for an operator of finite dimension; this one is %s.\n",
                   "  A generalized problem on a sequence space has no truncation rule here."),
            fmt_dim(A$dim))
    }
  }
  ## An operator with no matrix is not refused where its node type knows how to
  ## make one. Everything the caller supplied is forwarded unchanged, because the
  ## width is the only thing they did not supply and so the only thing chosen
  ## here.
  if (!is_finite_dim(A)) {
    run <- get_node(A$node)$spectrum
    ## No route to a matrix: the ordinary refusal, naming this verb.
    if (is.null(run)) require_finite_dim(A, "eigs")
    return(run(A, k, section, list(
      which = which, sigma = sigma, B = B, tol = tol, maxit = maxit, ncv = ncv,
      v0 = v0, seed = seed, method = method, inner = inner,
      floor_const = floor_const, norm_control = norm_control)))
  }
  if (length(section)) {
    stopf(paste0("section applies to an operator on a sequence space, and this one is %s.\n",
                 "  Its truncation is already fixed, so there is no width to choose."),
          fmt_dim(A$dim))
  }
  require_finite_dim(A, "eigs")
  n <- A$dim[1L]
  if (A$dim[1L] != A$dim[2L]) {
    stopf(paste0("eigs() needs a square operator; this one is %s.\n",
                 "  The singular triplets of a rectangular operator are svds()."),
          fmt_dim(A$dim))
  }
  if (!is.null(B)) {
    require_finite_dim(B, "eigs")
    if (!identical(as.numeric(B$dim), as.numeric(A$dim))) {
      stopf("B is %s and A is %s; a pencil needs the two the same shape",
            fmt_dim(B$dim), fmt_dim(A$dim))
    }
    if (!is.null(get_node(A$node)$certify)) {
      stopf("eigs() takes B for an operator certified as an eigenproblem of its own; a %s node certifies its own rows",
            A$node)
    }
  }
  herm <- isTRUE(capv(A, "hermitian"))
  definite <- !is.null(B) && isTRUE(capv(B, "hermitian")) &&
              isTRUE(capv(B, "positive_definite"))
  real_spectrum <- herm && (is.null(B) || definite)
  which_code <- normalise_which(which, EIGEN_WHICH, "eigs")
  if (real_spectrum && which_code %in% c("LI", "SI")) {
    stopf(paste0("which = '%s' orders by imaginary part, and this %s is established\n",
                 "  hermitian%s, so every eigenvalue is real and every imaginary part is zero."),
          which, if (is.null(B)) "operator" else "pencil",
          if (is.null(B)) "" else " with B positive definite")
  }
  s <- eigen_setup(A, k, maxit, "eigs", n)
  k <- s$k; maxit <- s$maxit

  ## Selection happens on the spectrum of the operator the recurrence runs on,
  ## which under a shift is the solve; measurement always happens on A.
  if (is.null(sigma)) {
    sel <- which_code
  } else {
    if (!identical(which, "largest")) {
      stopf(paste0("sigma and which cannot both be given.\n",
                   "  A shift asks for the eigenvalues nearest sigma, which fixes the end of the\n",
                   "  transformed spectrum the recurrence works at. Drop which, or drop sigma and\n",
                   "  ask for an end of the spectrum of A directly."))
    }
    if (!(is.numeric(sigma) || is.complex(sigma)) || length(sigma) != 1L ||
        !is.finite(sigma)) {
      stopf("sigma must be a single finite number")
    }
    ## A complex shift with no imaginary part is a real one.
    if (is.complex(sigma) && Im(sigma) == 0) sigma <- Re(sigma)
    sel <- SIGMA_SELECT
  }

  ## Which engine runs the recurrence. Every engine returns the same thing, the
  ## pairs with their Rayleigh quotients and residuals measured on A, so nothing
  ## below this line knows which one ran.
  pick <- select_backend("eigs", A,
                         list(k = k, which = sel, sigma = sigma, tol = tol, ncv = ncv, maxit = maxit, B = B,
                              interior = sel %in% c("LM", "SM")),
                         method, "eigs")
  ## The shift-invert operator is the hermitian one exactly when the engine
  ## requires a hermitian operator, since that engine runs on it.
  herm_engine <- "hermitian" %in% pick$backend$capabilities
  Op <- if (!is.null(sigma)) shift_invert_operator(A, sigma, inner, herm_engine, B)
        else if (!is.null(B)) generalized_operator(A, B, inner)
        else A
  ncv <- eigen_ncv(ncv %||% pick$backend$default_ncv(A, k), k, n)
  run <- pick$backend$run(Op, A, list(k = k, sel = sel, tol = tol, maxit = maxit,
                                      ncv = ncv, v0 = v0, seed = seed,
                                      floor_const = floor_const, B = B), pick$engine)
  label <- if (is.null(sigma)) {
             if (is.null(B)) pick$backend$label
             else sprintf("%s on B^-1 A", pick$backend$label)
           } else sprintf("%s on (A - %s %s)^-1%s", pick$backend$label, format(sigma),
                          if (is.null(B)) "I" else "B", if (is.null(B)) "" else " B")
  chosen <- if (is.null(sigma)) pick$backend$name
            else paste(pick$backend$name, "shift-invert")

  ## Which certificate the result carries is the node type's to say. A finite
  ## section answers a question about the operator underneath it, so its rows are
  ## a different set; everything else gets the eigenpair certificate.
  ctx <- list(A = A, values = run$values, vectors = run$vectors, tol = tol,
              iterations = run$steps, maxit = maxit, floor_const = floor_const,
              norm_control = norm_control, norm_lower_bound = run$anorm_lb,
              requested = k, stop_reason = run$stop_reason, B = B,
              norm_lower_bound_B = run$bnorm_lb %||% 0)
  cert <- (get_node(A$node)$certify %||% eigen_certify)(ctx)

  new_eigen_result(values = run$values, vectors = run$vectors,
                   certificate = cert, k = k, which = which, method = chosen,
                   iterations = run$steps, maxit = maxit, rounds = run$rounds,
                   converged = cert$values$converged, residuals = run$residuals,
                   dispatch = backend_dispatch(method, pick, paste0(
                     if (is.null(B)) {
                       if (herm) "hermitian operator, "
                       else "operator not established hermitian, "
                     } else if (real_spectrum) "hermitian pencil with B positive definite, "
                     else "pencil not established hermitian with B positive definite, ",
                     label)),
                   engine = run$engine)
}

## The default: the eigenpair certificate, on the operator the run measured
## against. The ||A|| estimate is taken here rather than in eigs() because a node
## certifying its own results may know the norm in closed form and have no use for
## a power iteration.
eigen_certify <- function(ctx) {
  A <- ctx$A; B <- ctx$B
  norm_est <- do.call(norm2, c(list(A = A), ctx$norm_control))
  ## The forward bound and the orthogonality row rest on A being hermitian, and
  ## for a pencil on B being hermitian positive definite as well; that evidence
  ## goes under the row with the theorem.
  definite <- is.null(B) || (isTRUE(capv(B, "hermitian")) &&
                             isTRUE(capv(B, "positive_definite")))
  herm_ev <- if (isTRUE(capv(A, "hermitian")) && definite) cape(A, "hermitian")
  eigen_certificate(A, ctx$values, ctx$vectors, tol = ctx$tol,
                    norm_estimate = norm_est,
                    iterations = ctx$iterations, maxit = ctx$maxit,
                    floor_const = ctx$floor_const,
                    norm_lower_bound = ctx$norm_lower_bound,
                    hermitian_evidence = herm_ev,
                    requested = ctx$requested,
                    stop_reason = ctx$stop_reason, subject = "eigen",
                    B = B,
                    norm_estimate_B = if (!is.null(B))
                      do.call(norm2, c(list(A = B), ctx$norm_control)),
                    norm_lower_bound_B = ctx$norm_lower_bound_B %||% 0,
                    definite_evidence = if (!is.null(B) && !is.null(herm_ev))
                      list(cape(B, "hermitian"), cape(B, "positive_definite")))
}

## (A - sigma I)^-1 as an operator, built out of the package's own solvers.
##
## Section 4.2 has this the right way round: shift-invert holds a solve object,
## not an operator, and what makes it usable as one here is that the inner method
## is a fixed iteration to a fixed tolerance. Whatever it gets wrong shows up in
## the residual the outer loop measures on A, which is why nothing downstream has
## to trust it.
##
## For a hermitian engine, hermitian by construction given that A is: the shifted
## operator inherits the declaration through caps_sum() and caps_scale(), and its
## inverse is hermitian with it. That is also why the adjoint callback is the
## apply -- a hermitian operator is its own adjoint, so norm2() gets its power
## iteration for nothing. For Arnoldi nothing is declared and nothing is needed:
## the inner method is GMRES, and no adjoint is supplied because Arnoldi applies
## only mode "N". A complex shift makes the operator complex whatever A is.
##
## For a pencil the operator is (A - sigma B)^-1 B, whose eigenvalues are
## 1 / (lambda - sigma) on the pencil's eigenvectors (ARPACK's mode 3). It is
## self-adjoint in the B inner product rather than the ordinary one, so nothing
## is declared about it; the inner solve is MINRES when A - sigma B is
## hermitian, which a hermitian engine has required of both.
shift_invert_operator <- function(A, sigma, inner, hermitian = TRUE, B = NULL) {
  n <- A$dim[1L]
  if (hermitian && is.complex(sigma)) {
    stopf(paste0("sigma must be real for a hermitian engine.\n",
                 "  A complex shift leaves A - sigma I non-hermitian; method = 'arnoldi' takes it."))
  }
  dtype <- if (is.complex(sigma) || (!is.null(B) && B$dtype == "complex")) "complex"
           else A$dtype
  shifted <- A - sigma * (B %||% linop_eye(n, dtype))
  ## Tighter than the outer tolerance by default: what the inner solve leaves
  ## behind is amplified by the transformation, and the outer loop can only
  ## report it, not remove it.
  ## GMRES restarted at its default of 30 stalls on the indefinite operator an
  ## interior shift makes: on convection-diffusion at n = 100, sigma = -1.3, it
  ## spends 1000 applies and stops at a relative residual of 6.4e-2. A longer
  ## cycle is the default here, bounded so that its basis stays a modest multiple
  ## of the outer one; inner = list(control = list(restart = )) overrides it.
  defaults <- if (hermitian) list(method = "minres", tol = 1e-12)
              else list(method = "gmres", tol = 1e-12,
                        control = list(restart = min(n, SHIFT_INVERT_RESTART)))
  args <- utils::modifyList(defaults, inner)
  apply_inverse <- function(X) {
    if (!is.null(B)) X <- linop_apply(B, X, "N")
    z <- do.call(solve, c(list(shifted, X), args))
    attr(z, "certificate") <- NULL
    z
  }
  if (!hermitian || !is.null(B)) {
    return(new_linop("fun", c(n, n), dtype, new_caps(),
                     list(apply = apply_inverse, adjoint = NULL)))
  }
  new_linop("fun", c(n, n), dtype,
            do.call(new_caps, list(
              hermitian = capability(TRUE, ev_construction(
                depends_on = Filter(Negate(is.null),
                                    list(cape(shifted, "hermitian"))))))),
            list(apply = apply_inverse, adjoint = apply_inverse))
}

## B^-1 A as an operator, for a pencil with no shift (ARPACK's mode 2): its
## eigenpairs are the pencil's, and each apply is one apply of A and one solve
## with B by the package's own solvers. inner carries that solve's arguments.
##
## The solve's error is absolute at the scale of ||B^-1 A v||, set by the large
## end of the spectrum, and lands on every eigencomponent, so the small end sees
## it relative to its own eigenvalue. On the linear-element mass matrix pencil of
## laplacian_1d(100), three smallest (lambda near 1e-3, the largest near 12), an
## inner tolerance of 1e-12 plateaus the residual row at 6.1e-10 and converges
## none in 1691 applies; 1e-14 converges all three in 310. Hence the tighter
## default here than under a shift, where the transformation favours the wanted
## end instead.
generalized_operator <- function(A, B, inner) {
  n <- A$dim[1L]
  dtype <- if (A$dtype == "complex" || B$dtype == "complex") "complex" else "double"
  args <- utils::modifyList(list(tol = 1e-14), inner)
  apply_solve <- function(X) {
    z <- do.call(solve, c(list(B, linop_apply(A, X, "N")), args))
    attr(z, "certificate") <- NULL
    z
  }
  new_linop("fun", c(n, n), dtype, new_caps(),
            list(apply = apply_solve, adjoint = NULL))
}

## The outer loop. A round builds a subspace, extracts, measures on A, locks what
## converged and restarts thickly from what did not.
##
## Locking is what full reorthogonalisation buys beyond stability: a converged
## eigenvector joins the set every later basis is orthogonalised against, so the
## next round searches the orthogonal complement and cannot return the same pair
## twice. A degenerate eigenvalue is the case that shows this is deflation rather
## than exclusion -- a second, orthogonal eigenvector for the same eigenvalue is
## still reachable, and is still correct to return.
##
## Restarting thickly rather than from one vector is not a refinement, it is what
## makes the method converge at all. A restart from the sum of the unconverged
## Ritz vectors throws away the subspace that produced them: on laplacian_1d(60)
## at ncv = 24, asking for four pairs, that spends 240 of a 300-iteration budget
## over 9 rounds, converges none of them and stalls at a backward error of 5.7e-4,
## where this one reaches 5.1e-13 in 96. Each round was rediscovering from one
## direction what the last one spent 24 applies learning.
##
## Wu and Simon's thick restart keeps the wanted Ritz vectors themselves as the
## start of the next basis. They span an approximately invariant subspace already,
## so in exact arithmetic the projected problem keeps them diagonal at the values
## they were extracted with and one arrow row ties them to the new directions;
## lanczos_run() measures both rather than assuming them.
## Where a request fits in one round there is no restart to differ over and the
## two agree to the last digit, which is the control that confines the difference
## to this block; dev_notes/spikes/restart-comparison.R runs both.
lanczos_core <- function(Op, A, k, sel, tol, maxit, ncv, v0, seed, floor_const,
                         B = NULL, prelocked = NULL) {
  n <- Op$dim[1L]
  eps <- .Machine$double.eps
  metric <- !is.null(B)

  Qv <- NULL
  BQv <- NULL
  lock_theta <- numeric(0)
  lock_value <- numeric(0)
  lock_res <- numeric(0)

  v <- eigen_start_vector(n, Op$dtype, v0, seed)
  keep <- NULL
  steps <- 0L
  rounds <- 0L
  anorm_lb <- 0
  bnorm_lb <- if (metric) 0 else NULL
  op_is_a <- identical(Op, A)
  stop_reason <- NULL
  progress <- NULL

  ## What comes back when the budget runs out before k pairs converge: the best
  ## pairs the last round produced, with their measured residuals, so a failed
  ## run reports a result and a certificate rather than an error.
  best_theta <- numeric(0); best_value <- numeric(0)
  best_vec <- NULL; best_res <- numeric(0)

  repeat {
    nl <- if (is.null(Qv)) 0L else ncol(Qv)
    want <- k - nl
    if (want <= 0L) break
    if (steps >= maxit) { stop_reason <- "budget exhausted"; break }
    npre <- if (is.null(prelocked)) 0L else ncol(prelocked)
    if (nl + npre >= n) { stop_reason <- "the locked set spans the operator"; break }
    rounds <- rounds + 1L

    ## The subspace this round may occupy, and how much of the budget is left to
    ## fill it with. Retained vectors take room in the first and cost nothing
    ## against the second, which is the whole economy of a thick restart.
    m <- min(ncv - nl, n - nl - npre)
    if (!is.null(keep)) {
      p <- min(ncol(keep$Y), m - 2L)
      keep <- if (p >= 1L) list(Y = keep$Y[, seq_len(p), drop = FALSE],
                                BY = if (metric) keep$BY[, seq_len(p), drop = FALSE],
                                H = keep$H[seq_len(p), seq_len(p), drop = FALSE])
              else NULL
    }

    ## The start of a round lives in the complement of what is locked. A start
    ## that lies entirely inside the locked space carries no new direction, and a
    ## fresh random one does.
    Qd <- cbind(prelocked, Qv)
    v <- orth_against(Qd, v, BQv)
    if (!is.null(keep)) v <- orth_against(keep$Y, v, keep$BY)
    nv <- col_norms(v)
    if (nv <= eps * sqrt(n)) {
      v <- eigen_random_vector(n, Op$dtype, seed, rounds)
      v <- orth_against(Qd, v, BQv)
      if (!is.null(keep)) v <- orth_against(keep$Y, v, keep$BY)
      nv <- col_norms(v)
      if (nv <= eps * sqrt(n)) {
        stop_reason <- "no direction left outside the locked subspace"
        break
      }
    }
    ## Unit in the inner product the basis is orthonormal in.
    if (metric) nv <- sqrt(col_dot(v, linop_apply(B, v, "N")))
    v <- v / nv

    run <- lanczos_run(Op, v, Qd, m, maxit - steps, keep, B, BQv)
    steps <- steps + run$applies
    if (run$applies == 0L) { stop_reason <- "the recurrence made no step"; break }
    ## Each step's ||A v_j|| for a unit v_j is a lower bound on ||A|| at no cost,
    ## and at an end of the spectrum near zero the measured pairs alone give one
    ## orders of magnitude short. Under a shift the step applies the solve.
    if (op_is_a) anorm_lb <- max(anorm_lb, run$apply_norm)

    ## The projected problem, hermitian, and real for a real operator.
    ee <- eigen(run$Tm, symmetric = TRUE)
    ord <- ritz_order(ee$values, sel)
    take <- ord[seq_len(min(want, run$size))]
    theta <- ee$values[take]
    Xr <- run$V %*% ee$vectors[, take, drop = FALSE]
    BXr <- if (metric) run$BV %*% ee$vectors[, take, drop = FALSE]

    ## Measured on A, in every case. Under a shift the recurrence has been
    ## working on a different operator entirely, and this is where that stops
    ## mattering.
    meas <- eigen_measure(A, Xr, B = B)
    mu <- meas$mu
    rr <- meas$rr
    anorm_lb <- max(anorm_lb, meas$anorm_lb)
    if (metric) bnorm_lb <- max(bnorm_lb, meas$bnorm_lb)

    ## A pair locks when it meets both rows of the certificate it will be read
    ## by, the backward error and the residual, each with the arithmetic floor.
    ## The backward error alone lets a pair with |mu| well below ||A|| lock and
    ## then fail the residual row.
    conv <- eigen_pair_bounds(rr, mu, anorm_lb, 1, tol, floor_const,
                              bnorm_lb, meas$bx)$converged
    if (any(conv)) {
      Qv <- cbind(Qv, Xr[, conv, drop = FALSE])
      if (metric) BQv <- cbind(BQv, BXr[, conv, drop = FALSE])
      lock_theta <- c(lock_theta, theta[conv])
      lock_value <- c(lock_value, mu[conv])
      lock_res <- c(lock_res, rr[conv])
    }
    ## The pairs this round did not lock, which are what a run ending here adds
    ## after the locked set; a locked pair is already in it.
    open <- !conv
    best_theta <- theta[open]; best_value <- mu[open]
    best_vec <- Xr[, open, drop = FALSE]; best_res <- rr[open]

    if (sum(conv) + nl >= k) break

    ## A round that locked nothing and moved no wanted pair, twice running, has
    ## nothing left to recover. Two rather than one, because a thick restart can
    ## spend a round rebuilding the directions it dropped.
    progress <- restart_progress(progress, rr, mu, anorm_lb, floor_const, conv)
    if (progress$stalled >= 2L) {
      stop_reason <- "the subspace stopped improving"
      break
    }

    ## The thick restart. Everything the round produced that was not locked is
    ## available to retain, in the order the request puts it; the wanted ones plus
    ## about half the subspace. That count is at the optimum of a fixed rule: on
    ## laplacian_1d(2000), one largest at ncv = 21, keeping 1, 3, 5, 8, 10, 13, 16
    ## or 19 costs 148561, 45957, 26037, 15634, 14255, 14517, 17606 or 34243
    ## applies, and a per-round choice maximising (m - p) sqrt(gap ratio) costs
    ## 32388. The retained block of the projected matrix is Z^H T Z, carried
    ## rather than set to the Ritz values, and the next run measures its coupling
    ## to the new directions instead of taking beta_last * s_i[last] for it.
    left_want <- want - sum(conv)
    avail <- setdiff(ord, take[conv])
    p <- min(length(avail), max(left_want + 3L, m %/% 2L), m - 2L)
    if (is.null(run$vnext) || p < 1L) {
      ## An invariant subspace, or no room to retain anything: there is no
      ## direction to continue from and the next round starts fresh.
      keep <- NULL
      v <- eigen_random_vector(n, Op$dtype, seed, rounds)
    } else {
      Z <- ee$vectors[, avail[seq_len(p)], drop = FALSE]
      H <- cross_adjoint(Z, run$Tm %*% Z)
      keep <- list(Y = run$V %*% Z, BY = if (metric) run$BV %*% Z,
                   H = (H + t(Conj(H))) / 2)
      v <- run$vnext
    }
  }

  ## What is returned: everything locked, then whatever the last round had left,
  ## up to k. An unconverged pair comes back with its residual and is reported as
  ## unconverged by the certificate rather than withheld.
  values <- lock_value
  vectors <- Qv
  residuals <- lock_res
  thetas <- lock_theta
  short <- k - length(values)
  if (short > 0L && length(best_value)) {
    keep <- seq_len(min(short, length(best_value)))
    values <- c(values, best_value[keep])
    vectors <- cbind(vectors, best_vec[, keep, drop = FALSE])
    residuals <- c(residuals, best_res[keep])
    thetas <- c(thetas, best_theta[keep])
  }
  if (!length(values)) {
    stopf("eigs() produced no pair at all; the operator applied to the starting vector is zero")
  }

  ord <- ritz_order(thetas, sel)
  list(values = values[ord], vectors = vectors[, ord, drop = FALSE],
       residuals = residuals[ord], thetas = thetas[ord], steps = steps, rounds = rounds,
       anorm_lb = anorm_lb, bnorm_lb = bnorm_lb, stop_reason = stop_reason)
}

## A Krylov space grown from one vector reaches one vector of each eigenspace in
## exact arithmetic, so the other copies of a repeated eigenvalue come only from
## rounding, and a request for k values can come back with a copy missing and
## the next distinct value in its place: the circle Laplacian's three smallest
## returned 0, 1, 4 where the truth is 0, 1, 1. A restarted block Krylov space
## would reach the copies directly, but at this package's default subspace a
## block of k leaves each restart a polynomial of degree about ncv / k: on
## laplacian_1d(400), four largest, a block of four converged one pair in 4703
## applies where one vector converges four in 764.
##
## So the search is finished and then checked. Deflate everything found, start
## from a fresh random direction, which has an O(1) component along a missed
## copy, and ask for the one best pair left. If it ranks above the worst pair
## returned, it was missed: it takes that place, the displaced pair stays in the
## deflation set, and the check repeats. If not, nothing better was reached from
## a random start and the search stops. One request for k = 1 per check, and
## none for k = 1, where any copy is a right answer.
deflation_sweep <- function(core, run, k, sel, tol, maxit, seed) {
  if (k <= 1L || !is.null(run$stop_reason) || length(run$values) < k) return(run)
  better <- function(a, b) switch(sel,
    LA = Re(a) > Re(b), SA = Re(a) < Re(b), LM = Mod(a) > Mod(b), SM = Mod(a) < Mod(b),
    LI = abs(Im(a)) > abs(Im(b)), SI = abs(Im(a)) < abs(Im(b)))
  found <- run$vectors
  ## A check may spend what the search spent and no more. At the small end of a
  ## spectrum the k = 1 request can stall on the residual row's plateau, and
  ## uncapped it spent the whole remaining budget there on OpenBLAS (300 of 300
  ## applies on laplacian_1d(60), against 120 with the reference BLAS); a check
  ## that does not converge swaps nothing, which is the answer when nothing was
  ## missed.
  cap <- run$steps
  for (s in seq_len(k)) {
    budget <- min(maxit - run$steps, cap)
    if (budget <= 0L) break
    cand <- core(k = 1L, maxit = budget, seed = seed + 9973L * s,
                 prelocked = qr.Q(qr(found)))
    run$steps <- run$steps + cand$steps
    run$rounds <- run$rounds + cand$rounds
    if (!is.null(cand$stop_reason) || !length(cand$values)) break
    worst <- ritz_order(run$thetas, sel)[k]
    margin <- 10 * tol * max(Mod(run$thetas[worst]), 1)
    if (!better(cand$thetas[1L], run$thetas[worst]) ||
        Mod(cand$thetas[1L] - run$thetas[worst]) <= margin) break
    found <- cbind(found, cand$vectors)
    if (is.complex(cand$values) && !is.complex(run$values)) {
      run$values <- as.complex(run$values); run$thetas <- as.complex(run$thetas)
    }
    if (is.complex(cand$vectors) && !is.complex(run$vectors)) storage.mode(run$vectors) <- "complex"
    run$values[worst] <- cand$values[1L]
    run$thetas[worst] <- cand$thetas[1L]
    run$residuals[worst] <- cand$residuals[1L]
    run$vectors[, worst] <- cand$vectors[, 1L]
    ord <- ritz_order(run$thetas, sel)
    run$values <- run$values[ord]; run$thetas <- run$thetas[ord]
    run$residuals <- run$residuals[ord]; run$vectors <- run$vectors[, ord, drop = FALSE]
    run$anorm_lb <- max(run$anorm_lb, cand$anorm_lb)
  }
  run
}


## One run of the recurrence, from a unit v orthogonal to the locked set Q and to
## the retained block Y, filling a basis of at most m vectors and spending at most
## max_applies of them on the operator. Each new direction A v_j is orthogonalised
## against Q and then against the whole of V, and the coefficients that second
## projection removes are kept as column j of the projected matrix,
##
##     T[1:j, j] = V[, 1:j]^H A v_j,
##
## mirrored to keep it hermitian. In exact arithmetic that is Lanczos: the column
## is (0, ..., beta_{j-1}, alpha_j), and after a thick restart the retained block
## is diagonal with one arrow row coupling it to v_{p+1}. In floating point the
## retained Ritz vectors are formed as V Z at every restart and carry a rounding
## error each time, which an assumed arrow form treats as exact pairs; the
## measured column is what they are. On the smallest end of laplacian_1d(n),
## whose residual row needs a vector accurate to a few eps in the high-frequency
## directions, the assumed form plateaus at a residual of 10.8 eps ||A|| at
## n = 500 (170 restarts), 20.5 at 1000 and 30 at 2000, above the 8.4 the request
## allows at n = 500, and the measured one at 7.2, 13.6 and 26. The plateau still
## grows with the square root of the restart count, as independent roundings do,
## from a lower start; an exact eigenvector evaluates to 0.5 eps ||A|| on the same
## rows. The reorthogonalisation computes the coefficients either way, so keeping
## them costs no apply.
lanczos_run <- function(Op, v, Q, m, max_applies, keep = NULL, B = NULL, BQ = NULL) {
  n <- nrow(v)
  eps <- .Machine$double.eps
  metric <- !is.null(B)
  p <- if (is.null(keep)) 0L else ncol(keep$Y)
  cplx <- is.complex(v) || (!is.null(keep) && is.complex(keep$Y)) ||
          Op$dtype == "complex" || (metric && B$dtype == "complex")
  zero <- if (cplx) complex(1) else 0
  V <- matrix(zero, n, m)
  Tm <- matrix(zero, m, m)
  ## In the B inner product of a generalized problem, B V beside V, so that every
  ## coefficient is (B v_i)^H w with one apply of B per step, as ARPACK does for
  ## bmat = "G".
  BV <- if (metric) matrix(zero, n, m)
  if (p > 0L) {
    V[, seq_len(p)] <- keep$Y
    if (metric) BV[, seq_len(p)] <- keep$BY
    Tm[seq_len(p), seq_len(p)] <- keep$H
  }
  V[, p + 1L] <- v
  if (metric) BV[, p + 1L] <- linop_apply(B, v, "N")
  bcols <- function(j) if (metric) BV[, j, drop = FALSE] else V[, j, drop = FALSE]

  i <- p
  applies <- 0L
  apply_norm <- 0
  beta_last <- 0
  vnext <- NULL

  repeat {
    i <- i + 1L
    applies <- applies + 1L
    W <- linop_apply(Op, V[, i, drop = FALSE], "N")
    ## The size the direction had before anything was taken out of it, which is
    ## the scale the breakdown test is relative to. The same test GMRES makes,
    ## with the same constant, and for the same reason: what it guards is a
    ## division, and normalising by a norm at rounding level would put an
    ## overflowed vector into the basis.
    w0 <- col_norms(W)
    apply_norm <- max(apply_norm, w0)
    ## The three-term part is subtracted first, so that the projection against
    ## the whole basis removes only what rounding and the retained block left;
    ## what it removes is added back into the column.
    a <- col_dot(bcols(i), W)
    W <- W - scale_cols(V[, i, drop = FALSE], a)
    if (i > p + 1L) W <- W - scale_cols(V[, i - 1L, drop = FALSE], b)
    W <- orth_against(Q, W, BQ)
    h <- orth_project(V[, seq_len(i), drop = FALSE], W,
                      if (metric) BV[, seq_len(i), drop = FALSE])
    W <- h$W
    h <- h$h[, 1L]
    h[i] <- h[i] + a
    if (i > p + 1L) h[i - 1L] <- h[i - 1L] + b
    Tm[seq_len(i), i] <- h
    Tm[i, seq_len(i)] <- Conj(h)
    ## A Rayleigh quotient of a hermitian operator is real; the imaginary part of
    ## the computed one is rounding.
    Tm[i, i] <- Re(h[i])

    b <- col_norms(W)
    ## An invariant subspace, which is convergence rather than failure: every
    ## Ritz pair of it is an exact eigenpair of the operator restricted to it.
    ## There is no direction left to restart from either, and vnext says so.
    if (b <= eps * w0) { beta_last <- b; break }
    if (metric) {
      BW <- linop_apply(B, W, "N")
      b <- sqrt(col_dot(W, BW))
    }
    beta_last <- b
    if (i >= m || applies >= max_applies) { vnext <- W / b; break }
    V[, i + 1L] <- W / b
    if (metric) BV[, i + 1L] <- BW / b
  }

  Tm <- Tm[seq_len(i), seq_len(i), drop = FALSE]
  if (!cplx) Tm <- Re(Tm)
  list(V = V[, seq_len(i), drop = FALSE],
       BV = if (metric) BV[, seq_len(i), drop = FALSE], Tm = Tm,
       size = i, applies = applies, vnext = vnext, beta_last = beta_last,
       apply_norm = apply_norm)
}
