## Section 7.1, the second of the seven Krylov methods. Preconditioned MINRES.
##
## CG needs positive definiteness. MINRES needs only that A is hermitian, which
## makes it the method the three-valued capability set of section 5.3 routes to
## when positive_definite is FALSE, and equally when it is NA: unknown is not
## false, and an operator nobody has proved definite is a MINRES problem rather
## than a refusal.
##
## The preconditioned recurrence. For a hermitian positive definite M = L L^H,
## MINRES runs on the split-preconditioned operator L^-1 A L^-H, hermitian
## because A is, where
##
##     || L^-1 (b - A x) ||_2^2 = (b - A x)^H M^-1 (b - A x),
##
## so what the method minimises is the M^-1 norm of the true residual. Writing
## v_j = L^-H q_j for the Lanczos vectors q_j of the split operator turns the
## recurrence into one that never forms L:
##
##     M^-1 A v_j = beta_{j+1} v_{j+1} + alpha_j v_j + beta_j v_{j-1},
##
## with the v_j orthonormal in the M inner product, alpha_j = v_j^H A v_j real,
## and beta_{j+1} = ||w_j||_{M^-1} for w_j = A v_j - alpha_j M v_j - beta_j M v_{j-1}.
## Carrying u_j = M v_j alongside v_j keeps M itself out of it, so only M^-1 is
## ever applied, and the update x_m = x_0 + V_m y_m uses the v_j directly.
##
## The least-squares problem min ||beta_1 e_1 - T_m y|| is reduced by plane
## rotations, and the solution is carried in direction vectors d_j = column j of
## V_m R_m^-1 rather than by storing V_m, which is what makes the method a short
## recurrence with fixed storage.
##
## Where this departs from CG. CG's preconditioned iterates coincide on all three
## sides and its recurrence residual is measured in the same norm the certificate
## reports, so nothing downstream depends on the preconditioner. MINRES minimises
## ||r||_{M^-1}, which is not ||r||_2 unless M is the identity. The quantity the
## recurrence tracks and the quantity the certificate reports are therefore
## different quantities whenever a preconditioner is present, and the outer loop
## is what reconciles them: the inner loop runs to a target converted into its own
## currency by the ratio measured at the start of the round, and the decision to
## stop is taken on a recomputed b - A x in the norm the certificate reports.
##
## Several right-hand sides run in lockstep, as in CG. Each column carries its own
## Lanczos scalars, its own rotations and its own direction vectors, so the
## iterates are exactly those of per-column MINRES, at one block apply per step
## instead of k.

## Plan 7.1: MINRES requires hermitian at declared minimum evidence, on the same
## reading of the evidence fields as CG's positive-definiteness requirement. An
## exact check on data the operator already holds is an identity rather than a
## probe, and probes are excluded by the guarantee field rather than by the source
## list. This filters method = "auto" and does not gate method = "minres".
MINRES_HERMITIAN_REQUIREMENT <- requirement(
  sources = c("construction", "adapter_contract", "theorem", "computation"),
  guarantees = "identity", min_confidence = 1)

## Relative slack allowed in the identities a hermitian operator satisfies before
## the declaration is treated as contradicted. Unlike CG's p^H A p <= 0 this is a
## threshold rather than a sign, because non-hermitian is not a one-sided
## condition, so where it sits was measured rather than chosen: over the ill
## conditioned, clustered and near-singular hermitian fixtures the worst relative
## violation observed is 5.3e-10, and a skew part of relative size 1e-6 produces
## 1.3e-06. The default sits three orders above the observed noise and catches
## every asymmetry from there upward.
##
## Below that an operator is hermitian to within 1e-6 and the contradiction is
## not the right instrument: what happens instead is that the solve fails to
## converge and the certificate reports the true residual it actually reached.
MINRES_SYMMETRY_TOL <- 1e-6

#' Solve A x = b by preconditioned MINRES
#'
#' @param A A square hermitian `linop`, definite or not.
#' @param b A vector or block of right-hand sides.
#' @param preconditioner A `preconditioner`, or `NULL`. Hermitian positive
#'   definite even though `A` is not; the method minimises the `M^-1` norm of the
#'   residual and that is a norm only for a definite `M`.
#' @param tol Relative residual tolerance, `||b - A x|| <= tol * ||b||`, in the
#'   euclidean norm whatever the preconditioner.
#' @param maxit Iteration budget. Defaults to `10 * n`.
#' @param x0 Starting iterate, or `NULL` for zero.
#' @param history Record every column's recurrence residual at every iteration.
#'   That quantity is `||r||_{M^-1}`, which equals `||r||_2` only without a
#'   preconditioner. Off by default.
#' @param floor_const `c` in the arithmetic floor of the certificate.
#' @param norm_control Arguments for the `||A||` estimate.
#' @param symmetry_tol Relative slack in the hermitian identities before the
#'   declaration is reported as contradicted.
#' @return A list with the solution, its certificate, and the iteration record.
#' @noRd
minres_solve <- function(A, b, preconditioner = NULL, tol = 1e-8, maxit = NULL,
                         x0 = NULL, history = FALSE,
                         floor_const = SOLVE_FLOOR_CONST, norm_control = list(),
                         symmetry_tol = MINRES_SYMMETRY_TOL) {
  if (!is_linop(A)) stopf("minres() expects a linop")
  n <- A$dim[2L]
  if (A$dim[1L] != n) {
    stopf(paste0("minres() needs a square operator; this one is %d x %d.\n",
                 "  A rectangular system is a least-squares problem and takes a different method."),
          A$dim[1L], A$dim[2L])
  }
  require_capability(A, "hermitian", "minres")
  check_preconditioner(preconditioner, "minres")

  was_vector <- is.null(dim(b))
  B <- as_block(b)
  if (nrow(B) != n) {
    stopf("non-conformable: operator is %d x %d, right-hand side has %d rows",
          n, n, nrow(B))
  }
  k <- ncol(B)
  maxit <- as.integer(maxit %||% min(10 * as.numeric(n), .Machine$integer.max))
  if (is.na(maxit) || maxit < 1L) stopf("maxit must be a positive integer")

  X <- if (is.null(x0)) {
    matrix(0, n, k)
  } else {
    x <- as_block(x0)
    if (!identical(dim(x), c(n, k))) {
      stopf("x0 is %d x %d; the right-hand side is %d x %d", nrow(x), ncol(x), n, k)
    }
    x
  }
  if (A$dtype == "complex" || is.complex(B) || is.complex(X)) {
    storage.mode(X) <- "complex"
  }

  apply_precond <- precond_applier(preconditioner)

  bn <- col_norms(B)
  ## A zero right-hand side has the exact solution 0, where a relative test is
  ## undefined; that column gets an absolute one.
  target <- tol * ifelse(bn > 0, bn, 1)

  iterations <- 0L
  rounds <- 0L
  krylov_lb <- 0
  hist <- list()
  prev_rn <- NULL

  repeat {
    ## The outer loop measures the true residual in the norm the certificate
    ## reports; the inner loop trusts its own recurrence in the norm it minimises.
    R <- B - linop_apply(A, X, "N")
    rn <- col_norms(R)
    if (all(rn <= target)) break
    if (iterations >= maxit) break
    ## A restart that did not reduce the true residual has nothing left to
    ## recover. Restarting costs more here than it does in CG, because it discards
    ## the Krylov space a minimal-residual method minimises over, which is the
    ## reason the inner loop returns only when it believes the target is met.
    if (!is.null(prev_rn) && all(rn >= prev_rn)) break
    prev_rn <- rn

    step <- minres_recurrence(A, X, R, rn, target, apply_precond,
                              maxit - iterations, history, symmetry_tol)
    X <- step$X
    iterations <- iterations + step$iterations
    krylov_lb <- max(krylov_lb, step$krylov_lb)
    if (history) hist <- c(hist, step$history)
    rounds <- rounds + 1L
    if (step$iterations == 0L) break
  }

  norm_est <- do.call(norm2, c(list(A = A), norm_control))
  cert <- solve_certificate(A, B, X, tol = tol, norm_estimate = norm_est,
                            iterations = iterations, maxit = maxit,
                            floor_const = floor_const,
                            norm_lower_bound = krylov_lb)

  list(x = undo_block(X, was_vector),
       certificate = cert,
       method = "minres",
       iterations = iterations,
       restarts = max(0L, rounds - 1L),
       converged = cert_status(cert, "residual") != "fail",
       residual = cert$values$residual,
       history = if (history && length(hist)) do.call(rbind, hist) else NULL)
}

## One run of the recurrence, from the true residual R, until every column meets
## its target or the budget runs out. Converged columns leave the active set, so
## the block narrows as the solve proceeds.
##
## `rn` is the euclidean norm of R, measured by the caller. It is needed here only
## to convert the caller's target into the M^-1 currency the recurrence works in;
## without a preconditioner the conversion is the identity.
minres_recurrence <- function(A, X, R, rn, target, apply_precond, maxit, history,
                              symmetry_tol) {
  eps <- .Machine$double.eps
  nc <- ncol(X)
  active <- which(rn > target)
  Xa <- X[, active, drop = FALSE]
  Ra <- R[, active, drop = FALSE]

  ## beta_1 = ||r_0||_{M^-1}, v_1 = M^-1 r_0 / beta_1, u_1 = M v_1 = r_0 / beta_1.
  Z <- apply_precond(Ra)
  b1sq <- col_dot(Ra, Z)
  if (any(b1sq <= 0)) precond_not_hpd("<r, M^-1 r>")
  beta1 <- sqrt(b1sq)
  U <- scale_cols(Ra, 1 / beta1)
  V <- scale_cols(Z, 1 / beta1)
  U_old <- zero_block(U)
  V_old <- zero_block(V)
  P_old <- zero_block(V)
  D1 <- zero_block(U)
  D2 <- zero_block(U)

  ## The recurrence minimises ||r||_{M^-1} and the caller asked in ||r||_2. The
  ## ratio between the two is known exactly at the start of the round, for the
  ## residual the round starts from, so the target is carried across by it and the
  ## outer loop checks the result in the norm that was actually requested.
  inner_target <- beta1 * (target[active] / rn[active])

  phibar <- beta1
  beta <- rep(0, length(active))
  cs_p <- rep(1, length(active)); sn_p <- rep(0, length(active))
  cs_p2 <- cs_p; sn_p2 <- sn_p

  krylov_lb <- 0
  hist <- list()
  it <- 0L

  shrink <- function(keep) {
    active <<- active[keep]
    Xa <<- Xa[, keep, drop = FALSE]
    U <<- U[, keep, drop = FALSE]; U_old <<- U_old[, keep, drop = FALSE]
    V <<- V[, keep, drop = FALSE]; V_old <<- V_old[, keep, drop = FALSE]
    P_old <<- P_old[, keep, drop = FALSE]
    D1 <<- D1[, keep, drop = FALSE]; D2 <<- D2[, keep, drop = FALSE]
    cs_p <<- cs_p[keep]; sn_p <<- sn_p[keep]
    cs_p2 <<- cs_p2[keep]; sn_p2 <<- sn_p2[keep]
    phibar <<- phibar[keep]; beta <<- beta[keep]
    inner_target <<- inner_target[keep]
  }

  while (it < maxit && length(active)) {
    it <- it + 1L
    P <- linop_apply(A, V, "N")

    ## ||A v|| / ||v|| <= ||A||_2 for every v, so the iteration hands the
    ## certificate a lower bound on the norm at no extra apply.
    vn <- col_norms(V)
    if (any(vn > 0)) krylov_lb <- max(krylov_lb, max(col_norms(P)[vn > 0] / vn[vn > 0]))

    alpha <- col_dot(V, P)
    minres_check_hermitian(V, V_old, P, P_old, it, symmetry_tol)

    W <- P - scale_cols(U, alpha) - scale_cols(U_old, beta)
    Znew <- apply_precond(W)
    b2sq <- col_dot(W, Znew)

    ## w^H M^-1 w > 0 for every nonzero w when M is hermitian positive definite.
    ## A w that has collapsed to rounding level is the Krylov space running out,
    ## which is convergence rather than a contradiction, and is separated by size
    ## exactly as CG separates its benign case.
    spent <- col_norms(W) <= sqrt(eps) * col_norms(P)
    if (any(b2sq <= 0 & !spent)) precond_not_hpd("<w, M^-1 w>")
    beta_next <- sqrt(pmax(b2sq, 0))

    ## Two previous rotations act on this column of the tridiagonal, then a new
    ## one eliminates beta_{j+1} against the diagonal.
    eps_j <- sn_p2 * beta
    dbar <- cs_p2 * beta
    delta <- cs_p * dbar + sn_p * alpha
    gbar <- cs_p * alpha - sn_p * dbar
    gamma <- sqrt(gbar^2 + beta_next^2)

    ## gamma = 0 is gbar = 0 and beta_{j+1} = 0 together: the projected system is
    ## singular and no step can be taken. The column keeps what it has and leaves.
    singular <- gamma <= 0
    gsafe <- ifelse(singular, 1, gamma)
    cs <- ifelse(singular, 1, gbar / gsafe)
    sn <- ifelse(singular, 0, beta_next / gsafe)
    phi <- ifelse(singular, 0, cs * phibar)
    phibar <- ifelse(singular, phibar, -sn * phibar)

    D <- scale_cols(V - scale_cols(D1, delta) - scale_cols(D2, eps_j), 1 / gsafe)
    Xa <- Xa + scale_cols(D, phi)

    D2 <- D1; D1 <- D
    U_old <- U; V_old <- V; P_old <- P
    ## beta_{j+1} = 0 leaves no direction to continue in, and the rotation above
    ## has already set phibar to zero there, so the column is converged and the
    ## zero vector it carries is never read again.
    inv <- ifelse(beta_next > 0, 1 / beta_next, 0)
    U <- scale_cols(W, inv)
    V <- scale_cols(Znew, inv)

    cs_p2 <- cs_p; sn_p2 <- sn_p
    cs_p <- cs; sn_p <- sn
    beta <- beta_next

    if (history) {
      row <- rep(NA_real_, nc)
      row[active] <- abs(phibar)
      hist[[length(hist) + 1L]] <- row
    }

    done <- abs(phibar) <= inner_target | singular
    if (any(done)) {
      X[, active[done]] <- Xa[, done, drop = FALSE]
      if (all(done)) { active <- active[FALSE]; break }
      shrink(!done)
    }
  }

  if (length(active)) X[, active] <- Xa
  list(X = X, iterations = it, krylov_lb = krylov_lb, history = hist)
}

## Is A actually hermitian? The identity to test is the definition, <x, A y> =
## <A x, y>, on the two vectors the iteration already has: A v_j is this step's
## apply and A v_{j-1} was the previous one, so the test costs no apply.
##
## The recurrence supplies a second candidate, v_{j-1}^H A v_j = beta_j, which is
## the wrong instrument. It holds only while the v_j remain orthogonal, and
## classical Lanczos loses orthogonality as Ritz values converge, so it reaches a
## relative violation of 1e0 on hermitian operators with clustered spectra and on
## most complex hermitian ones. It tests orthogonality, not symmetry. The adjoint
## identity does not depend on orthogonality at all and stays at 1e-10 across the
## same fixtures.
##
## The imaginary part of v^H A v is a third signal, independent of the first and
## equally free, and vacuous in real arithmetic.
minres_check_hermitian <- function(V, V_old, P, P_old, it, symmetry_tol) {
  scale <- col_norms(V) * col_norms(P)
  ok <- scale > 0
  if (!any(ok)) return(invisible(TRUE))
  bad <- ok & abs(col_dot_im(V, P)) > symmetry_tol * scale
  if (it >= 2L) {
    gap <- sqrt((col_dot(V_old, P) - col_dot(P_old, V))^2 +
                (col_dot_im(V_old, P) - col_dot_im(P_old, V))^2)
    cross <- col_norms(V_old) * col_norms(P)
    bad <- bad | (cross > 0 & gap > symmetry_tol * cross)
  }
  if (any(bad)) operator_not_hermitian(sum(bad))
  invisible(TRUE)
}

operator_not_hermitian <- function(ncols) {
  stopf(paste0("minres() found <x, A y> and <A x, y> disagreeing in %d column(s).\n",
               "  The operator declares hermitian = TRUE and the iteration contradicts it.\n",
               "  MINRES applies only to a hermitian operator: it minimises over a Krylov space\n",
               "  built by a three-term recurrence that a non-hermitian operator does not satisfy."),
        ncols)
}

precond_not_hpd <- function(what) {
  stopf(paste0("minres() reached %s <= 0.\n",
               "  The preconditioner declares positive_definite = TRUE and the iteration contradicts it.\n",
               "  MINRES minimises the M^-1 norm of the residual, which is a norm only for a\n",
               "  hermitian positive definite M, so definiteness is required even though A is indefinite."),
        what)
}
