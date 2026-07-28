## Section 7.1, the sixth of the seven Krylov methods. LSMR, on the same
## rectangular least-squares problem LSQR takes, Fong and Saunders 2011.
##
## Same bidiagonalisation, different minimisation. LSQR minimises ||b - A x||
## over the Krylov space, which is CG applied to the normal equations; LSMR
## minimises ||A^H (b - A x)||, which is MINRES applied to them. Both quantities
## converge to the same x, and the two methods differ in which one falls
## monotonically on the way.
##
## That difference is why this method is here rather than being a variant of the
## last one. The certificate reports ||A^H r|| / (||A|| ||r||) as the
## least-squares backward error, because Stewart's perturbation makes it one, and
## LSMR is the method that minimises the numerator of exactly that quantity at
## every step. Its stopping estimate is therefore monotone in the number the
## certificate will report, and LSQR's is not: LSQR's ||A^H r|| can rise from one
## step to the next while ||r|| falls, so a budget spent on LSQR can end on a
## worse backward error than a step earlier reached.
##
## Two QR factorisations rather than one. The first is LSQR's, reducing the lower
## bidiagonal B_k to an upper bidiagonal R_k with diagonal rho and superdiagonal
## theta. Applying MINRES to the normal equations means minimising over the same
## space the norm of
##
##     B_k^H (beta_1 e_1 - B_k y) = alpha_1 beta_1 e_1 - B_k^H B_k y,
##
## and B_k^H B_k = R_k^H R_k, so the projected problem is a least-squares problem
## in R_k^H, which is lower bidiagonal. A second sequence of rotations reduces it,
## giving the bar-quantities below. zetabar is that problem's residual, and it is
## ||A^H r_k|| exactly.
##
## Two direction vectors rather than one. LSQR carries w, the columns of
## V_k R_k^-1; LSMR carries h, the same thing, and hbar, its image under the
## second factorisation, and the iterate moves along hbar.
##
## The residual estimate is Fong and Saunders section 3.2. ||r_k|| is not a
## quantity the second projected problem carries, so it is recovered from a third
## rotation sequence applied to the same scalars; the alternative is an extra
## apply per step to measure it, which is what the loop in solvers-bidiag.R does
## once per round rather than once per step.
##
## Every scalar here is real and non-negative for the same reason as in LSQR:
## alpha and beta are norms, so B_k is real, and both factorisations of a real
## matrix are real.

#' Solve min ||b - A x|| by LSMR
#'
#' @param A A `linop` of any shape, square or rectangular. Both modes `"N"` and
#'   `"C"` are applied, so the operator has to supply an adjoint.
#' @param b A vector or block of right-hand sides.
#' @param preconditioner A `preconditioner` with `side = "right"`, or `NULL`. The
#'   iteration runs on `A M^-1` and recovers `x = M^-1 y`, so `M^-H` is applied
#'   as well as `M^-1`.
#' @param tol Relative tolerance, applied to `||b - A x|| / ||b||` where the
#'   system is compatible and to `||A^H r|| / (||A|| ||r||)` where it is not.
#'   Both are backward errors; see [solve_certificate()].
#' @param maxit Iteration budget. Defaults to `10 * n`.
#' @param x0 Starting iterate, or `NULL` for zero.
#' @param history Record every column's estimated residual at every iteration.
#'   Off by default.
#' @param conlim Stop a column once the estimated condition of the projected
#'   bidiagonal exceeds this, measured on the same triangular factor LSQR uses.
#' @param floor_const `c` in the arithmetic floor of the certificate.
#' @param norm_control Arguments for the `||A||` estimate.
#' @return A list with the solution, its certificate, and the iteration record.
#' @noRd
lsmr_solve <- function(A, b, preconditioner = NULL, tol = 1e-8, maxit = NULL,
                       x0 = NULL, history = FALSE,
                       conlim = KRYLOV_CONDITION_LIMIT,
                       floor_const = SOLVE_FLOOR_CONST, norm_control = list()) {
  bidiag_solve(A, b, "lsmr", lsmr_recurrence, preconditioner, tol, maxit,
               x0, history, conlim, floor_const, norm_control)
}

lsmr_recurrence <- function(A, X, R, active, target, tol, apply_precond,
                            apply_precond_adj, maxit, history, conlim,
                            anorm_lb) {
  eps <- .Machine$double.eps
  nc <- ncol(X)
  nothing <- function() list(X = X, iterations = 0L, anorm_lb = anorm_lb,
                             limited = FALSE, history = list())

  start <- bidiag_start(A, R, active, apply_precond_adj)
  if (is.null(start)) return(nothing())
  active <- start$active
  U <- start$U; V <- start$V
  alpha <- start$alpha; beta <- start$beta

  ## A real operator under a complex preconditioner reaches here with a complex
  ## direction and a real iterate waiting for it, and assigning one into the
  ## other would discard the imaginary part without saying so.
  if (is.complex(V) && !is.complex(X)) storage.mode(X) <- "complex"

  ## The two directions and the update. D is accumulated in the space the
  ## recurrence works in, which with a right preconditioner is the space of
  ## y = M x, so M^-1 is applied once to the assembled update rather than once
  ## per direction.
  H <- V
  Hbar <- zero_block(V)
  D <- zero_block(V)

  ## The first factorisation carries alphabar, the diagonal the next rotation
  ## eliminates against; the second carries cbar, sbar and rhobar, and zetabar,
  ## which is ||A^H r_k||. rho and rhobar are seeded at 1 so the first step's
  ## hbar update reduces to hbar = h, which is what a factorisation of nothing
  ## leaves behind.
  alphabar <- alpha
  zetabar <- alpha * beta
  zeta <- rep(0, length(active))
  rho <- rep(1, length(active))
  rhobar <- rep(1, length(active))
  cbar <- rep(1, length(active))
  sbar <- rep(0, length(active))

  ## The third sequence, for ||r_k||. Fong and Saunders section 3.2.
  betadd <- beta
  betad <- rep(0, length(active))
  rhodold <- rep(1, length(active))
  tautildeold <- rep(0, length(active))
  thetatilde <- rep(0, length(active))

  ## Extreme diagonals of the first triangular factor, per column, for the
  ## condition estimate. Seeded so the first step cannot trip it.
  dmax <- rep(0, length(active))
  dmin <- rep(Inf, length(active))

  hist <- list()
  it <- 0L
  limited <- FALSE

  ## Everything carried per column narrows here, and nothing narrows anywhere
  ## else: a vector left out of this list keeps the width the round started with
  ## and is silently recycled against the ones that did narrow.
  shrink <- function(keep) {
    active <<- active[keep]
    U <<- U[, keep, drop = FALSE]
    V <<- V[, keep, drop = FALSE]
    H <<- H[, keep, drop = FALSE]
    Hbar <<- Hbar[, keep, drop = FALSE]
    D <<- D[, keep, drop = FALSE]
    alpha <<- alpha[keep]; beta <<- beta[keep]
    alphabar <<- alphabar[keep]; zetabar <<- zetabar[keep]; zeta <<- zeta[keep]
    rho <<- rho[keep]; rhobar <<- rhobar[keep]
    cbar <<- cbar[keep]; sbar <<- sbar[keep]
    betadd <<- betadd[keep]; betad <<- betad[keep]
    rhodold <<- rhodold[keep]; tautildeold <<- tautildeold[keep]
    thetatilde <<- thetatilde[keep]
    dmax <<- dmax[keep]; dmin <<- dmin[keep]
  }

  retire <- function(done) {
    cols <- active[done]
    X[, cols] <<- X[, cols, drop = FALSE] + apply_precond(D[, done, drop = FALSE])
  }

  while (it < maxit && length(active)) {
    it <- it + 1L

    gk <- bidiag_step(A, U, V, alpha, apply_precond, apply_precond_adj, anorm_lb)
    U <- gk$U; Vt <- gk$Vt; beta <- gk$beta
    alpha_new <- gk$alpha; anorm_lb <- gk$anorm_lb

    ## ---------------------------------------------- the first factorisation --
    rhoold <- rho
    rho <- sqrt(alphabar^2 + beta^2)
    ill <- it > 1L & (pmin(dmin, rho) <= 0 | pmax(dmax, rho) > conlim * pmin(dmin, rho))
    if (any(ill)) {
      ## The step about to be taken divides by this rho, so a column whose
      ## projected problem has gone singular keeps the iterate it had before the
      ## step rather than the amplification of one.
      limited <- TRUE
      retire(which(ill))
      if (all(ill)) { active <- active[FALSE]; break }
      keep <- !ill
      shrink(keep)
      Vt <- Vt[, keep, drop = FALSE]; alpha_new <- alpha_new[keep]
      rho <- rho[keep]; rhoold <- rhoold[keep]
    }
    dmax <- pmax(dmax, rho)
    dmin <- pmin(dmin, rho)

    cs <- alphabar / rho
    sn <- beta / rho
    thetanew <- sn * alpha_new
    alphabar <- cs * alpha_new

    ## --------------------------------------------- the second factorisation --
    ## thetabar and rhotemp are formed from the rotation the previous step left
    ## behind, so both are read before cbar and sbar are replaced.
    rhobarold <- rhobar
    zetaold <- zeta
    thetabar <- sbar * rho
    rhotemp <- cbar * rho
    rhobar <- sqrt(rhotemp^2 + thetanew^2)
    cbar <- rhotemp / rhobar
    sbar <- thetanew / rhobar
    zeta <- cbar * zetabar
    zetabar <- -sbar * zetabar

    ## ------------------------------------------------- directions and update --
    Hbar <- H - scale_cols(Hbar, thetabar * rho / (rhoold * rhobarold))
    D <- D + scale_cols(Hbar, zeta / (rho * rhobar))
    ## alpha_{j+1} = 0 is the exhausted Krylov space: A^H u_{j+1} already lies in
    ## the span built so far, so there is no new direction. The zero column that
    ## leaves here is never read again, and dividing by alpha to produce it would
    ## put NaN into the iterate instead.
    Vn <- scale_cols(Vt, ifelse(alpha_new > 0, 1 / alpha_new, 0))
    H <- Vn - scale_cols(H, thetanew / rho)
    V <- Vn
    alpha <- alpha_new

    ## ----------------------------------------------------- the ||r|| estimate --
    ## A third rotation sequence, applied to the scalars the other two produced.
    ## Nothing here touches a vector, so the estimate costs no apply.
    betahat <- cs * betadd
    betadd <- -sn * betadd
    thetatildeold <- thetatilde
    rhotildeold <- sqrt(rhodold^2 + thetabar^2)
    ctildeold <- rhodold / rhotildeold
    stildeold <- thetabar / rhotildeold
    thetatilde <- stildeold * rhobar
    rhodold <- ctildeold * rhobar
    betad <- -stildeold * betad + ctildeold * betahat
    tautildeold <- (zetaold - thetatildeold * tautildeold) / rhotildeold
    taud <- (zeta - thetatilde * tautildeold) / rhodold

    ## ||r_j|| from the third sequence, and ||A^H r_j|| = |zetabar_{j+1}| from
    ## the second, both exact for the projected problem and both drifting from
    ## the true quantities as the recurrence loses orthogonality, which is why
    ## the outer loop re-measures rather than trusting them.
    rest <- sqrt((betad - taud)^2 + betadd^2)
    atest <- abs(zetabar)

    if (history) {
      row <- rep(NA_real_, nc)
      row[active] <- rest
      hist[[length(hist) + 1L]] <- row
    }

    ## A column leaves on whichever test applies to it. The optimality test uses
    ## the running lower bound on ||A||, which can only make it harder to meet.
    done <- rest <= target[active] |
            (rest > 0 & atest <= tol * anorm_lb * rest) |
            alpha <= eps * anorm_lb
    if (any(done)) {
      retire(which(done))
      if (all(done)) { active <- active[FALSE]; break }
      shrink(!done)
    }
  }

  if (length(active)) retire(seq_along(active))
  list(X = X, iterations = it, anorm_lb = anorm_lb, limited = limited,
       history = hist)
}
