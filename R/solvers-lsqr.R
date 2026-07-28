## Section 7.1, the fifth of the seven Krylov methods. LSQR, on the rectangular
## least-squares problem min ||b - A x||, Paige and Saunders 1982.
##
## This is the first method in the package whose problem is a minimisation rather
## than an equation, and almost everything that distinguishes it follows from
## that. The iterate lives in the domain and the residual in the codomain, which
## for a rectangular A are different spaces; b need not be in the range of A, so
## b - A x does not go to zero and its size is a fact about the problem rather
## than a measure of the solve; and what does go to zero is A^H r, which the
## certificate reads as a backward error through Stewart's perturbation.
##
## The bidiagonalisation is in solvers-bidiag.R, shared with LSMR. What is here
## is what LSQR does with it: the iterate x_k = x_0 + V_k y_k minimising
## ||b - A x|| over the Krylov space solves
##
##     min_y || beta_1 e_1 - B_k y ||,
##
## which plane rotations reduce as they do in MINRES and GMRES. The difference
## from both is that no basis is stored: the QR of a bidiagonal is a three-term
## recurrence, so x, one search direction w and the two current bidiagonal
## vectors are the whole state.
##
## Several right-hand sides run in lockstep, as in CG, MINRES and GMRES. Each
## column carries its own bidiagonal scalars, its own rotation and its own search
## direction, so the iterates are exactly those of per-column LSQR at two block
## applies per step instead of 2k.

#' Solve min ||b - A x|| by LSQR
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
#'   bidiagonal exceeds this. The singular values of the bidiagonal lie inside
#'   those of `A`, so the estimate is a lower bound on `cond(A)` and the test
#'   stops later than a sharp one would, never earlier.
#' @param floor_const `c` in the arithmetic floor of the certificate.
#' @param norm_control Arguments for the `||A||` estimate.
#' @return A list with the solution, its certificate, and the iteration record.
#' @noRd
lsqr_solve <- function(A, b, preconditioner = NULL, tol = 1e-8, maxit = NULL,
                       x0 = NULL, history = FALSE,
                       conlim = KRYLOV_CONDITION_LIMIT,
                       floor_const = SOLVE_FLOOR_CONST, norm_control = list()) {
  bidiag_solve(A, b, "lsqr", lsqr_recurrence, preconditioner, tol, maxit,
               x0, history, conlim, floor_const, norm_control)
}

## The QR of the bidiagonal, taken one rotation at a time. phibar carries the
## residual of the projected problem and rhobar the diagonal the next rotation
## will eliminate against.
lsqr_recurrence <- function(A, X, R, active, target, tol, apply_precond,
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
  W <- V

  ## A real operator under a complex preconditioner reaches here with a complex
  ## direction and a real iterate waiting for it, and assigning one into the
  ## other would discard the imaginary part without saying so.
  if (is.complex(V) && !is.complex(X)) storage.mode(X) <- "complex"

  ## The update, accumulated in the space the recurrence works in. With a right
  ## preconditioner that is the space of y = M x, so M^-1 is applied once to the
  ## assembled update rather than once per direction.
  D <- zero_block(V)

  phibar <- beta
  rhobar <- alpha
  ## Extreme diagonals of the triangular factor, per column, for the condition
  ## estimate. Seeded so the first step cannot trip it.
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
    W <<- W[, keep, drop = FALSE]
    D <<- D[, keep, drop = FALSE]
    alpha <<- alpha[keep]; beta <<- beta[keep]
    phibar <<- phibar[keep]; rhobar <<- rhobar[keep]
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

    ## ------------------------------------------------------------ rotation --
    ## Eliminates beta_{j+1} against the current diagonal. Both are norms, so the
    ## rotation is real whatever the operator, and so is every scalar below it.
    rho <- sqrt(rhobar^2 + beta^2)
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
      rho <- rho[keep]
    }
    dmax <- pmax(dmax, rho)
    dmin <- pmin(dmin, rho)

    cs <- rhobar / rho
    sn <- beta / rho
    theta <- sn * alpha_new
    rhobar <- -cs * alpha_new
    phi <- cs * phibar
    phibar <- sn * phibar

    ## alpha_{j+1} = 0 is the exhausted Krylov space: A^H u_{j+1} already lies in
    ## the span built so far, so there is no new direction. The zero column that
    ## leaves here is never read again, and dividing by alpha to produce it would
    ## put NaN into the iterate instead.
    Vn <- scale_cols(Vt, ifelse(alpha_new > 0, 1 / alpha_new, 0))
    D <- D + scale_cols(W, phi / rho)
    W <- Vn - scale_cols(W, theta / rho)
    V <- Vn
    alpha <- alpha_new

    ## --------------------------------------------------------- the estimates --
    ## ||r_j|| = phibar_{j+1} and ||A^H r_j|| = phibar_{j+1} alpha_{j+1} |c_j|,
    ## both exact for the projected problem and both drifting from the true
    ## quantities as the recurrence loses orthogonality, which is why the outer
    ## loop re-measures rather than trusting them.
    rest <- phibar
    atest <- phibar * alpha * abs(cs)

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
