## Section 7.1, the seventh of the seven Krylov methods. BiCGSTAB, van der Vorst
## 1992, on a general square system.
##
## Like GMRES it requires nothing of the operator and applies no adjoint: every
## apply below is mode "N". Unlike GMRES it stores no basis. The two methods sit
## at opposite ends of the same trade: GMRES minimises the residual over the whole
## Krylov space and pays for it in storage and orthogonalisation that grow with
## the round, and restarts to bound them; BiCGSTAB runs a short recurrence in
## fixed storage with two applies per step and minimises nothing, so its residual
## is not monotone and its recurrence can break down.
##
## The iteration. With a shadow vector rhat0 fixed for the round,
##
##     rho_j     = <rhat0, r_{j-1}>,
##     beta      = (rho_j / rho_{j-1}) (alpha / omega_{j-1}),
##     p_j       = r_{j-1} + beta (p_{j-1} - omega_{j-1} v_{j-1}),
##     v_j       = A p_j,
##     alpha     = rho_j / <rhat0, v_j>,
##     s         = r_{j-1} - alpha v_j,
##     t         = A s,
##     omega_j   = <t, s> / <t, t>,
##     x_j       = x_{j-1} + alpha p_j + omega_j s,
##     r_j       = s - omega_j t.
##
## The BiCG half is the alpha step and the stabilising half is the omega step,
## which is one step of GMRES(1) on the intermediate residual s: omega minimises
## ||s - omega t||, and that is the whole of the minimisation this method does.
##
## Complex arithmetic. rho, alpha, beta and omega are all genuinely complex for a
## complex operator, since none of the products above is one a hermitian
## recurrence forms, so this file uses col_cdot() throughout as GMRES does.
##
## Several right-hand sides run in lockstep, as in every other method here. Each
## column carries its own shadow vector, its own scalars and its own directions,
## so the iterates are exactly those of per-column BiCGSTAB at two block applies
## per step instead of 2k.

## ------------------------------------------------------------ the three sides

## Section 4.3 leaves BiCGSTAB unrestricted, and as with GMRES that is not a
## statement that the sides agree. Each is a different iteration on a different
## operator and they produce different iterates.
##
##   right   BiCGSTAB on A M^-1, which is Saad's Algorithm 7.7. r stays the true
##           b - A x, so the recurrence measures the quantity the certificate
##           reports and nothing has to be converted. M^-1 is applied to p and to
##           s, and the iterate moves along those two preconditioned directions.
##
##   left    BiCGSTAB on M^-1 A with right-hand side M^-1 b. The recurrence
##           residual is M^-1 (b - A x), so the stopping test needs the currency
##           conversion MINRES introduced. The iterate moves along p and s
##           themselves.
##
##   split   BiCGSTAB on L^-1 A L^-H for M = L L^H, reached without ever forming
##           L. Under the uniform change of variable w~ = L^-1 w every vector of
##           the iteration maps to one of the right-preconditioned iteration, and
##           the euclidean inner product maps to <u, M^-1 v>. So the split form is
##           the right form with one substitution: the inner product. It measures
##           sqrt(<r, M^-1 r>), which again needs the conversion.
##
## Split would cost three further preconditioner applies per step if <u, M^-1 v>
## were formed on demand. It costs none, because M^-1 r, M^-1 s and M^-1 p obey
## the same linear recurrences their unpreconditioned counterparts do, so carrying
## them alongside leaves only M^-1 v and M^-1 t to be applied. That is the trick
## GMRES uses to carry u_j = M v_j, one level up: two applies per step on every
## side.
##
## The split form is defined only for a hermitian positive definite M. The
## BiCGSTAB row of PRECOND_REQUIREMENTS asks only for `fixed`, so that is checked
## at run time on <r, M^-1 r>, through the same message GMRES uses.

## How small an inner product has to be, relative to the two vectors that formed
## it, before the scalar it is about to divide is noise rather than a number.
##
## The three breakdowns of BiCGSTAB are all of this shape. rho = <rhat0, r> going
## to zero is the BiCG breakdown, <rhat0, v> going to zero is the alpha
## breakdown, and <t, s> going to zero is the failure of the stabilising step,
## which leaves omega at zero and the p recurrence with a division by it.
##
## The threshold was measured rather than assumed, on operators where no
## breakdown is present. Over 12 seeds on six fixtures, the smallest relative
## value each of the three reached across a full solve was
##
##                     rho        <rhat0, v>   <t, s>
##   convdiff mu=.3    2.8e-05    3.4e-04      9.5e-02
##   convdiff mu=.7    1.2e-04    1.7e-05      9.1e-02
##   laplacian         1.0e-08    1.5e-08      8.8e-02
##   kms rho=.7        1.1e-05    1.3e-04      3.7e-01
##   shifted lap       8.8e-12    2.3e-11      1.6e-04
##   complex dense     7.0e-03    1.5e-02      2.8e-01
##
## so a healthy solve came within a factor of about 900 of 1e-14 at its worst, on
## a nearly singular indefinite operator, and stayed six orders above it
## everywhere else. 1e-14 is also about 45 eps, which is where a relative inner
## product stops carrying information at all, so the threshold is at the
## arithmetic floor rather than at a level tuned to the fixtures.
BICGSTAB_BREAKDOWN_TOL <- 1e-14

#' Solve A x = b by BiCGSTAB
#'
#' @param A A square `linop`. No capability is required, and no adjoint.
#' @param b A vector or block of right-hand sides.
#' @param preconditioner A `preconditioner`, or `NULL`. Its `side` selects the
#'   algorithm; see the note at the top of this file.
#' @param tol Relative residual tolerance, `||b - A x|| <= tol * ||b||`, in the
#'   euclidean norm whatever the preconditioner.
#' @param maxit Iteration budget. Defaults to `10 * n`. Each iteration is two
#'   applies of `A`.
#' @param x0 Starting iterate, or `NULL` for zero.
#' @param history Record every column's recurrence residual at every iteration.
#'   That quantity is the norm the chosen side works in, which is `||r||_2` only
#'   for right preconditioning. Off by default.
#' @param breakdown_tol Relative size below which an inner product about to
#'   become a denominator is treated as a breakdown.
#' @param floor_const `c` in the arithmetic floor of the certificate.
#' @param norm_control Arguments for the `||A||` estimate.
#' @return A list with the solution, its certificate, and the iteration record.
#' @noRd
bicgstab_solve <- function(A, b, preconditioner = NULL, tol = 1e-8, maxit = NULL,
                           x0 = NULL, history = FALSE,
                           breakdown_tol = BICGSTAB_BREAKDOWN_TOL,
                           floor_const = SOLVE_FLOOR_CONST, norm_control = list()) {
  s <- solver_setup(A, b, x0, maxit, "bicgstab")
  B <- s$B; X <- s$X; maxit <- s$maxit
  was_vector <- s$was_vector

  check_preconditioner(preconditioner, "bicgstab")
  apply_precond <- precond_applier(preconditioner)
  ## Without a preconditioner every side is the identity, and "right" is the one
  ## that converts nothing, so it is also the cheapest way to say "none".
  side <- if (is.null(preconditioner)) "right" else preconditioner$side
  pd_declared <- isTRUE(preconditioner$positive_definite)

  bn <- col_norms(B)
  ## A zero right-hand side has the exact solution 0, where a relative test is
  ## undefined; that column gets an absolute one.
  target <- tol * ifelse(bn > 0, bn, 1)

  iterations <- 0L
  rounds <- 0L
  krylov_lb <- 0
  broke <- FALSE
  hist <- list()
  prev_rn <- NULL

  repeat {
    ## The outer loop measures the true residual in the norm the certificate
    ## reports; the inner loop trusts its own recurrence in the norm it works in.
    R <- B - linop_apply(A, X, "N")
    rn <- col_norms(R)
    if (all(rn <= target)) break
    if (iterations >= maxit) break
    ## A round that did not reduce the true residual has nothing left to recover.
    ## This is also what ends a solve whose recurrence keeps breaking down: the
    ## restart below re-seeds the shadow vector, which is the standard cure, and
    ## a cure that recovers nothing twice running is reported rather than looped
    ## on.
    if (!is.null(prev_rn) && all(rn >= prev_rn)) break
    prev_rn <- rn

    step <- bicgstab_round(A, X, R, rn, target, apply_precond, side,
                           maxit - iterations, history, breakdown_tol,
                           pd_declared)
    X <- step$X
    iterations <- iterations + step$iterations
    krylov_lb <- max(krylov_lb, step$krylov_lb)
    broke <- broke || step$broke
    if (history) hist <- c(hist, step$history)
    rounds <- rounds + 1L
    if (step$iterations == 0L) break
  }

  norm_est <- do.call(norm2, c(list(A = A), norm_control))
  cert <- solve_certificate(A, B, X, tol = tol, norm_estimate = norm_est,
                            iterations = iterations, maxit = maxit,
                            floor_const = floor_const,
                            norm_lower_bound = krylov_lb,
                            stop_reason = if (broke && iterations < maxit)
                              "the recurrence broke down and a fresh shadow residual did not recover it"
                            else NULL)

  list(x = undo_block(X, was_vector),
       certificate = cert,
       method = "bicgstab",
       iterations = iterations,
       restarts = max(0L, rounds - 1L),
       converged = cert_status(cert, "residual") != "fail",
       residual = cert$values$residual,
       history = if (history && length(hist)) do.call(rbind, hist) else NULL)
}

## One round from the true residual R, of at most `steps` steps. Returns when
## every active column has met its target, broken down, or the round is spent.
##
## The shadow vector is the residual this round starts from. Any vector not
## orthogonal to r_0 will do, and this one makes rho_1 = <r_0, r_0> exactly the
## squared norm the round begins with, so the first step cannot break down; a
## breakdown later ends the round for that column and the loop above re-seeds
## from a residual that has moved.
bicgstab_round <- function(A, X, R, rn, target, apply_precond, side, steps,
                           history, breakdown_tol, pd_declared) {
  nc <- ncol(X)
  active <- which(rn > target)
  Xa <- X[, active, drop = FALSE]
  Ra <- R[, active, drop = FALSE]

  ## The residual this side works in, the norm it measures, and for the split
  ## form the preconditioned residual carried alongside it.
  ZR <- NULL
  if (side == "left") {
    Rr <- apply_precond(Ra)
    beta0 <- col_norms(Rr)
    if (any(beta0 <= 0)) left_precond_singular("bicgstab")
  } else if (side == "split") {
    Rr <- Ra
    ZR <- apply_precond(Ra)
    q <- col_dot(Ra, ZR)
    if (any(q <= 0)) split_precond_not_hpd("bicgstab", "<r, M^-1 r>", pd_declared)
    beta0 <- sqrt(q)
  } else {
    Rr <- Ra
    beta0 <- rn[active]
  }

  ## The recurrence works in its own norm and the caller asked in the euclidean
  ## one. The ratio is exact for the residual this round starts from, and the
  ## outer loop re-measures afterwards, which is what makes an approximate
  ## conversion safe rather than what makes it right.
  inner_target <- beta0 * (target[active] / rn[active])

  ## <u, M^-1 v> for the split form and <u, v> for the other two, with the
  ## preconditioned second argument supplied by the caller because it is always
  ## either carried or already needed for an apply.
  ip <- if (side == "split") function(U, V, ZV) col_cdot(U, ZV)
        else function(U, V, ZV) col_cdot(U, V)

  Rhat <- Rr
  P <- Rr
  ZP <- ZR
  rho <- ip(Rhat, Rr, ZR)

  ## A complex preconditioner over a real operator: the scalars below are complex
  ## whenever any vector is, and a real iterate waiting for a complex update
  ## would take its real part without saying so. Both blocks have to move, since
  ## the retired columns are written back into X.
  if (is.complex(Rr) && !is.complex(Xa)) {
    storage.mode(Xa) <- "complex"
    storage.mode(X) <- "complex"
  }

  krylov_lb <- 0
  hist <- list()
  it <- 0L
  broke <- FALSE

  ## Everything carried per column narrows here, and nothing narrows anywhere
  ## else: a vector left out of this list keeps the width the round started with
  ## and is silently recycled against the ones that did narrow.
  shrink <- function(keep) {
    active <<- active[keep]
    Xa <<- Xa[, keep, drop = FALSE]
    Rr <<- Rr[, keep, drop = FALSE]
    Rhat <<- Rhat[, keep, drop = FALSE]
    P <<- P[, keep, drop = FALSE]
    if (!is.null(ZR)) ZR <<- ZR[, keep, drop = FALSE]
    if (!is.null(ZP)) ZP <<- ZP[, keep, drop = FALSE]
    rho <<- rho[keep]
    inner_target <<- inner_target[keep]
  }

  retire <- function(done) {
    X[, active[done]] <<- Xa[, done, drop = FALSE]
  }

  ## The norm this side measures, on the residual it carries.
  measure <- function(Rr, ZR) {
    if (side == "split") sqrt(pmax(col_dot(Rr, ZR), 0)) else col_norms(Rr)
  }

  ## One history row per iteration, opened when the iteration starts and written
  ## into wherever a column's recurrence residual becomes known. A column that
  ## leaves at the early exit below has its value known half a step before the
  ## ones that stay, and a step that ends there would otherwise leave no row at
  ## all for an iteration that happened.
  open_row <- function() hist[[length(hist) + 1L]] <<- rep(NA_real_, nc)
  note <- function(cols, vals) hist[[length(hist)]][cols] <<- vals

  ## An inner product that is about to become a denominator, against the two
  ## vectors that formed it.
  degenerate <- function(v, U, V) Mod(v) <= breakdown_tol * col_norms(U) * col_norms(V)

  ## A column whose recurrence has broken down keeps the iterate it had before
  ## the step rather than the amplification of a division by noise, and leaves.
  ## The caller re-measures and starts again from a residual that has moved.
  drop_broken <- function(bad) {
    broke <<- TRUE
    retire(which(bad))
    if (all(bad)) { active <<- active[FALSE]; return(TRUE) }
    shrink(!bad)
    FALSE
  }

  while (it < steps && length(active)) {
    it <- it + 1L
    if (history) open_row()

    ## ------------------------------------------------------ the BiCG half --
    ## p_j, and the direction the iterate moves along, which is p itself for the
    ## left form and M^-1 p for the other two.
    if (side == "left") {
      ZPj <- P
      V <- linop_apply(A, P, "N")
      lb_src <- P; lb_img <- V
      V <- apply_precond(V)
      ZV <- NULL
    } else {
      ZPj <- if (side == "split") ZP else apply_precond(P)
      V <- linop_apply(A, ZPj, "N")
      lb_src <- ZPj; lb_img <- V
      ZV <- if (side == "split") apply_precond(V) else NULL
    }
    ## ||A z|| / ||z|| <= ||A||_2 for every z, so the iteration hands the
    ## certificate a lower bound on the norm at no extra apply. Measured on the
    ## raw operator, before any preconditioner touches the result.
    zn <- col_norms(lb_src)
    if (any(zn > 0)) {
      an <- col_norms(lb_img)
      krylov_lb <- max(krylov_lb, max(an[zn > 0] / zn[zn > 0]))
    }

    denom <- ip(Rhat, V, ZV)
    bad <- degenerate(denom, Rhat, if (side == "split") ZV else V)
    if (any(bad)) {
      if (drop_broken(bad)) break
      keep <- !bad
      V <- V[, keep, drop = FALSE]; ZPj <- ZPj[, keep, drop = FALSE]
      if (!is.null(ZV)) ZV <- ZV[, keep, drop = FALSE]
      denom <- denom[keep]
    }
    alpha <- rho / denom

    S <- Rr - scale_cols(V, alpha)
    ZS <- if (side == "split") ZR - scale_cols(ZV, alpha) else NULL

    ## The published early exit. A column whose intermediate residual already
    ## meets the target takes the BiCG half of the step and stops there, which is
    ## also what keeps the stabilising half from dividing by a t formed out of
    ## nothing.
    sn <- measure(S, ZS)
    done <- sn <= inner_target
    if (any(done)) {
      if (history) note(active[done], sn[done])
      Xa[, done] <- Xa[, done, drop = FALSE] +
        scale_cols(ZPj[, done, drop = FALSE], alpha[done])
      retire(which(done))
      if (all(done)) { active <- active[FALSE]; break }
      keep <- !done
      shrink(keep)
      V <- V[, keep, drop = FALSE]; ZPj <- ZPj[, keep, drop = FALSE]
      S <- S[, keep, drop = FALSE]; alpha <- alpha[keep]; sn <- sn[keep]
      if (!is.null(ZV)) ZV <- ZV[, keep, drop = FALSE]
      if (!is.null(ZS)) ZS <- ZS[, keep, drop = FALSE]
    }

    ## ----------------------------------------------- the stabilising half --
    ## One step of GMRES(1) on s: omega minimises ||s - omega t||, in whichever
    ## norm this side works in.
    if (side == "left") {
      ZSj <- S
      Tt <- apply_precond(linop_apply(A, S, "N"))
      ZT <- NULL
    } else {
      ZSj <- if (side == "split") ZS else apply_precond(S)
      Tt <- linop_apply(A, ZSj, "N")
      ZT <- if (side == "split") apply_precond(Tt) else NULL
    }

    ts <- ip(Tt, S, ZS)
    bad <- degenerate(ts, Tt, if (side == "split") ZS else S)
    if (any(bad)) {
      if (drop_broken(bad)) break
      keep <- !bad
      Tt <- Tt[, keep, drop = FALSE]; S <- S[, keep, drop = FALSE]
      ZSj <- ZSj[, keep, drop = FALSE]; ZPj <- ZPj[, keep, drop = FALSE]
      V <- V[, keep, drop = FALSE]; alpha <- alpha[keep]
      ts <- ts[keep]
      if (!is.null(ZV)) ZV <- ZV[, keep, drop = FALSE]
      if (!is.null(ZS)) ZS <- ZS[, keep, drop = FALSE]
      if (!is.null(ZT)) ZT <- ZT[, keep, drop = FALSE]
    }
    omega <- ts / ip(Tt, Tt, ZT)

    Xa <- Xa + scale_cols(ZPj, alpha) + scale_cols(ZSj, omega)
    Rr <- S - scale_cols(Tt, omega)
    if (side == "split") ZR <- ZS - scale_cols(ZT, omega)

    ## ------------------------------------------------------- the next step --
    rho_new <- ip(Rhat, Rr, ZR)
    bad <- degenerate(rho_new, Rhat, if (side == "split") ZR else Rr)
    ## A residual that has gone to zero reaches the same test as a shadow vector
    ## that has gone orthogonal to it, and only one of the two is a breakdown.
    rrn <- measure(Rr, ZR)
    bad <- bad & rrn > inner_target
    if (history) note(active, rrn)

    done <- rrn <= inner_target
    if (any(done | bad)) {
      broke <- broke || any(bad)
      retire(which(done | bad))
      if (all(done | bad)) { active <- active[FALSE]; break }
      keep <- !(done | bad)
      shrink(keep)
      V <- V[, keep, drop = FALSE]
      rho_new <- rho_new[keep]; omega <- omega[keep]; alpha <- alpha[keep]
      if (!is.null(ZV)) ZV <- ZV[, keep, drop = FALSE]
    }

    beta <- (rho_new / rho) * (alpha / omega)
    P <- Rr + scale_cols(P - scale_cols(V, omega), beta)
    if (side == "split") ZP <- ZR + scale_cols(ZP - scale_cols(ZV, omega), beta)
    rho <- rho_new
  }

  if (length(active)) retire(seq_along(active))
  list(X = X, iterations = it, krylov_lb = krylov_lb, broke = broke,
       history = hist)
}
