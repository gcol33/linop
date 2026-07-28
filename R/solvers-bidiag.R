## The Golub-Kahan bidiagonalisation, and the loop the two least-squares methods
## of section 7.1 run around it.
##
## From beta_1 u_1 = b - A x_0 and alpha_1 v_1 = A^H u_1,
##
##     beta_{j+1} u_{j+1}  = A v_j       - alpha_j u_j,
##     alpha_{j+1} v_{j+1} = A^H u_{j+1} - beta_{j+1} v_j,
##
## with alpha and beta the norms that normalise each new direction, so that
## A V_k = U_{k+1} B_k for a lower bidiagonal B_k and A^H U_{k+1} = V_k B_k^H +
## alpha_{k+1} v_{k+1} e_{k+1}^H.
##
## LSQR and LSMR build the same B_k. They differ in what they minimise over the
## space it spans: LSQR takes ||b - A x||, which is CG applied to the normal
## equations, and LSMR takes ||A^H (b - A x)||, which is MINRES applied to them.
## Everything before that choice is shared, and that is this file; the choice
## itself is a recurrence each method owns.
##
## Arithmetic. alpha and beta are norms by construction, so B_k is real and
## non-negative whether A is real or complex, and every scalar either recurrence
## forms from it is real. The bidiagonal vectors carry whatever storage mode the
## operator and the right-hand side imply, and nothing else has to.
##
## Two applies per step, one in mode "N" and one in mode "C", so a least-squares
## method cannot run on an operator that supplies only a forward action.

## One run of a recurrence, from the true residual R, until every column meets
## whichever of the two tests applies to it or the budget runs out. The contract
## a recurrence keeps with the loop below:
##
##   in:  A, X, R, active, target, tol, the two preconditioner appliers, the
##        remaining budget, history, conlim, the running lower bound on ||A||
##   out: X with the converged columns updated, the iterations spent, the lower
##        bound as it stands, whether conlim ended a column, and the history rows
##
## Converged columns leave the active set, so the block narrows as the solve
## proceeds and no apply is spent on a column that is already done.
##
## Neither method has a restart parameter: the recurrences are short and store no
## basis, so there is nothing whose growth a round would bound. What ends a round
## is the recurrence claiming convergence on its own estimates, and the loop then
## re-measures. That is the division of labour CG established, for the same
## reason: the estimates drift, and they drift where the answer is most
## converged.
bidiag_solve <- function(A, b, method, recurrence, preconditioner, tol, maxit,
                         x0, history, conlim, floor_const, norm_control) {
  s <- solver_setup(A, b, x0, maxit, method, square = FALSE)
  B <- s$B; X <- s$X; maxit <- s$maxit
  was_vector <- s$was_vector

  check_preconditioner(preconditioner, method)
  apply_precond <- precond_applier(preconditioner)
  apply_precond_adj <- precond_adjoint_applier(preconditioner, method)

  bn <- col_norms(B)
  ## A zero right-hand side has the exact solution 0, where a relative test is
  ## undefined; that column gets an absolute one.
  target <- tol * ifelse(bn > 0, bn, 1)

  iterations <- 0L
  rounds <- 0L
  anorm_lb <- 0
  limited <- FALSE
  hist <- list()
  prev <- NULL

  repeat {
    ## Both quantities the certificate can test are measured here, because which
    ## of them applies is a fact about b and the range of A that the recurrence
    ## is in no position to decide. One apply each, against the two per step the
    ## recurrence spends.
    R <- B - linop_apply(A, X, "N")
    rn <- col_norms(R)
    AtR <- linop_apply(A, R, "C")
    atn <- col_norms(AtR)

    ## Both readings, in the same currency the certificate uses. anorm_lb is a
    ## lower bound on ||A||, which enlarges the optimality ratio, so this test
    ## stops later than the certificate's own and never earlier. Before the first
    ## step it is zero and the test is simply never met, which is correct: an
    ## iteration that has not run has proved nothing.
    met <- rn <= target | (anorm_lb > 0 & atn <= tol * anorm_lb * rn)
    if (all(met)) break
    if (iterations >= maxit) break
    ## Neither quantity improved on the last pass, so the recurrence has nothing
    ## left to recover and looping would only spend the budget. Both are measured
    ## in the currency their own test uses, so the comparison is with the pair
    ## that would have stopped the solve.
    now <- pmin(rn / ifelse(bn > 0, bn, 1),
                ifelse(anorm_lb > 0 & rn > 0, atn / (anorm_lb * rn), Inf))
    if (!is.null(prev) && all(now >= prev)) break
    prev <- now

    step <- recurrence(A, X, R, which(!met), target, tol, apply_precond,
                       apply_precond_adj, maxit - iterations, history,
                       conlim, anorm_lb)
    X <- step$X
    iterations <- iterations + step$iterations
    anorm_lb <- max(anorm_lb, step$anorm_lb)
    limited <- limited || step$limited
    if (history) hist <- c(hist, step$history)
    rounds <- rounds + 1L
    if (step$iterations == 0L) break
  }

  norm_est <- do.call(norm2, c(list(A = A), norm_control))
  cert <- solve_certificate(A, B, X, tol = tol, norm_estimate = norm_est,
                            iterations = iterations, maxit = maxit,
                            floor_const = floor_const,
                            norm_lower_bound = anorm_lb,
                            least_squares = TRUE,
                            stop_reason = if (limited)
                              "the projected bidiagonal reached the condition limit" else NULL)

  list(x = undo_block(X, was_vector),
       certificate = cert,
       method = method,
       iterations = iterations,
       restarts = max(0L, rounds - 1L),
       converged = cert_status(cert, "convergence") == "pass",
       residual = cert$values$residual,
       optimality = cert$values$optimality,
       history = if (history && length(hist)) do.call(rbind, hist) else NULL)
}

## beta_1 u_1 = r_0 and alpha_1 v_1 = (A M^-1)^H u_1 = M^-H A^H u_1, the start
## both recurrences make.
##
## alpha_1 = 0 means A^H r_0 = 0, so x_0 already satisfies the normal equations
## and there is no direction to move in. Such a column is dropped here rather
## than carried as a zero one, which the first division would turn into NaN.
## NULL comes back when nothing is left to run.
bidiag_start <- function(A, R, active, apply_precond_adj) {
  if (!length(active)) return(NULL)
  R <- R[, active, drop = FALSE]
  beta <- col_norms(R)
  if (!any(beta > 0)) return(NULL)
  U <- scale_cols(R, 1 / beta)
  V <- apply_precond_adj(linop_apply(A, U, "C"))
  alpha <- col_norms(V)

  live <- alpha > 0
  if (!any(live)) return(NULL)
  if (!all(live)) {
    active <- active[live]
    U <- U[, live, drop = FALSE]; V <- V[, live, drop = FALSE]
    alpha <- alpha[live]; beta <- beta[live]
  }
  list(active = active, U = U, V = scale_cols(V, 1 / alpha),
       alpha = alpha, beta = beta)
}

## One bidiagonalisation step. Both applies are of A itself, with M^-1 on the way
## in and M^-H on the way out, which is exactly A M^-1 and its adjoint without
## either being formed.
##
## v_{j+1} comes back unnormalised, because the two recurrences reach the point
## where they need it at different places in their own bookkeeping and both need
## its norm before they need it.
##
## ||A z|| / ||z|| <= ||A||_2 for every z, and the u are unit vectors, so
## ||A^H u|| is a second such ratio. Both are lower bounds on ||A|| itself rather
## than on the preconditioned operator, which is what the certificate needs, and
## both are read off applies the step was making anyway.
bidiag_step <- function(A, U, V, alpha, apply_precond, apply_precond_adj,
                        anorm_lb) {
  Z <- apply_precond(V)
  AZ <- linop_apply(A, Z, "N")
  Ut <- AZ - scale_cols(U, alpha)
  beta <- col_norms(Ut)
  U <- scale_cols(Ut, ifelse(beta > 0, 1 / beta, 0))

  AtU <- linop_apply(A, U, "C")
  Vt <- apply_precond_adj(AtU) - scale_cols(V, beta)

  zn <- col_norms(Z)
  if (any(zn > 0)) {
    anorm_lb <- max(anorm_lb, max(col_norms(AZ)[zn > 0] / zn[zn > 0]))
  }
  anorm_lb <- max(anorm_lb, max(col_norms(AtU)))

  list(U = U, Vt = Vt, beta = beta, alpha = col_norms(Vt), anorm_lb = anorm_lb)
}
