## Section 7.1, the first of the seven Krylov methods. Preconditioned conjugate
## gradients.
##
## Sides. The table in 4.3 leaves CG unrestricted, and that is a statement rather
## than a gap: for a hermitian positive definite M the left, split and right
## preconditioned forms generate the same iterates. Left PCG is CG in the M inner
## product, split PCG is CG on L^-1 A L^-H in the euclidean one, right PCG is CG
## on A M^-1 in the M^-1 one, and the three sequences coincide (Saad, Iterative
## Methods for Sparse Linear Systems, section 9.2). One implementation therefore
## serves all three declared sides, and the residual it reports is always the
## true b - A x rather than a preconditioned one, so nothing downstream depends
## on which side was declared.
##
## Several right-hand sides run in lockstep rather than one after another. Each
## column's recurrence is independent of the others, so the iterates are exactly
## those of per-column CG, and one block apply per step replaces k of them. This
## is the two-tier apply of section 5.2 paying for itself. It is not block CG,
## which couples the columns and is a different method with different iterates.

## Plan 7.1: CG requires positive_definite at declared minimum evidence. Naming
## the method is the caller asserting their own declaration, so an explicit
## method = "cg" accepts any evidence behind a TRUE value; method = "auto" is the
## package choosing on the caller's behalf, and applies this.
##
## computation is admitted alongside construction because an exact check on data
## the operator already holds is an identity, not a probe: linop_diag() proves
## definiteness from the signs of d. Probes are excluded by the guarantee field,
## which they never satisfy, rather than by the source list.
CG_PD_REQUIREMENT <- requirement(
  sources = c("construction", "adapter_contract", "theorem", "computation"),
  guarantees = "identity", min_confidence = 1)

#' Solve A x = b by preconditioned conjugate gradients
#'
#' @param A A square hermitian positive definite `linop`.
#' @param b A vector or block of right-hand sides.
#' @param preconditioner A `preconditioner`, or `NULL`.
#' @param tol Relative residual tolerance, `||b - A x|| <= tol * ||b||`.
#' @param maxit Iteration budget. Defaults to `10 * n`.
#' @param x0 Starting iterate, or `NULL` for zero.
#' @param history Record the residual norm of every column at every iteration.
#'   Off by default: the record is one row per iteration and nothing truncates
#'   it.
#' @param floor_const `c` in the arithmetic floor of the certificate.
#' @param norm_control Arguments for the `||A||` estimate.
#' @return A list with the solution, its certificate, and the iteration record.
#' @noRd
cg_solve <- function(A, b, preconditioner = NULL, tol = 1e-8, maxit = NULL,
                     x0 = NULL, history = FALSE,
                     floor_const = SOLVE_FLOOR_CONST, norm_control = list()) {
  if (!is_linop(A)) stopf("cg() expects a linop")
  n <- A$dim[2L]
  if (A$dim[1L] != n) {
    stopf(paste0("cg() needs a square operator; this one is %d x %d.\n",
                 "  A rectangular system is a least-squares problem and takes a different method."),
          A$dim[1L], A$dim[2L])
  }
  require_capability(A, "hermitian", "cg")
  require_capability(A, "positive_definite", "cg")
  check_preconditioner(preconditioner, "cg")

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

  apply_precond <- function(R) {
    if (is.null(preconditioner)) return(R)
    Z <- as_block(preconditioner$apply_inverse(R))
    if (!identical(dim(Z), dim(R))) {
      stopf("the preconditioner returned a %d x %d block for a %d x %d residual",
            nrow(Z), ncol(Z), nrow(R), ncol(R))
    }
    Z
  }

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
    ## The outer loop measures the true residual; the inner loop trusts the
    ## recurrence. The recurrence drifts, and it drifts exactly where the answer
    ## is most converged, so the decision to stop is never taken on it.
    R <- B - linop_apply(A, X, "N")
    rn <- col_norms(R)
    if (all(rn <= target)) break
    if (iterations >= maxit) break
    ## A restart that did not reduce the true residual has nothing left to
    ## recover. Stopping there is honest; looping is not.
    if (!is.null(prev_rn) && all(rn >= prev_rn)) break
    prev_rn <- rn

    step <- cg_recurrence(A, X, R, target, apply_precond,
                          maxit - iterations, history)
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

  ## A residual line of "qualified" counts as converged. It means the tolerance
  ## was met once the arithmetic floor is allowed for, which is the strongest
  ## statement floating point supports; treating it as failure would reinstate
  ## exactly the S0.6 behaviour the floor exists to remove. The certificate keeps
  ## the distinction for a reader who needs it.
  list(x = undo_block(X, was_vector),
       certificate = cert,
       method = "cg",
       iterations = iterations,
       restarts = max(0L, rounds - 1L),
       converged = cert_status(cert, "residual") != "fail",
       residual = cert$values$residual,
       history = if (history && length(hist)) do.call(rbind, hist) else NULL)
}

## One run of the recurrence, from the true residual R, until every column meets
## its target or the budget runs out. Converged columns leave the active set, so
## the block narrows as the solve proceeds and no apply is spent on a column that
## is already done.
cg_recurrence <- function(A, X, R, target, apply_precond, maxit, history) {
  nc <- ncol(X)
  active <- which(col_norms(R) > target)
  Xa <- X[, active, drop = FALSE]
  Ra <- R[, active, drop = FALSE]
  r0 <- col_norms(Ra)

  Z <- apply_precond(Ra)
  rho <- col_dot(Ra, Z)
  if (any(rho <= 0)) precond_not_pd()
  Pd <- Z

  krylov_lb <- 0
  hist <- list()
  it <- 0L

  while (it < maxit && length(active)) {
    it <- it + 1L
    Q <- linop_apply(A, Pd, "N")

    ## ||A p|| / ||p|| <= ||A||_2 for every p, so the iteration hands the
    ## certificate a lower bound on the norm at no extra apply.
    pn <- col_norms(Pd)
    if (any(pn > 0)) krylov_lb <- max(krylov_lb, max(col_norms(Q)[pn > 0] / pn[pn > 0]))

    pq <- col_dot(Pd, Q)
    bad <- pq <= 0
    if (any(bad)) {
      ## p^H A p <= 0 contradicts positive definiteness. Near convergence p can
      ## underflow, so a column whose residual has already fallen by sqrt(eps)
      ## from where this run started is frozen rather than read as evidence
      ## against the declaration. Anywhere else the declaration is wrong, and
      ## saying so is the point of declaring it.
      benign <- bad & (col_norms(Ra) <= sqrt(.Machine$double.eps) * r0)
      if (any(bad & !benign)) operator_not_pd(sum(bad & !benign))
      X[, active[benign]] <- Xa[, benign, drop = FALSE]
      keep <- !benign
      active <- active[keep]
      if (!length(active)) break
      Xa <- Xa[, keep, drop = FALSE]; Ra <- Ra[, keep, drop = FALSE]
      Pd <- Pd[, keep, drop = FALSE]; Q <- Q[, keep, drop = FALSE]
      rho <- rho[keep]; r0 <- r0[keep]; pq <- pq[keep]
    }

    alpha <- rho / pq
    Xa <- Xa + scale_cols(Pd, alpha)
    Ra <- Ra - scale_cols(Q, alpha)
    rn <- col_norms(Ra)

    if (history) {
      row <- rep(NA_real_, nc)
      row[active] <- rn
      hist[[length(hist) + 1L]] <- row
    }

    done <- rn <= target[active]
    if (any(done)) {
      X[, active[done]] <- Xa[, done, drop = FALSE]
      keep <- !done
      active <- active[keep]
      if (!length(active)) break
      Xa <- Xa[, keep, drop = FALSE]; Ra <- Ra[, keep, drop = FALSE]
      Pd <- Pd[, keep, drop = FALSE]; rho <- rho[keep]; r0 <- r0[keep]
    }

    Z <- apply_precond(Ra)
    rho_new <- col_dot(Ra, Z)
    if (any(rho_new <= 0)) precond_not_pd()
    Pd <- Z + scale_cols(Pd, rho_new / rho)
    rho <- rho_new
  }

  if (length(active)) X[, active] <- Xa
  list(X = X, iterations = it, krylov_lb = krylov_lb, history = hist)
}

operator_not_pd <- function(ncols) {
  stopf(paste0("cg() reached p^H A p <= 0 in %d column(s) while the residual was still large.\n",
               "  The operator declares positive_definite = TRUE and the iteration contradicts it.\n",
               "  CG applies only to a positive definite operator."), ncols)
}

precond_not_pd <- function() {
  stopf(paste0("cg() reached <r, M^-1 r> <= 0.\n",
               "  The preconditioner declares positive_definite = TRUE and the iteration contradicts it.\n",
               "  A preconditioner for CG must be hermitian positive definite on any side."))
}
