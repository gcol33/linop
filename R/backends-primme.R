## The PRIMME engine behind method = "primme", in the shape the reference cores
## return, so that eigs() and svds() certify its pairs exactly as they certify
## their own.
##
## What crosses from the engine is a set of vectors and nothing else that is
## reported. The value reported for an eigenpair is the Rayleigh quotient
## x^H A x / (x^H x) measured here on A, not PRIMME's Ritz value: under a shift
## the engine iterates on (A - sigma I)^-1 and its value is a value of that
## operator, and in every case the certificate's exhibited perturbation
## E = -r x^H - x r^H is exact only when the value is the Rayleigh quotient of x.
## For a singular triplet it is the Rayleigh quotient of the augmented vector,
## Re(u^H A v) for unit u and v. The residuals, and from them convergence, are
## measured on A by the certificate; PRIMME's own residual norms are not used.
##
## An eigenvalue run stops on the rule lanczos_core() locks on, the two rows of
## the eigenpair certificate with its floor, through the engine's convergence
## test. Iterating on A, PRIMME's value and residual norm are already the pair's
## on A, and ||A|| is bounded below by the larger of PRIMME's estimate and every
## |value| seen. Under a shift they belong to (A - sigma I)^-1, whose residual
## bounds the one on A only by ||A - sigma I|| times it, so the test applies A
## to the vector and measures there, bounding ||A|| below as the reference does.
## A check PRIMME makes without a vector cannot be measured on A and does not
## converge the pair.
##
## An eigenvalue run starts from the vector the reference would, so the two
## differ in the method and not in where they began.
##
## The block is k, the number of pairs asked for. Grown from one vector, a
## Krylov-type method reaches one vector of each eigenspace in exact arithmetic
## and only rounding supplies the rest. Asked for the k copies of a k-fold
## eigenvalue, over orders 30 to 200, multiplicities 2 to 8 and ten seeds each,
## a block of k returned every copy every time; a block of one returned 6 or 7
## of an 8-fold one at orders 30 and 60, and the reference, which restarts from
## a random vector, 3 to 6. On the n = 1e6 well the block costs 32.6 s and 166
## applies against 15.1 s and 121. dev_notes/spikes/primme-block-probe.R and
## primme-variants-probe.R.
##
## A PRIMME return code other than 0 (converged) or a budget exhausted is an
## engine failure and is raised as an error naming the code. The eigensolver
## reports an exhausted budget as -3; the singular value solver adds -100 or
## -200 to an eigensolver code according to which of its two stages returned
## it, so there it is -103 or -203. An exhausted budget is not a failure: it
## comes back as the pairs reached, and the certificate says which of them did
## not converge.

PRIMME_EIGS_BUDGET <- -3L
PRIMME_SVDS_BUDGET <- c(-103L, -203L)

primme_failure <- function(info, verb) {
  stopf("%s(method = 'primme'): PRIMME returned error code %d", verb, info)
}

## PRIMME refuses a tolerance below machine epsilon, and reads a tolerance of
## zero as a request for its own default of 1e4 eps. Neither is what tol means
## here, where a tolerance the arithmetic cannot reach is allowed and the
## certificate reports it against the floor.
primme_supports <- function(request) {
  tol <- request$tol
  if (!is.null(tol) && tol < .Machine$double.eps) {
    return(sprintf("it refuses tol = %g, below machine epsilon", tol))
  }
  ## With a block of k, PRIMME's restart has to leave room for the k it keeps
  ## and a block of k more, and below that its own parameter check fails: the
  ## smallest basis it runs with was measured at 2k + 1 for k = 1 to 8 at every
  ## target of both solvers, except a single eigenpair at largest or smallest
  ## magnitude, an interior target for PRIMME, where it is 4.
  ## dev_notes/spikes/primme-block-probe.R.
  ncv <- request$ncv
  if (!is.null(ncv)) {
    need <- if (request$k == 1L && isTRUE(request$interior)) 4L else 2L * request$k + 1L
    if (ncv < need) {
      return(sprintf("with a block of k = %d it needs ncv >= %d", request$k, need))
    }
  }
  TRUE
}

## Columns of an exhausted run that the engine never filled are not pairs and
## are not returned. Nothing is padded, and a run that filled none returns none
## with a certificate saying so.
filled_columns <- function(X) which(col_norms(X) > 0)

primme_eigs_run <- function(Op, A, ctl, engine) {
  n <- Op$dim[1L]
  sel <- ctl$sel
  v <- eigen_start_vector(n, Op$dtype, ctl$v0, ctl$seed)
  on_A <- identical(Op, A)
  lb <- 0
  convtest <- function(value, x, rnorm, anorm) {
    if (on_A) {
      lb <<- max(lb, anorm, abs(value))
      return(eigen_pair_bounds(rnorm, value, lb, 1, ctl$tol, ctl$floor_const)$converged)
    }
    if (is.null(x)) return(FALSE)
    meas <- eigen_measure(A, matrix(x, ncol = 1L))
    lb <<- max(lb, meas$anorm_lb)
    eigen_pair_bounds(meas$rr, meas$mu, lb, 1, ctl$tol, ctl$floor_const)$converged
  }
  out <- engine$fn("primme_eigsh")(
    function(X) linop_apply(Op, X, "N"), n, ctl$k, which = sel, tol = ctl$tol,
    maxit = ctl$maxit, ncv = ctl$ncv, block = ctl$k, v0 = v,
    complex = Op$dtype == "complex", seed = ctl$seed, floor = ctl$floor_const,
    convtest = convtest)
  if (!out$info %in% c(0L, PRIMME_EIGS_BUDGET)) primme_failure(out$info, "eigs")

  keep <- filled_columns(out$vectors)
  X <- out$vectors[, keep, drop = FALSE]
  theta <- out$values[keep]

  ## The same two lower bounds on ||A|| the reference reads, the ones the
  ## convergence test collected, and PRIMME's largest |Ritz value| when the
  ## engine iterated on A itself. Under a shift it iterated on the inverse,
  ## whose values bound nothing about A.
  anorm_lb <- if (on_A) max(out$stats$estimateLargestSVal, lb) else lb
  mu <- rr <- numeric(0)
  if (length(keep)) {
    meas <- eigen_measure(A, X)
    mu <- meas$mu
    rr <- meas$rr
    anorm_lb <- max(anorm_lb, meas$anorm_lb)
  }

  ord <- ritz_order(theta, sel)
  list(values = mu[ord], vectors = X[, ord, drop = FALSE],
       residuals = rr[ord], steps = as.integer(out$stats$numMatvecs),
       rounds = as.integer(out$stats$numRestarts) + 1L, anorm_lb = anorm_lb,
       stop_reason = if (out$info == PRIMME_EIGS_BUDGET) "budget exhausted" else NULL,
       engine = list(stats = out$stats, params = out$params, info = out$info))
}

## Each pair of columns of U and V, normalised, measured on A as a triplet: its
## value Re(u^H A v), the residual of the augmented pair [u; v] / sqrt(2), which
## combines ||A v - d u|| and ||A^H u - d v|| in the norm that vector is unit in,
## and the lower bound on ||A|| the applies give, max ||A v|| and |d|.
svd_measure <- function(A, U, V) {
  U <- scale_cols(U, 1 / col_norms(U))
  V <- scale_cols(V, 1 / col_norms(V))
  AV <- linop_apply(A, V, "N")
  AtU <- linop_apply(A, U, "C")
  d <- col_dot(U, AV)
  ru <- col_norms(AV - scale_cols(U, d))
  rv <- col_norms(AtU - scale_cols(V, d))
  list(U = U, V = V, d = d, rr = sqrt(ru^2 + rv^2) / sqrt(2),
       anorm_lb = max(col_norms(AV), abs(d)))
}

## A singular value run stops on the same rule, the two rows of the eigenpair
## certificate of the augmented pair, through the engine's convergence test.
## PRIMME passes both vectors on the augmented problem and one on the normal
## equations of its first stage, where the other side is the operator applied
## to the one given; the triplet is measured on A either way.
primme_svds_run <- function(A, ctl, engine) {
  m <- A$dim[1L]; n <- A$dim[2L]
  sel <- ctl$sel
  lb <- 0
  convtest <- function(value, u, v, rnorm, anorm) {
    if (is.null(u) && is.null(v)) return(FALSE)
    if (is.null(u)) u <- linop_apply(A, matrix(v, ncol = 1L), "N")
    if (is.null(v)) v <- linop_apply(A, matrix(u, ncol = 1L), "C")
    if (col_norms(as.matrix(u)) == 0 || col_norms(as.matrix(v)) == 0) return(FALSE)
    meas <- svd_measure(A, matrix(u, ncol = 1L), matrix(v, ncol = 1L))
    lb <<- max(lb, meas$anorm_lb)
    eigen_pair_bounds(meas$rr, meas$d, lb, 1, ctl$tol, ctl$floor_const)$converged
  }
  out <- engine$fn("primme_svds")(
    function(X, trans) linop_apply(A, X, if (trans) "C" else "N"),
    m, n, ctl$k, which = sel, tol = ctl$tol, maxit = ctl$maxit, ncv = ctl$ncv,
    block = ctl$k, complex = A$dtype == "complex", seed = ctl$seed,
    convtest = convtest)
  if (!out$info %in% c(0L, PRIMME_SVDS_BUDGET)) primme_failure(out$info, "svds")

  keep <- intersect(filled_columns(out$u), filled_columns(out$v))
  U <- out$u[, keep, drop = FALSE]
  V <- out$v[, keep, drop = FALSE]
  d <- rr <- numeric(0)
  anorm_lb <- lb
  if (length(keep)) {
    meas <- svd_measure(A, U, V)
    U <- meas$U; V <- meas$V; d <- meas$d; rr <- meas$rr
    anorm_lb <- max(anorm_lb, meas$anorm_lb)
  }

  ord <- ritz_order(d, sel)
  list(d = d[ord], u = U[, ord, drop = FALSE], v = V[, ord, drop = FALSE],
       residuals = rr[ord], steps = as.integer(out$stats$numMatvecs),
       rounds = as.integer(out$stats_stage1$numRestarts +
                           out$stats_stage2$numRestarts) + 1L,
       anorm_lb = anorm_lb,
       stop_reason = if (out$info %in% PRIMME_SVDS_BUDGET) "budget exhausted" else NULL,
       engine = list(stats = out$stats, stage1 = out$stats_stage1,
                     stage2 = out$stats_stage2, params = out$params,
                     info = out$info))
}
