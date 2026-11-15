## RSpectra behind method = "rspectra": the Spectra library's implicitly
## restarted Arnoldi and Lanczos, the method of ARPACK written in C++, in the
## shape the reference cores return, so that eigs() and svds() certify its pairs
## exactly as they certify their own.
##
## Only a real operator reaches it. Handed a complex vector, RSpectra coerces it
## to double, and a complex operator behind a callback comes back with the
## spectrum of Re(A), converged and wrong (6.946 on a hermitian operator of order
## 40); a complex dense matrix is refused at its door, a callback is not. The
## registry's dtype check refuses a complex operator by name before any call. A
## complex shift makes the operator the iteration runs on complex and is refused
## for the same reason.
##
## A shift reaches it as an operator, never as sigma: RSpectra refuses sigma on a
## function, and the iteration runs on eigs()'s own (A - sigma I)^-1 instead, the
## same solve the native engines run on.
##
## What crosses from the engine is vectors. The value reported is the Rayleigh
## quotient measured on A, as for every engine. The engine stops on its own rule,
## ||r|| < tol max(eps^(2/3), |theta|) on the operator it iterates, since it has
## no convergence hook; on A that is the certificate's residual row without the
## floor, so a converged run passes the residual and backward error rows, and
## under a shift the certificate measures on A what the engine accepted on the
## inverse. A run that does not converge returns only the pairs that did, possibly
## none, and the certificate counts what is missing.
##
## The budget is applies, and Spectra's is restarts. A run spends ncv steps on its
## first factorisation and at most ncv - k on each restart, the vectors a restart
## discards, so maxit applies allow (maxit - ncv) / (ncv - k) restarts; a step of
## the singular value solver is two applies, one in each mode. Measured, the bound
## is met with equality for the symmetric and singular value drivers and one step
## short per restart for the general one, which keeps a vector more against a
## split conjugate pair. A budget below the first factorisation is refused, since
## nothing can stop the run inside it. The applies are counted in the callback
## and that count is reported.
##
## auto never picks it: Arnoldi and Lanczos take every real operator it takes,
## and they stop on the certificate's rule.

RSPECTRA_GENERAL_WHICH <- c(LM = "LM", SM = "SM", LA = "LR", SA = "SR",
                            LI = "LI", SI = "SI")

## The symmetric driver is the one to use where the operator the iteration runs
## on is established hermitian, which under a shift it is not, because the shift
## is applied through GMRES with nothing declared.
rspectra_symmetric <- function(A, request) {
  isTRUE(capv(A, "hermitian")) && is.null(request$sigma)
}

rspectra_eigs_supports <- function(A, request) {
  if (is.complex(request$sigma)) {
    return("a complex shift makes the operator complex, and it computes in real arithmetic")
  }
  n <- A$dim[1L]
  if (n < 3L) return("it needs an operator of order at least 3")
  sym <- rspectra_symmetric(A, request)
  kmax <- if (sym) n - 1L else n - 2L
  if (request$k > kmax) {
    return(sprintf("it computes at most %d eigenpairs of an operator of order %d", kmax, n))
  }
  need <- if (sym) request$k + 1L else request$k + 2L
  if (!is.null(request$ncv) && request$ncv < need) {
    return(sprintf("it needs ncv >= %d for k = %d", need, request$k))
  }
  rspectra_budget_refusal(A, request, n, 1L)
}

rspectra_budget_refusal <- function(A, request, dim_limit, per_step) {
  ncv <- min(request$ncv %||% default_ncv(A, request$k), dim_limit)
  if (!is.null(request$maxit) && request$maxit < per_step * ncv) {
    return(sprintf("its first factorisation spends %d applies, more than maxit = %d",
                   per_step * ncv, request$maxit))
  }
  TRUE
}

rspectra_eigs_run <- function(Op, A, ctl, engine) {
  n <- Op$dim[1L]
  sym <- isTRUE(capv(Op, "hermitian"))
  applies <- 0L
  f <- function(x, args) {
    applies <<- applies + 1L
    as.vector(linop_apply(Op, matrix(x, ncol = 1L), "N"))
  }
  v <- as.numeric(eigen_start_vector(n, "double", ctl$v0, ctl$seed))
  opts <- list(ncv = ctl$ncv, tol = ctl$tol, retvec = TRUE, initvec = v,
               maxitr = rspectra_restarts(ctl$maxit, ctl$ncv, ctl$k, 1L))
  which <- if (sym) ctl$sel else RSPECTRA_GENERAL_WHICH[[ctl$sel]]
  out <- rspectra_call(if (sym) engine$fn("eigs_sym") else engine$fn("eigs"),
                       list(f, ctl$k, n = n, which = which, opts = opts))

  theta <- out$values
  X <- out$vectors
  ## A non-real pair comes back in complex storage with the real eigenvalues
  ## beside it, and their vectors carry an imaginary part of exactly zero.
  if (is.complex(theta) && all(Im(theta) == 0) && all(Im(X) == 0)) {
    theta <- Re(theta); X <- Re(X)
  }
  anorm_lb <- 0
  mu <- rr <- numeric(0)
  if (length(theta)) {
    meas <- eigen_measure(A, X, hermitian = isTRUE(capv(A, "hermitian")))
    mu <- meas$mu; rr <- meas$rr
    anorm_lb <- meas$anorm_lb
  }
  ord <- ritz_order(theta, ctl$sel)
  list(values = mu[ord], vectors = X[, ord, drop = FALSE], residuals = rr[ord],
       steps = applies, rounds = as.integer(out$niter),
       anorm_lb = anorm_lb,
       stop_reason = if (out$nconv < ctl$k) "budget exhausted" else NULL,
       engine = list(nconv = out$nconv, niter = out$niter, nops = out$nops,
                     version = engine$info$version))
}

rspectra_svds_supports <- function(A, request) {
  if (isTRUE(request$v0)) return("it takes no starting vector")
  w <- min(A$dim)
  if (w < 3L) return("it needs both dimensions at least 3")
  if (request$k >= w) {
    return(sprintf("it computes at most %d singular triplets of a %s operator",
                   w - 1L, fmt_dim(A$dim)))
  }
  if (!is.null(request$ncv) && request$ncv <= request$k) {
    return(sprintf("it needs ncv > k = %d", request$k))
  }
  rspectra_budget_refusal(A, request, w, 2L)
}

## Spectra's singular values are the symmetric eigenproblem on A^T A or A A^T,
## so every step of it is two applies, one in each mode.
rspectra_svds_run <- function(A, ctl, engine) {
  m <- A$dim[1L]; n <- A$dim[2L]
  applies <- 0L
  mode_fn <- function(mode) function(x, args) {
    applies <<- applies + 1L
    as.vector(linop_apply(A, matrix(x, ncol = 1L), mode))
  }
  opts <- list(ncv = ctl$ncv, tol = ctl$tol,
               maxitr = rspectra_restarts(ctl$maxit, ctl$ncv, ctl$k, 2L))
  out <- rspectra_call(engine$fn("svds"),
                       list(mode_fn("N"), ctl$k, nu = ctl$k, nv = ctl$k, opts = opts,
                            Atrans = mode_fn("C"), dim = c(m, n)))

  d <- rr <- numeric(0)
  U <- matrix(0, m, 0L); V <- matrix(0, n, 0L)
  anorm_lb <- 0
  nconv <- length(out$d)
  if (nconv) {
    meas <- svd_measure(A, out$u, out$v)
    U <- meas$U; V <- meas$V; d <- meas$d; rr <- meas$rr
    anorm_lb <- meas$anorm_lb
  }
  ord <- ritz_order(d, ctl$sel)
  list(d = d[ord], u = U[, ord, drop = FALSE], v = V[, ord, drop = FALSE],
       residuals = rr[ord], steps = applies, rounds = as.integer(out$niter),
       anorm_lb = anorm_lb,
       stop_reason = if (nconv < ctl$k) "budget exhausted" else NULL,
       engine = list(nconv = nconv, niter = out$niter, nops = out$nops,
                     version = engine$info$version))
}

rspectra_restarts <- function(maxit, ncv, k, per_step) {
  max(0L, as.integer((maxit %/% per_step - ncv) %/% max(1L, ncv - k)))
}

## A run that does not converge warns and returns the pairs that did. That is
## the stop_reason the certificate reports, not a warning to pass on.
rspectra_call <- function(fn, args) {
  withCallingHandlers(do.call(fn, args), warning = function(w) {
    if (grepl("converged, less than k", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}
