## The finite section, and the certificate that carries a statement about an
## operator no computation ever touched.
##
## P_n H P_n on the indices |j| <= n, as a node holding the operator it truncates.
## The relationship is structural rather than annotated: print() and explain()
## walk to the child and show an Inf x Inf operator underneath a finite one, and
## the certificate reads the limiting coefficients off it.

section_apply <- function(op, X, mode) {
  a <- op$args$a
  b <- op$args$b
  m <- length(a)
  ## Real and symmetric, so all four modes are the same action. The recycling is
  ## down each column, which is what a block stored column-major wants.
  Y <- a * X
  Y[-m, ] <- Y[-m, ] + b * X[-1, , drop = FALSE]
  Y[-1, ] <- Y[-1, ] + b * X[-m, , drop = FALSE]
  Y
}

section_materialize <- function(op) {
  m <- op$dim[1L]
  M <- diag(op$args$a, nrow = m)
  M[cbind(seq_len(m - 1L), seq.int(2L, m))] <- op$args$b
  M[cbind(seq.int(2L, m), seq_len(m - 1L))] <- op$args$b
  M
}

#' Truncate an operator on a sequence space
#'
#' `P_n A P_n` on the indices `-n` to `n`, as an operator of dimension
#' `2n + 1`. The result is an ordinary `linop`: it has a matrix, [eigs()] runs on
#' it, and the arithmetic composes with anything else.
#'
#' What it returns holds the operator it truncates, so the printed tree shows the
#' infinite operator underneath the finite one. That is also where the certificate
#' comes from: [eigs()] on a finite section certifies a statement about `A`, not
#' about the truncation, because the truncation is not the thing anyone wanted to
#' know about.
#'
#' `n` has to be large enough that the window sits strictly inside it. Outside the
#' window the coefficients are at their limits and the eigenvalue equation is the
#' free recurrence, and every closed form in the certificate reads the tail that
#' way; a truncation that cuts through the perturbation has none of them.
#'
#' @param A An operator on a sequence space, from [linop_jacobi()].
#' @param n Half-width of the window. The result is `2n + 1` by `2n + 1`.
#' @return A `linop`.
#' @examples
#' H <- linop_jacobi(diagonal = c(1.5, -0.5, 2))
#' F <- finite_section(H, n = 12)
#' dim(F)
#' print(F)
#' @export
finite_section <- function(A, n) {
  P <- jacobi_of(A, "finite_section")
  if (!identical(A, P)) {
    stopf(paste0("finite_section() takes the operator on the sequence space, not a section of it.\n",
                 "  Truncating a truncation is an ordinary submatrix and says nothing further\n",
                 "  about the operator underneath; call finite_section() on the operator itself."))
  }
  n <- as.integer(n)
  if (is.na(n) || length(n) != 1L) stopf("n must be a single whole number")
  if (n <= P$args$radius) {
    stopf(paste0("n = %d cuts through the window; it has to exceed %d.\n",
                 "  Outside the window the coefficients are at their limits and the eigenvalue\n",
                 "  equation is the free recurrence, which is what makes the truncation error a\n",
                 "  closed form. A section that keeps part of the perturbation outside itself\n",
                 "  has no such statement behind it."),
          n, P$args$radius)
  }
  j <- seq.int(-n, n)
  m <- 2L * n + 1L
  caps <- list(hermitian = capability(TRUE, ev_construction()),
               symmetric = capability(TRUE, ev_construction()),
               real = capability(TRUE, ev_construction()))
  new_linop("section", c(m, m), "double", do.call(new_caps, caps),
            list(A = P, n = n, a = jacobi_a(P, j), b = jacobi_b(P, j[-m])),
            cost = 5 * m, finite = TRUE)
}

## ------------------------------------------------------- the certificate ----

## Every row of the certificate below measures an absolute distance on the
## spectrum, so the caller's relative tolerance is made absolute the way
## eigen_certificate() makes it, against the size of the operator. Written once
## because the width search targets the same number the certificate grades
## against, and two spellings of it would let the search stop just short of what
## the certificate wants.
section_target <- function(P, tol) tol * jacobi_norm_bound(P)

## The fifth shape. Every earlier one certifies a computation against the operator
## the computation ran on; this one certifies a statement about an operator that
## has no matrix, from a computation on one that does.
##
## The whole of it rests on one exact identity. Let u solve the finite problem on
## |j| <= n and let u~ be u extended by zero. Then (H - q) u~ agrees with
## (H_n - q) u at every |j| <= n, because the couplings the truncation removed
## multiply entries that are zero, and it vanishes at every |j| > n + 1. It is
## nonzero at exactly two sites,
##
##     r_{n+1} = b_n u_n,      r_{-(n+1)} = b_{-(n+1)} u_{-n},
##
## and both coefficients are b_inf, since n exceeds the window. So
##
##     ||(H - q) u~||^2 = ||(H_n - q) u||^2 + b_inf^2 (u_n^2 + u_{-n}^2)
##
## exactly. The two terms are separately reportable and separately actionable: the
## first is what the eigensolve left behind and answers to tol and ncv, the second
## is what the truncation costs and answers to n. Folding them into one number
## would report the size of the error without saying which knob moves it.
##
## H is self-adjoint, so the residual bound applies to the sum, and
##
##     dist(q, sigma(H)) <= sqrt(residual^2 + truncation^2)
##
## which is the forward-error row. It rests on nothing declared: the class is
## self-adjoint by construction and the identity above is arithmetic. That makes
## it the first unconditional deterministic bound in the package, and the only one
## about an object that was never in memory.
##
## @param A A finite section.
## @param values The Rayleigh quotients eigs() measured, one per pair.
## @param vectors The Ritz vectors.
## @noRd
section_certificate <- function(A, values, vectors, tol, iterations, maxit,
                                floor_const = SOLVE_FLOOR_CONST,
                                requested = NULL, stop_reason = NULL) {
  r <- cert_rows()
  eps <- .Machine$double.eps
  P <- A$args$A
  n <- A$args$n
  m <- A$dim[1L]
  X <- as_block(vectors)
  q <- as.numeric(values)
  k <- ncol(X)
  requested <- requested %||% k

  hn <- jacobi_norm_bound(P)
  band <- jacobi_band(P)
  b_inf <- abs(P$args$b_inf)
  flr <- floor_const * eps * hn
  want <- section_target(P, tol)

  ## Recomputed here rather than taken from the iteration, for the reason every
  ## certificate in the package recomputes: what a recurrence believes about
  ## itself drifts, and it drifts where the answer is most converged.
  xn <- col_norms(X)
  sx <- ifelse(xn > 0, xn, 1)
  res <- col_norms(linop_apply(A, X, "N") - scale_cols(X, q)) / sx
  trunc <- b_inf * sqrt(Mod(X[1L, ])^2 + Mod(X[m, ])^2) / sx
  total <- sqrt(res^2 + trunc^2)
  bound <- total + flr
  met <- bound <= want

  ## ------------------------------------------------------------- floor row --
  r$add("arithmetic floor", "pass",
        sprintf("||H|| <= %.4g by row sums; c = %g, eps = %.3g, floor %.3e",
                hn, floor_const, eps, flr),
        source = "construction", guarantee = "deterministic_bound", confidence = 1)

  ## --------------------------------------------------------- residual rows --
  r$add("finite residual",
        cert_level(max(res), want, flr),
        sprintf("worst ||H_n u - q u|| / ||u|| = %.3e against %.3e; %d of %d Lanczos steps",
                max(res), want, iterations, maxit))

  cert_add_orthogonality(r, X, floor_const)

  ## The decay rate is what says whether the truncation term is still measuring
  ## truncation. It falls like rho^n until it reaches the level at which the
  ## finite eigensolve stores the tail at all, and past that it is a plateau: at
  ## v = 1 and n = 80 the analytic value is 1.1e-17 and the measured one 7.7e-17.
  ## The bound still holds there, because the floor covers exactly that gap, and
  ## the row says which regime it is in rather than changing its verdict.
  rho <- decay_rate(P, q)
  worst <- which.max(trunc)
  r$add("truncation bound",
        cert_level(max(trunc), want, flr),
        sprintf("worst |b| sqrt(u_{-n}^2 + u_n^2) / ||u|| = %.3e against %.3e at n = %d; decay rate %s%s",
                max(trunc), want, n,
                if (is.na(rho[worst])) "undefined inside the band" else sprintf("%.6g", rho[worst]),
                if (max(trunc) <= flr)
                  "; the tail is at the arithmetic floor, so this no longer measures truncation"
                else ""),
        source = "theorem", guarantee = "deterministic_bound", confidence = 1)

  ## ------------------------------------------------------------ isolation ---
  ## sigma_ess is known exactly, so this is the one inequality that separates a
  ## discrete eigenvalue from a discretisation of the continuum. A value with
  ## spectrum within `bound` of it has shown a genuine eigenvalue only if that
  ## whole neighbourhood clears the band; inside it, the finite section is
  ## resolving continuous spectrum and there is no eigenvalue to certify.
  gap <- pmax(q - band[2L], band[1L] - q)
  isolated <- gap > bound
  r$add("isolation", if (all(isolated)) "pass" else "fail",
        sprintf("sigma_ess(H) = [%.6g, %.6g] by Weyl; the closest value clears it by %.3e against the error bound %.3e%s",
                band[1L], band[2L], min(gap), max(bound),
                if (all(isolated)) ""
                else sprintf("; %d of %d values have not separated from the continuum",
                             sum(!isolated), k)),
        source = "theorem", guarantee = "deterministic_bound", confidence = 1)

  ## ------------------------------------------------------------ enclosure ---
  ## Temple's inequality against the band edge, which is quadratic in the residual
  ## where the residual bound is linear: at v = 1 and n = 20 the half-width is
  ## 6.3e-9 against 3.9e-5.
  ##
  ## It needs no spectrum of H between the band edge and the eigenvalue, which is
  ## why it applies to one value per side and not to the rest: the value adjacent
  ## to the edge is the only one with nothing that could be in the way, and for
  ## any other the near end of the interval would be a neighbouring eigenvalue no
  ## residual locates. Even for that one the hypothesis is not provable from a
  ## residual, so the source is the theorem and the guarantee is what the theorem
  ## yields with a hypothesis nobody checked. The three evidence fields being
  ## independent is what lets one row say both.
  ##
  ## The far end is free. The Rayleigh quotient of the extended vector against H
  ## is the Rayleigh quotient of u against H_n, exactly, because the couplings the
  ## truncation removed multiply zeros, so the variational principle puts the
  ## extreme eigenvalue of H beyond q with nothing assumed at all -- beyond it up
  ## to the floor, since q is computed and an exact-arithmetic inequality at the
  ## last bit is not one a certificate may state.
  above <- q > band[2L]
  adjacent <- vapply(seq_len(k), function(i) {
    if (!isolated[i]) return(FALSE)
    others <- seq_len(k) != i & isolated & above == above[i]
    if (above[i]) !any(others & q < q[i]) else !any(others & q > q[i])
  }, logical(1))
  if (!any(adjacent)) {
    r$add("enclosure", "not_checked",
          "no value has separated from the essential spectrum, so there is no isolated eigenvalue to bracket",
          source = "theorem", guarantee = "estimate", confidence = NA_real_)
  } else {
    half <- ifelse(adjacent, bound^2 / gap, NA_real_)
    lo <- ifelse(above, q - flr, q - half)
    hi <- ifelse(above, q + half, q + flr)
    r$add("enclosure", cert_level(max(half[adjacent]), want, flr),
          sprintf("[%.12g, %.12g] for the %d of %d values adjacent to the band, half-width %.3e, assuming no spectrum between the band edge and the eigenvalue",
                  min(lo[adjacent]), max(hi[adjacent]), sum(adjacent), k,
                  max(half[adjacent])),
          source = "theorem", guarantee = "estimate", confidence = NA_real_)
  }

  ## A finite section can be solved to the last digit and still leave a truncation
  ## error six orders above the tolerance, and the two send the reader to different
  ## knobs. Where the eigensolve met its target and the statement did not, the
  ## section is what is short, and saying "stopped early" there would name the
  ## wrong one.
  cert_add_convergence(r, met, requested, iterations, maxit,
                       if (!all(met) && all(res <= want))
                         sprintf("the eigensolve met its target and the section did not; the truncation term is %.3e at n = %d",
                                 max(trunc), n)
                       else stop_reason,
                       "the subspace had stopped improving")

  ## -------------------------------------------------------- forward error ---
  r$add("forward error", if (all(met)) "pass" else "fail",
        sprintf(paste0("min_j |q - lambda_j(H)| <= %.3e = sqrt(%.3e^2 + %.3e^2) + %.3e, ",
                       "against %.3e; the eigenvectors of H are not covered"),
                max(bound), max(res), max(trunc), flr, want),
        evidence = evidence("theorem", "deterministic_bound", 1))

  build_certificate(r$collect(), subject = "finite section",
                    values = list(residual = res, truncation = trunc,
                                  bound = bound, gap = gap, decay = rho,
                                  norm = hn, floor = flr, n = n,
                                  band = band, converged = met),
                    evidence = r$collect_evidence())
}

## eigs() hands every certificate builder the same context, so a node type that
## certifies its own results supplies one function and nothing else.
section_certify <- function(ctx) {
  section_certificate(ctx$A, ctx$values, ctx$vectors, tol = ctx$tol,
                      iterations = ctx$iterations, maxit = ctx$maxit,
                      floor_const = ctx$floor_const, requested = ctx$requested,
                      stop_reason = ctx$stop_reason)
}

register_section_node <- function() {
  linop_register_node("section", section_apply, section_materialize,
                      function(op) 5 * op$dim[1L], overwrite = TRUE,
                      certify = section_certify)
}

## --------------------------------------------------- the width, chosen here --

## eigs() on the operator itself. There is no matrix, so there is a truncation,
## and the point of this is that the caller does not have to pick it.
##
## No width is available before the solve. The tail of an eigenvector falls at
## rho = decay_rate(H, lambda), and lambda is what the solve produces; all that
## is known in advance is that lambda lies outside the band and inside ||H||,
## which puts rho anywhere in (0, 1). So the width is measured rather than
## derived.
##
## One measurement is enough, because outside the window the eigenvector is
## exactly geometric rather than asymptotically so. The certificate's truncation
## term
##
##     eta(n) = |b_inf| sqrt(u_{-n}^2 + u_n^2) / ||u||
##
## therefore satisfies eta(n') = eta(n) rho^(n' - n), and the width that brings
## it to a target is
##
##     n' = n + log(target / eta(n)) / log(rho),
##
## a division rather than a search. Where a doubling scheme spends log2 of the
## overshoot in solves, this spends one solve to measure and one to confirm.
## Each pair has its own rho and its own eta, so the width taken is the largest
## of the k predictions.
##
## Three things stop it, and each is reported rather than retried. The tail can
## stop falling with n, which is the plateau: past the level at which the finite
## eigensolve stores the tail at all, eta measures the eigenvector's own error
## and no width improves it. Nothing may separate from the essential spectrum, in
## which case there is no rate to predict with and no eigenvalue to find, and the
## width doubles until n_max rather than pretending to converge. And n_max itself
## is a budget: a rho within 1e-4 of 1 is an eigenvalue that close to the band
## edge, and the honest answer there is the width it would have taken.

## How much free tail the first section leaves beyond the window, which is what
## the first measurement is taken on: eta is read at the section's edge, and a
## section that barely contains the window has no free region for the eigenvector
## to have become geometric in, so the rate it is predicted with is read off a
## vector that is not yet the right one. A width rather than a free tail would
## mean something different for every window. On a 121-site well the two cost
## three widths (62 -> 80 -> 96) against one (101); on a single site they differ
## by one index.
SECTION_N_START <- 40L

## The widest. m = 2n + 1 is 20001 columns, which at the default subspace is a
## basis of 6 MB, so this is a budget on patience rather than on memory: needing
## more means rho > 0.9995, an eigenvalue within about 1e-7 of the band edge.
SECTION_N_MAX <- 10000L

## What a prediction aims at, as a fraction of the tolerance. eta is predicted
## from a vector that is itself only accurate to its own residual, so aiming
## exactly at the acceptance level would miss about half the time and spend a
## whole extra solve doing it.
SECTION_N_AIM <- 1 / 8

section_eigs <- function(A, k, control, args) {
  P <- jacobi_of(A, "eigs")
  if (!is.null(args$v0)) {
    stopf(paste0("eigs() on the operator itself does not take v0.\n",
                 "  The width of the section is what this call chooses, so a starting vector has\n",
                 "  no length to be yet. Truncate with finite_section(A, n) and pass v0 to\n",
                 "  eigs() on that."))
  }
  bad <- setdiff(names(control), c("n_start", "n_max"))
  if (length(bad)) {
    stopf("section takes n_start and n_max; it does not take: %s",
          paste(bad, collapse = ", "))
  }
  ctl <- utils::modifyList(
    list(n_start = SECTION_N_START, n_max = SECTION_N_MAX), control)
  n_start <- as.integer(ctl$n_start)
  n_max <- as.integer(ctl$n_max)
  if (is.na(n_start) || is.na(n_max) || n_start < 1L) {
    stopf("n_start and n_max must be whole numbers, and n_start at least 1")
  }
  n <- P$args$radius + n_start
  if (n_max < n) {
    stopf(paste0("n_max = %d is below the first width %d.\n",
                 "  The window reaches %d, and n_start = %d of free tail beyond it is what the\n",
                 "  first measurement is taken on."),
          n_max, n, P$args$radius, n_start)
  }

  want <- section_target(P, args$tol)
  ## Half the tolerance rather than all of it, because the residual term is
  ## added in quadrature with this one and answers to a different knob. A search
  ## that stopped at exactly the tolerance would leave the truncation binding.
  accept <- want / 2
  widths <- integer(0)
  prev <- Inf

  repeat {
    widths <- c(widths, n)
    fit <- do.call(eigs, c(list(A = finite_section(P, n), k = k), args))
    v <- fit$certificate$values
    eta <- max(v$truncation)
    if (eta <= accept) { why <- "the tail is below half the tolerance"; break }
    if (n >= n_max) { why <- sprintf("n_max = %d", n_max); break }
    if (eta >= prev) { why <- "widening stopped shrinking the tail"; break }
    prev <- eta
    rate <- is.finite(v$decay) & v$decay < 1 & v$truncation > 0
    n <- if (any(rate)) {
      max(n + 1L, as.integer(ceiling(max(
        v$n + log(want * SECTION_N_AIM / v$truncation[rate]) / log(v$decay[rate])))))
    } else {
      2L * n
    }
    n <- min(n, n_max)
  }

  fit$method <- paste(fit$method, "on a finite section")
  fit$n <- widths[length(widths)]
  fit$widths <- widths
  fit$certificate$values$widths <- widths
  ## The search is a dispatcher choosing, so it is recorded where every other
  ## dispatch in the package is recorded and printed the same way.
  fit$certificate$dispatch <- list(
    requested = args$method, chosen = fit$method,
    reason = sprintf("operator on l^2(Z); width %s, stopped because %s",
                     paste(widths, collapse = " -> "), why))
  fit
}
