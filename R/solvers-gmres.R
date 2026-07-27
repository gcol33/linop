## Section 7.1, the third and fourth of the seven Krylov methods. Restarted
## GMRES, and FGMRES, which is the same iteration with the preconditioned basis
## stored rather than recomputed.
##
## GMRES requires nothing of the operator. CG needs positive definiteness and
## MINRES needs hermitian; this one needs a square operator and an apply, which
## makes it the method that answers when the capability set says nothing at all.
## It does not need an adjoint either: every apply below is mode "N", so an
## operator that supplies only a forward action is solvable here and nowhere else
## in the package.
##
## The Arnoldi relation. Modified Gram-Schmidt on the Krylov space produces
##
##     A V_m = V_{m+1} H_bar_m,      H_bar_m of size (m+1) x m, upper Hessenberg,
##
## and the iterate x_m = x_0 + V_m y minimising ||r_0 - A V_m y|| is found from
##
##     min_y || beta e_1 - H_bar_m y ||,
##
## which plane rotations reduce as they do in MINRES. The difference is that a
## Hessenberg matrix has no short recurrence behind it: every basis vector has to
## be stored and orthogonalised against, which is what restarting exists to bound.
##
## Complex arithmetic is where this departs from the two hermitian methods. Their
## Krylov scalars are real because A is hermitian, so aaa-utils.R computes them in
## real arithmetic. The Gram-Schmidt coefficient <v_i, A v_j> of a general complex
## operator is complex, and so is the sine of the rotation that eliminates
## h_{j+1,j}, so this file uses col_cdot() throughout and carries a rotation
## defined for complex data.
##
## Several right-hand sides run in lockstep, as in CG and MINRES. Each column
## carries its own Hessenberg matrix, its own rotations and its own basis, so the
## iterates are exactly those of per-column GMRES at one block apply per step
## instead of k. A column that converges inside a round is back-solved at the step
## it reached and leaves, so the block narrows mid-round rather than only at a
## restart.

## ------------------------------------------------------------ the three sides

## Section 4.3 leaves GMRES unrestricted, and unlike CG that is not a statement
## that the sides agree. They do not: left, right and split preconditioning
## generate different iterates and minimise different quantities, and this is the
## first method in the package for which `side` selects an algorithm rather than
## being carried along.
##
##   right   Arnoldi on A M^-1. Minimises ||b - A x||_2, the quantity the
##           certificate reports, so nothing has to be converted. The update is
##           x = x0 + M^-1 (V_m y), one extra apply per round and no extra basis.
##
##   left    Arnoldi on M^-1 A. Minimises ||M^-1 (b - A x)||_2. The recurrence is
##           measuring a preconditioned residual, so the stopping test needs the
##           conversion MINRES introduced.
##
##   split   Arnoldi on L^-1 A L^-H for M = L L^H. Minimises ||b - A x||_{M^-1},
##           again not the reported quantity. The preconditioner contract carries
##           M^-1 as a single map and never L, so the split form is reached by the
##           change of variable v_j = L^-H q_j that MINRES uses, which turns
##           euclidean orthogonality of the q_j into M-orthogonality of the v_j:
##
##               h_ij = <v_i, A v_j>,   w = A v_j - sum_i h_ij u_i,
##               h_{j+1,j} = sqrt(<w, M^-1 w>),
##               u_{j+1} = w / h_{j+1,j},   v_{j+1} = M^-1 w / h_{j+1,j},
##
##           with u_j = M v_j carried alongside v_j so that M is never applied.
##           Two bases instead of one, and only M^-1 is ever called.
##
## The split form is defined only for a hermitian positive definite M, since
## M = L L^H is what makes ||.||_{M^-1} a norm. The GMRES row of
## PRECOND_REQUIREMENTS asks only for `fixed`, so that is checked where it can be
## checked honestly, at run time on <w, M^-1 w>, and the message distinguishes a
## contradicted declaration from a property nobody claimed.

## The factor a direction has to lose in the first Gram-Schmidt pass before a
## second one is worth taking, from Daniel, Gragg, Kaufman and Stewart. At
## 1/sqrt(2) a direction that kept more than about 71% of its length is left
## alone, which is the usual setting and is where the measurements in
## dev_notes/gmres-and-the-second-pass.md were taken.
GMRES_REORTH_ETA <- 1 / sqrt(2)

## The ratio of extreme diagonals of the triangular factor at which the projected
## least-squares problem is treated as numerically singular and the Krylov space
## stops being extended.
##
## Plan 7.1 names condition estimation among the things a Krylov method has to get
## right, and this is where GMRES needs it. A Krylov space can go on admitting new
## directions long after those directions have stopped meaning anything: on a
## fixture with kappa(A) = 3.6e25 the recurrence residual falls monotonically past
## step 18 while the true residual climbs from 0.60 to 1.7e3 and ||x|| reaches
## 1.7e10, because the minimiser over a numerically meaningless space is
## numerically meaningless and the projected problem cannot tell.
##
## max|R_ii| / min|R_ii| is a lower bound on cond(R), so this stops later than a
## sharp estimate would and never earlier, which is the right direction for a test
## that ends an iteration. 1/eps is the usual limit, and the same one LSQR and
## LSMR will want.
GMRES_CONDITION_LIMIT <- 1 / .Machine$double.eps

#' Solve A x = b by restarted GMRES
#'
#' @param A A square `linop`. No capability is required, and no adjoint.
#' @param b A vector or block of right-hand sides.
#' @param preconditioner A `preconditioner`, or `NULL`. Its `side` selects the
#'   algorithm; see the note at the top of this file.
#' @param tol Relative residual tolerance, `||b - A x|| <= tol * ||b||`, in the
#'   euclidean norm whatever the preconditioner.
#' @param restart Krylov dimension per round, the `m` of GMRES(m). Storage and
#'   the work of orthogonalisation both grow with it; convergence does too.
#' @param maxit Total iteration budget across all rounds. Defaults to `10 * n`.
#' @param x0 Starting iterate, or `NULL` for zero.
#' @param history Record every column's recurrence residual at every iteration.
#'   That quantity is the norm the chosen side minimises, which is `||r||_2` only
#'   for right preconditioning. Off by default.
#' @param reorth Allow a second Gram-Schmidt pass where the first one cancelled.
#'   On by default. Off is plain modified Gram-Schmidt, which loses orthogonality
#'   on an ill-conditioned Krylov basis.
#' @param conlim Stop extending the Krylov space once the projected problem's
#'   estimated condition exceeds this. A space can keep admitting directions long
#'   after they have stopped meaning anything, and the projected residual cannot
#'   tell.
#' @param flexible Allow the preconditioner to change between applications, and
#'   store the preconditioned basis so the update stays valid. This is FGMRES.
#' @param floor_const `c` in the arithmetic floor of the certificate.
#' @param norm_control Arguments for the `||A||` estimate.
#' @return A list with the solution, its certificate, and the iteration record.
#' @noRd
gmres_solve <- function(A, b, preconditioner = NULL, tol = 1e-8, restart = 30L,
                        maxit = NULL, x0 = NULL, history = FALSE, reorth = TRUE,
                        conlim = GMRES_CONDITION_LIMIT, flexible = FALSE,
                        floor_const = SOLVE_FLOOR_CONST, norm_control = list()) {
  method <- if (flexible) "fgmres" else "gmres"
  if (!is_linop(A)) stopf("%s() expects a linop", method)
  n <- A$dim[2L]
  if (A$dim[1L] != n) {
    stopf(paste0("%s() needs a square operator; this one is %d x %d.\n",
                 "  A rectangular system is a least-squares problem and takes a different method."),
          method, A$dim[1L], A$dim[2L])
  }
  check_preconditioner(preconditioner, method)

  was_vector <- is.null(dim(b))
  B <- as_block(b)
  if (nrow(B) != n) {
    stopf("non-conformable: operator is %d x %d, right-hand side has %d rows",
          n, n, nrow(B))
  }
  k <- ncol(B)
  maxit <- as.integer(maxit %||% min(10 * as.numeric(n), .Machine$integer.max))
  if (is.na(maxit) || maxit < 1L) stopf("maxit must be a positive integer")
  restart <- as.integer(restart)
  if (is.na(restart) || restart < 1L) stopf("restart must be a positive integer")
  ## A Krylov space cannot exceed the dimension of the space it lives in, and an
  ## m above n only allocates arrays no step will ever reach.
  restart <- min(restart, n)

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
    ## recover. GMRES(m) stagnating on a fixed Krylov dimension is a documented
    ## outcome rather than an anomaly, so it is reported and not looped on.
    if (!is.null(prev_rn) && all(rn >= prev_rn)) break
    prev_rn <- rn

    step <- gmres_round(A, X, R, rn, target, apply_precond, side, flexible,
                        min(restart, maxit - iterations), history, reorth,
                        conlim, pd_declared)
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
       method = method,
       iterations = iterations,
       restarts = max(0L, rounds - 1L),
       converged = cert_status(cert, "residual") != "fail",
       residual = cert$values$residual,
       history = if (history && length(hist)) do.call(rbind, hist) else NULL)
}

#' Solve A x = b by flexible GMRES
#'
#' The same iteration as [gmres_solve()], allowed to see a preconditioner that
#' changes between applications. An inner solve run to a loose tolerance is the
#' canonical case, and it is the only method in v0.1 that accepts one.
#'
#' @inheritParams gmres_solve
#' @return A list with the solution, its certificate, and the iteration record.
#' @noRd
fgmres_solve <- function(A, b, preconditioner = NULL, tol = 1e-8, restart = 30L,
                         maxit = NULL, x0 = NULL, history = FALSE, reorth = TRUE,
                         conlim = GMRES_CONDITION_LIMIT,
                         floor_const = SOLVE_FLOOR_CONST, norm_control = list()) {
  gmres_solve(A, b, preconditioner = preconditioner, tol = tol, restart = restart,
              maxit = maxit, x0 = x0, history = history, reorth = reorth,
              conlim = conlim, flexible = TRUE, floor_const = floor_const,
              norm_control = norm_control)
}

## One round of Arnoldi from the true residual R, of at most `steps` steps.
## Returns when every active column has met its target or the round is spent; a
## column that meets it earlier is back-solved at the step it reached and dropped,
## so no apply is spent on a column that is already done.
##
## `rn` is the euclidean norm of R, measured by the caller. It is needed here only
## to convert the caller's target into the currency this side works in, which for
## right preconditioning is the identity.
gmres_round <- function(A, X, R, rn, target, apply_precond, side, flexible,
                        steps, history, reorth, conlim, pd_declared) {
  eps <- .Machine$double.eps
  nc <- ncol(X)
  active <- which(rn > target)
  Xa <- X[, active, drop = FALSE]
  Ra <- R[, active, drop = FALSE]
  ka <- length(active)

  ## The first basis vector, and with it the norm this side minimises.
  U <- NULL
  if (side == "split") {
    Zt <- apply_precond(Ra)
    q <- col_dot(Ra, Zt)
    if (any(q <= 0)) split_precond_not_hpd(pd_declared)
    beta <- sqrt(q)
    U <- list(scale_cols(Ra, 1 / beta))
    V <- list(scale_cols(Zt, 1 / beta))
  } else if (side == "left") {
    W0 <- apply_precond(Ra)
    beta <- col_norms(W0)
    if (any(beta <= 0)) left_precond_singular()
    V <- list(scale_cols(W0, 1 / beta))
  } else {
    beta <- rn[active]
    V <- list(scale_cols(Ra, 1 / beta))
  }
  Zst <- list()

  ## The recurrence minimises its own norm and the caller asked in the euclidean
  ## one. The ratio is exact for the residual this round starts from, and the
  ## outer loop re-measures afterwards, which is what makes an approximate
  ## conversion safe rather than what makes it right.
  inner_target <- beta * (target[active] / rn[active])

  cplx <- is.complex(Ra) || is.complex(V[[1L]])
  zed <- if (cplx) complex(1) else 0
  Hr <- array(zed, c(steps, steps, ka))     # the rotated Hessenberg, per column
  gv <- matrix(zed, steps + 1L, ka)         # the rotated right-hand side
  csm <- matrix(0, steps, ka)               # rotation cosines, always real
  snm <- matrix(zed, steps, ka)             # rotation sines, complex in general
  gv[1L, ] <- beta

  krylov_lb <- 0
  hist <- list()
  it <- 0L
  ## Extreme diagonals of the triangular factor, per column, for the condition
  ## estimate. Seeded so the first step cannot trip it.
  dmax <- rep(0, ka)
  dmin <- rep(Inf, ka)

  ## A column leaves as soon as its own least-squares residual meets the target.
  ## The update is formed from the basis this side builds its iterate out of:
  ## the preconditioned one when it was stored, the plain one otherwise, and for
  ## fixed right preconditioning by applying M^-1 once to the assembled step.
  retire <- function(idx, j) {
    basis <- if (flexible) Zst else V
    upd <- matrix(zed, nrow(Xa), length(idx))
    for (t in seq_along(idx)) {
      c_t <- idx[t]
      y <- back_substitute(matrix(Hr[seq_len(j), seq_len(j), c_t], j, j),
                           gv[seq_len(j), c_t])
      for (i in seq_len(j)) upd[, t] <- upd[, t] + basis[[i]][, c_t] * y[i]
    }
    if (side == "right" && !flexible) upd <- apply_precond(upd)
    X[, active[idx]] <<- Xa[, idx, drop = FALSE] + upd
  }

  ## Everything carried per column narrows here, and nothing narrows anywhere
  ## else: a vector left out of this list keeps the width the round started with
  ## and is silently recycled against the ones that did narrow.
  shrink <- function(keep) {
    active <<- active[keep]
    Xa <<- Xa[, keep, drop = FALSE]
    V <<- lapply(V, function(M) M[, keep, drop = FALSE])
    if (length(Zst)) Zst <<- lapply(Zst, function(M) M[, keep, drop = FALSE])
    if (!is.null(U)) U <<- lapply(U, function(M) M[, keep, drop = FALSE])
    Hr <<- Hr[, , keep, drop = FALSE]
    gv <<- gv[, keep, drop = FALSE]
    csm <<- csm[, keep, drop = FALSE]
    snm <<- snm[, keep, drop = FALSE]
    inner_target <<- inner_target[keep]
    dmax <<- dmax[keep]
    dmin <<- dmin[keep]
  }

  while (it < steps && length(active)) {
    it <- it + 1L
    j <- it

    ## ---------------------------------------------- the new Krylov direction --
    if (side == "right") {
      Zj <- apply_precond(V[[j]])
      if (flexible) Zst[[j]] <- Zj
      W <- linop_apply(A, Zj, "N")
      src <- Zj
    } else {
      W <- linop_apply(A, V[[j]], "N")
      src <- V[[j]]
    }
    ## ||A v|| / ||v|| <= ||A||_2 for every v, so the iteration hands the
    ## certificate a lower bound on the norm at no extra apply. Measured on the
    ## raw operator, before any preconditioner touches the result.
    sn_v <- col_norms(src)
    an <- col_norms(W)
    if (any(sn_v > 0)) krylov_lb <- max(krylov_lb, max(an[sn_v > 0] / sn_v[sn_v > 0]))
    if (side == "left") W <- apply_precond(W)
    ## The size the new direction had before orthogonalisation, which is the scale
    ## every test below is relative to.
    w0n <- col_norms(W)

    ## ------------------------------------------------------ orthogonalisation --
    ## Modified Gram-Schmidt against every previous basis vector, with a second
    ## pass taken only where the first one cancelled: the Daniel-Gragg-Kaufman-
    ## Stewart criterion, a second pass when the direction lost more than a factor
    ## eta of its length.
    ##
    ## Taking it unconditionally is worse than not taking it at all. Once the
    ## first pass has cancelled a direction down to rounding level, the
    ## coefficients the second pass computes are noise, and they are added to the
    ## Hessenberg column rather than discarded; on a Krylov basis that is already
    ## nearly dependent that noise is what the projected problem gets solved with.
    ##
    ## The criterion is per column and the mask below is exact, so a column that
    ## does not need the second pass is left bitwise untouched by it and a column
    ## that does gets what its own solve would have given. That is what keeps the
    ## test on this decision from costing the lockstep identity.
    ##
    ## The coefficient <v_i, w> is complex exactly when the basis or the new
    ## direction is, so the column that accumulates it is allocated from those two
    ## rather than from the storage mode the round started in.
    hcol <- matrix(if (is.complex(W) || is.complex(V[[1L]])) complex(1) else 0,
                   j, ncol(W))
    proj <- if (side == "split") U else V
    for (i in seq_len(j)) {
      ci <- col_cdot(V[[i]], W)
      hcol[i, ] <- hcol[i, ] + ci
      W <- W - scale_cols(proj[[i]], ci)
    }
    if (isTRUE(reorth)) {
      again <- col_norms(W) < GMRES_REORTH_ETA * w0n
      if (any(again)) {
        ## Zero where the first pass sufficed. Subtracting an exact zero multiple
        ## leaves those columns bit for bit as they were.
        mask <- as.numeric(again)
        for (i in seq_len(j)) {
          ci <- col_cdot(V[[i]], W) * mask
          hcol[i, ] <- hcol[i, ] + ci
          W <- W - scale_cols(proj[[i]], ci)
        }
      }
    }

    if (side == "split") {
      Zt <- apply_precond(W)
      q <- col_dot(W, Zt)
      ## A w that has collapsed to rounding level is the Krylov space running
      ## out, which is convergence rather than a contradiction, and is separated
      ## by size exactly as CG separates its benign case. This threshold is the
      ## generous one on purpose: it guards against accusing a correct
      ## preconditioner, so it errs towards silence.
      spent <- col_norms(W) <= sqrt(eps) * w0n
      if (any(q <= 0 & !spent)) split_precond_not_hpd(pd_declared)
      hnext <- sqrt(pmax(q, 0))
    } else {
      hnext <- col_norms(W)
    }
    ## Breakdown is the tight threshold, because what it guards is a division:
    ## after two Gram-Schmidt passes an h_{j+1,j} at rounding level is a genuinely
    ## dependent direction, and normalising by it would put an overflowed vector
    ## into the basis.
    broke <- hnext <= eps * w0n

    ## The Hessenberg column decides the storage mode of everything downstream. A
    ## real operator under a complex preconditioner reaches here with a complex
    ## coefficient and real arrays waiting for it, and assigning one into the
    ## other would discard the imaginary part without saying so.
    if (is.complex(hcol) && !is.complex(Hr)) {
      storage.mode(Hr) <- "complex"
      storage.mode(gv) <- "complex"
      storage.mode(snm) <- "complex"
      storage.mode(X) <- "complex"
      storage.mode(Xa) <- "complex"
      zed <- complex(1)
    }

    ## ------------------------------------------------------------- rotations --
    ## The previous rotations act on this column of the Hessenberg, then a new one
    ## eliminates h_{j+1,j} against the diagonal.
    if (j >= 2L) {
      for (i in seq_len(j - 1L)) {
        tmp <- csm[i, ] * hcol[i, ] + snm[i, ] * hcol[i + 1L, ]
        hcol[i + 1L, ] <- -Conj(snm[i, ]) * hcol[i, ] + csm[i, ] * hcol[i + 1L, ]
        hcol[i, ] <- tmp
      }
    }
    rot <- givens_rotation(hcol[j, ], hnext)

    ## The new diagonal of the triangular factor, and with it the condition of the
    ## projected problem this step would leave behind. A column whose projected
    ## problem has gone singular is retired at the last step that was still worth
    ## solving, before this one is written, so what it keeps is an iterate rather
    ## than the amplification of one.
    nd <- Mod(rot$r)
    ill <- j > 1L & (pmin(dmin, nd) <= 0 |
                     pmax(dmax, nd) > conlim * pmin(dmin, nd))
    if (any(ill)) {
      retire(which(ill), j - 1L)
      if (all(ill)) { active <- active[FALSE]; break }
      keep <- !ill
      shrink(keep)
      hcol <- hcol[, keep, drop = FALSE]
      hnext <- hnext[keep]; broke <- broke[keep]
      W <- W[, keep, drop = FALSE]
      if (side == "split") Zt <- Zt[, keep, drop = FALSE]
      rot <- lapply(rot, function(v) v[keep])
      nd <- nd[keep]
    }
    dmax <- pmax(dmax, nd)
    dmin <- pmin(dmin, nd)

    csm[j, ] <- rot$cs
    snm[j, ] <- rot$sn
    hcol[j, ] <- rot$r
    Hr[seq_len(j), j, ] <- hcol
    gv[j + 1L, ] <- -Conj(rot$sn) * gv[j, ]
    gv[j, ] <- rot$cs * gv[j, ]

    if (history) {
      row <- rep(NA_real_, nc)
      row[active] <- Mod(gv[j + 1L, ])
      hist[[length(hist) + 1L]] <- row
    }

    ## ------------------------------------------------- extend, or let a column go --
    ## h_{j+1,j} = 0 is the happy breakdown: the Krylov space is invariant and the
    ## iterate at this step is the exact solution of the projected problem. The
    ## rotation has already driven gv[j+1] to zero there, so it retires through
    ## the same test rather than through a special case, and the zero direction
    ## the basis would carry is never read again.
    done <- Mod(gv[j + 1L, ]) <= inner_target | broke
    if (j >= steps) done <- rep(TRUE, length(done))

    if (any(done)) {
      retire(which(done), j)
      if (all(done)) { active <- active[FALSE]; break }
      shrink(!done)
      broke <- broke[!done]
      hnext <- hnext[!done]
      W <- W[, !done, drop = FALSE]
      if (side == "split") Zt <- Zt[, !done, drop = FALSE]
    }

    inv <- ifelse(broke, 0, 1 / hnext)
    if (side == "split") {
      U[[j + 1L]] <- scale_cols(W, inv)
      V[[j + 1L]] <- scale_cols(Zt, inv)
    } else {
      V[[j + 1L]] <- scale_cols(W, inv)
    }
  }

  list(X = X, iterations = it, krylov_lb = krylov_lb, history = hist)
}

## The plane rotation that eliminates a real non-negative g against a possibly
## complex f, in the LAPACK convention: cs real, sn complex, cs^2 + |sn|^2 = 1 and
##
##     [  cs        sn ] [ f ]   [ r ]
##     [ -conj(sn)  cs ] [ g ] = [ 0 ].
##
## g is a norm here and so never negative, which is what lets cs stay real. The
## real case falls out of the same expressions, with the phase of f reducing to
## its sign, so one rotation serves both storage modes.
givens_rotation <- function(f, g) {
  af <- Mod(f)
  d <- sqrt(af^2 + g^2)
  cs <- rep(1, length(f))
  sn <- f * 0
  r <- f * 0
  live <- d > 0
  if (any(live)) {
    fl <- f[live]; al <- af[live]; dl <- d[live]
    ## The phase of f, and 1 where f has none, so the rotation is defined at
    ## f = 0 as well: there it is the swap that moves g onto the diagonal.
    ph <- fl * 0 + 1
    nz <- al > 0
    ph[nz] <- fl[nz] / al[nz]
    cs[live] <- al / dl
    sn[live] <- ph * (g[live] / dl)
    r[live] <- ph * dl
  }
  list(cs = cs, sn = sn, r = r)
}

## Back substitution on the small triangular system the rotations leave behind.
## base::backsolve() is real-only, and the projected system of a complex operator
## is complex, so it is done here for both.
##
## A zero on the diagonal is the degenerate breakdown where the projected problem
## has no component in that direction; contributing nothing is what the
## least-squares solution does there.
back_substitute <- function(Rm, g) {
  m <- length(g)
  y <- g
  for (i in seq.int(m, 1L)) {
    if (i < m) y[i] <- y[i] - sum(Rm[i, (i + 1L):m] * y[(i + 1L):m])
    y[i] <- if (Rm[i, i] == 0) 0 else y[i] / Rm[i, i]
  }
  y
}

## Two messages from one check. The GMRES row of the section 4.3 table asks only
## for `fixed`, so a split preconditioner reaching here may have declared
## definiteness or may have declared nothing, and those are different facts: one
## is a declaration the iteration has contradicted, the other is a requirement of
## the method the caller never claimed to meet. Reporting them the same way would
## accuse a caller who never made the claim.
split_precond_not_hpd <- function(declared) {
  stopf(paste0("gmres() reached <w, M^-1 w> <= 0 on a split-preconditioned solve.\n",
               "  %s\n",
               "  Split preconditioning runs on L^-1 A L^-H for M = L L^H, which exists only for a\n",
               "  hermitian positive definite M; without one the quantity the method minimises is\n",
               "  not a norm. Use side = 'left' or side = 'right', which ask nothing of M."),
        if (declared)
          "The preconditioner declares positive_definite = TRUE and the iteration contradicts it."
        else
          "The preconditioner declares no definiteness, and the split form needs it.")
}

left_precond_singular <- function() {
  stopf(paste0("gmres() found M^-1 r = 0 for a nonzero residual r.\n",
               "  A left preconditioner has to be invertible: it builds the Krylov space of M^-1 A,\n",
               "  and a singular M collapses the first basis vector to nothing."))
}
