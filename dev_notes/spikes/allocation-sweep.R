## Is the Conj copy the last of its kind, or the first one found?
##
## `dev_notes/spikes/eigs-compilable-ceiling.R` bounds what compiling linop into
## C could remove: 1.25x to 1.68x, because 44 to 78 per cent of self time is
## already inside BLAS and FFT. That formula sends everything outside linop's
## namespace to `compiled` and calls it irremovable, and the Conj copy was
## exactly that -- `Conj` is a base primitive, so its allocation was counted as
## already-compiled while it was 18 to 26 per cent of eigs() self time and a
## four-line change in R deleted it.
##
## So the ceiling bounds interpreter overhead in linop's own closures and does
## not bound redundant allocation. This script asks how much redundant
## allocation is left, on the paths the ceiling script did not profile: the six
## solver hot loops rather than the eigensolver's orthogonalisation.
##
## Two parts, because they answer different questions and only one of them
## depends on which state R/ is in.
##
##   A  the helpers, both variants defined here, so the measurement is the same
##      before and after any change to R/. Bitwise identity is checked, not
##      assumed: a candidate that is not bitwise identical is reported as such
##      and disqualified, because the package asserts lockstep identities.
##   B  the solvers, profiled through the installed code, so this half is a
##      before/after pair like eigs-orth-crossover.R.
##
## From the repository root: Rscript dev_notes/spikes/allocation-sweep.R
##
## LINOP_SWEEP_TAG=before names the output files, so the pair can sit side by
## side in results/.

if (!file.exists("DESCRIPTION")) stop("run this from the repository root")

devtools::load_all(".", quiet = TRUE)

TAG <- Sys.getenv("LINOP_SWEEP_TAG", "")
suffix <- if (nzchar(TAG)) paste0("-", TAG) else ""

OUT_DIR <- file.path("dev_notes", "spikes", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
CSV_H <- file.path(OUT_DIR, paste0("allocation-sweep-helpers", suffix, ".csv"))
CSV_S <- file.path(OUT_DIR, paste0("allocation-sweep-solvers", suffix, ".csv"))
LOG   <- file.path(OUT_DIR, paste0("allocation-sweep-detail", suffix, ".txt"))

cat("", file = LOG)
say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOG, append = TRUE)
  flush.console()
}

## Median of enough repetitions to clear the clock, rather than one timing.
## Returns seconds per call.
time_call <- function(expr, seconds = 0.4) {
  f <- function() eval(expr, envir = parent.frame(3L))
  t0 <- Sys.time()
  f()
  one <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  reps <- max(3L, min(2000L, ceiling(seconds / max(one, 1e-6))))
  ts <- numeric(reps)
  for (i in seq_len(reps)) {
    t1 <- Sys.time()
    f()
    ts[i] <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  }
  stats::median(ts)
}

## ======================================================================== A
##
## The helpers, current against candidate.
##
## `current` is copied out of R/aaa-utils.R rather than called through the
## namespace, so part A measures the same two expressions whichever state the
## package is in and the before/after pair of part B cannot move it.

cur_scale_cols <- function(X, v) X * rep(v, each = nrow(X))

## The candidate. `rep(v, each = n)` materialises a full n x k double array to
## hold k distinct values, and the multiply then allocates the result: two
## n x k allocations where one is needed. At k = 1 the expansion is an n-vector
## of one repeated value and a scalar multiply is the same arithmetic, so the
## guard costs a length check and removes an allocation the size of the block.
new_scale_cols <- function(X, v) {
  if (length(v) == 1L) X * v else X * rep(v, each = nrow(X))
}

## Per-column assignment for k > 1. One n x k allocation (the result, which is
## needed) and k length-n temporaries, against the current two n x k. Whether
## the interpreted loop pays for itself is the measurement.
loop_scale_cols <- function(X, v) {
  if (length(v) == 1L) return(X * v)
  for (j in seq_along(v)) X[, j] <- X[, j] * v[j]
  X
}

cur_col_norms <- function(X) {
  if (is.complex(X)) sqrt(colSums(Re(X)^2 + Im(X)^2)) else sqrt(colSums(X^2))
}
## Mod(X)^2 forms one intermediate where Re()^2 + Im()^2 forms four. Real
## storage is untouched, so this is the complex branch only.
mod_col_norms <- function(X) {
  if (is.complex(X)) sqrt(colSums(Mod(X)^2)) else sqrt(colSums(X^2))
}
## The BLAS route, listed because it is the obvious one and it has to be
## disqualified on evidence rather than left unmentioned. dsyrk sums in its own
## order, so this is not the same arithmetic.
blas_col_norms <- function(X) {
  if (!is.complex(X) && ncol(X) == 1L) sqrt(crossprod(X)[1L]) else cur_col_norms(X)
}

cur_col_dot <- function(X, Y) {
  if (is.complex(X) || is.complex(Y)) colSums(Re(X) * Re(Y) + Im(X) * Im(Y))
  else colSums(X * Y)
}
blas_col_dot <- function(X, Y) {
  if (!is.complex(X) && !is.complex(Y) && ncol(X) == 1L) crossprod(X, Y)[1L]
  else cur_col_dot(X, Y)
}

## Two allocations, the second only to change the storage mode of the first.
cur_zero_block <- function(X) {
  Z <- matrix(0, nrow(X), ncol(X))
  if (is.complex(X)) storage.mode(Z) <- "complex"
  Z
}
new_zero_block <- function(X) {
  Z <- if (is.complex(X)) complex(nrow(X) * ncol(X)) else numeric(nrow(X) * ncol(X))
  dim(Z) <- c(nrow(X), ncol(X))
  Z
}

## ------------------------------------------------------ bitwise, or not

## identical() on doubles is bit-level for finite values, which is the property
## the package's lockstep tests rest on. A candidate that fails this is not a
## drop-in whatever it measures, and saying so is the point of checking.
bitwise <- function(a, b) identical(a, b)

check_bitwise <- function(label, f_cur, f_new, cases) {
  ok <- TRUE
  for (nm in names(cases)) {
    args <- cases[[nm]]
    a <- do.call(f_cur, args)
    b <- do.call(f_new, args)
    same <- bitwise(a, b)
    if (!same) {
      worst <- max(abs(as.numeric(a) - as.numeric(b)) /
                     pmax(abs(as.numeric(a)), .Machine$double.xmin))
      say("    %-28s %-16s NOT bitwise (relative %.2e)", label, nm, worst)
      ok <- FALSE
    }
  }
  if (ok) say("    %-28s bitwise identical on all %d cases", label, length(cases))
  ok
}

say("======================================================== A. the helpers")
say("")

set.seed(11L)
mk <- function(n, k, cplx = FALSE) {
  X <- matrix(stats::rnorm(n * k), n, k)
  if (cplx) X <- X + 1i * matrix(stats::rnorm(n * k), n, k)
  X
}

## Cases span the shapes the solvers actually produce: k = 1 is a single right
## hand side, k > 1 is a block, and both storage modes occur.
B_small <- mk(1000L, 1L); B_blk <- mk(1000L, 8L); B_cplx <- mk(1000L, 4L, TRUE)
v1 <- 1.7; v8 <- stats::rnorm(8); v4 <- stats::rnorm(4)

say("  scale_cols")
ok_new <- check_bitwise("guard on length(v) == 1", cur_scale_cols, new_scale_cols, list(
  `k=1` = list(B_small, v1), `k=8` = list(B_blk, v8),
  `k=4 complex` = list(B_cplx, v4), `k=8 scalar v` = list(B_blk, v1)))
ok_loop <- check_bitwise("per-column loop", cur_scale_cols, loop_scale_cols, list(
  `k=1` = list(B_small, v1), `k=8` = list(B_blk, v8),
  `k=4 complex` = list(B_cplx, v4)))

say("  col_norms")
ok_mod <- check_bitwise("Mod(X)^2", cur_col_norms, mod_col_norms, list(
  `k=1` = list(B_small), `k=4 complex` = list(B_cplx)))
ok_blas_n <- check_bitwise("crossprod at k=1", cur_col_norms, blas_col_norms, list(
  `k=1` = list(B_small), `k=8` = list(B_blk)))

say("  col_dot")
ok_blas_d <- check_bitwise("crossprod at k=1", cur_col_dot, blas_col_dot, list(
  `k=1` = list(B_small, mk(1000L, 1L)), `k=8` = list(B_blk, mk(1000L, 8L))))

say("  zero_block")
ok_zero <- check_bitwise("typed allocation", cur_zero_block, new_zero_block, list(
  real = list(B_small), complex = list(B_cplx)))

say("")

## ------------------------------------------------------------- and the cost

NS <- c(1e4, 1e5, 1e6)
KS <- c(1L, 4L, 16L)

hrows <- list()
for (n in NS) {
  for (k in KS) {
    X <- mk(n, k)
    v <- stats::rnorm(k)
    Y <- mk(n, k)
    Xc <- mk(n, k, TRUE)

    entries <- list(
      list("scale_cols", "guard", quote(cur_scale_cols(X, v)), quote(new_scale_cols(X, v)), ok_new),
      list("scale_cols", "loop", quote(cur_scale_cols(X, v)), quote(loop_scale_cols(X, v)), ok_loop),
      list("col_norms complex", "Mod", quote(cur_col_norms(Xc)), quote(mod_col_norms(Xc)), ok_mod),
      list("col_norms", "crossprod", quote(cur_col_norms(X)), quote(blas_col_norms(X)), ok_blas_n),
      list("col_dot", "crossprod", quote(cur_col_dot(X, Y)), quote(blas_col_dot(X, Y)), ok_blas_d),
      list("zero_block complex", "typed", quote(cur_zero_block(Xc)), quote(new_zero_block(Xc)), ok_zero))

    for (e in entries) {
      t_cur <- time_call(e[[3L]])
      t_new <- time_call(e[[4L]])
      hrows[[length(hrows) + 1L]] <- data.frame(
        helper = e[[1L]], candidate = e[[2L]], n = n, k = k,
        cur_seconds = t_cur, new_seconds = t_new,
        speedup = t_cur / t_new, bitwise = e[[5L]],
        stringsAsFactors = FALSE)
    }
    utils::write.csv(do.call(rbind, hrows), CSV_H, row.names = FALSE)
    say("  timed n = %.0e k = %2d", n, k)
    gc(FALSE)
  }
}
H <- do.call(rbind, hrows)
utils::write.csv(H, CSV_H, row.names = FALSE)

say("")
say("%-20s %-10s %8s %4s %9s %8s", "helper", "candidate", "n", "k", "speedup", "bitwise")
for (i in seq_len(nrow(H))) {
  say("%-20s %-10s %8.0e %4d %8.2fx %8s",
      H$helper[i], H$candidate[i], H$n[i], H$k[i], H$speedup[i], H$bitwise[i])
}
say("")

## ======================================================================== B
##
## The solvers, through the package as installed.
##
## Matrix free with a closed form, so a profiled run is checked against truth
## rather than assumed, and validated before any timing runs.

## Dirichlet Laplacian, hermitian positive definite, eigenvalues
## 4 sin^2(k pi / (2(n+1))).
make_spd <- function(n) {
  force(n)
  function(X) {
    X <- as.matrix(X)
    Y <- 2.0001 * X
    Y[-n, ] <- Y[-n, , drop = FALSE] - X[-1L, , drop = FALSE]
    Y[-1L, ] <- Y[-1L, , drop = FALSE] - X[-n, , drop = FALSE]
    Y
  }
}
## The same stencil with an upwind first difference added, so it is
## nonsymmetric and neither CG nor MINRES may run on it.
make_nonsym <- function(n) {
  force(n)
  function(X) {
    X <- as.matrix(X)
    Y <- 2.0001 * X
    Y[-n, ] <- Y[-n, , drop = FALSE] - 1.3 * X[-1L, , drop = FALSE]
    Y[-1L, ] <- Y[-1L, , drop = FALSE] - 0.7 * X[-n, , drop = FALSE]
    Y
  }
}

ref_spd <- function(n) {
  M <- matrix(0, n, n); diag(M) <- 2.0001
  i <- seq_len(n - 1L)
  M[cbind(i, i + 1L)] <- -1; M[cbind(i + 1L, i)] <- -1
  M
}
ref_nonsym <- function(n) {
  M <- matrix(0, n, n); diag(M) <- 2.0001
  i <- seq_len(n - 1L)
  M[cbind(i, i + 1L)] <- -1.3; M[cbind(i + 1L, i)] <- -0.7
  M
}

## A rectangular operator for the two least-squares methods, built from a
## stencil on the first n rows so it stays matrix free.
make_rect <- function(m, n) {
  force(m); force(n)
  f <- make_spd(n)
  list(
    apply = function(X) {
      X <- as.matrix(X)
      rbind(f(X), 0.5 * X[seq_len(m - n), , drop = FALSE])
    },
    adjoint = function(Y) {
      Y <- as.matrix(Y)
      Z <- f(Y[seq_len(n), , drop = FALSE])
      Z[seq_len(m - n), ] <- Z[seq_len(m - n), , drop = FALSE] +
        0.5 * Y[n + seq_len(m - n), , drop = FALSE]
      Z
    })
}
ref_rect <- function(m, n) rbind(ref_spd(n), cbind(0.5 * diag(m - n), matrix(0, m - n, n - (m - n))))

## The callbacks apply the matrices their definitions describe. Checked at a
## small size, since no apply here depends on n.
validate_applies <- function(n = 200L) {
  X <- with_preserved_seed(3L, matrix(stats::rnorm(n * 2L), n, 2L))
  e1 <- max(abs(make_spd(n)(X) - ref_spd(n) %*% X)) / max(abs(ref_spd(n) %*% X))
  e2 <- max(abs(make_nonsym(n)(X) - ref_nonsym(n) %*% X)) / max(abs(ref_nonsym(n) %*% X))
  m <- n + n %/% 2L
  rr <- make_rect(m, n); Mr <- ref_rect(m, n)
  e3 <- max(abs(rr$apply(X) - Mr %*% X)) / max(abs(Mr %*% X))
  Yv <- with_preserved_seed(4L, matrix(stats::rnorm(m * 2L), m, 2L))
  e4 <- max(abs(rr$adjoint(Yv) - t(Mr) %*% Yv)) / max(abs(t(Mr) %*% Yv))
  say("  validate: spd %.2e  nonsym %.2e  rect %.2e  rect adjoint %.2e", e1, e2, e3, e4)
  worst <- max(e1, e2, e3, e4)
  if (!is.finite(worst) || worst > 1e-12) stop(sprintf("a callback disagrees with its definition: %.3e", worst))
  ## Symmetry of the one that claims it, and asymmetry of the one that does not,
  ## so the capability each operator declares is earned rather than asserted.
  Ms <- ref_spd(n); Mn <- ref_nonsym(n)
  if (!isTRUE(all.equal(Ms, t(Ms)))) stop("the spd fixture is not symmetric")
  if (isTRUE(all.equal(Mn, t(Mn)))) stop("the nonsym fixture is symmetric")
  invisible(TRUE)
}

## `maxit` fixed and `tol = 0` so every method does the same number of steps and
## the profiles are comparable. Not a cap on a reported quantity: nothing here
## reports an accuracy, and the accuracy that matters was validated above.
STEPS <- 40L

METHODS <- list(
  cg       = list(kind = "spd",    caps = c(hermitian = TRUE, positive_definite = TRUE)),
  minres   = list(kind = "spd",    caps = c(hermitian = TRUE)),
  gmres    = list(kind = "nonsym", caps = character(0)),
  bicgstab = list(kind = "nonsym", caps = character(0)),
  lsqr     = list(kind = "rect",   caps = character(0)),
  lsmr     = list(kind = "rect",   caps = character(0)))

build_op <- function(kind, n) {
  if (kind == "rect") {
    m <- n + n %/% 2L
    rr <- make_rect(m, n)
    return(linop(rr$apply, adjoint = rr$adjoint, dim = c(m, n)))
  }
  f <- if (kind == "spd") make_spd(n) else make_nonsym(n)
  caps <- METHODS[[if (kind == "spd") "cg" else "gmres"]]$caps
  linop(f, dim = c(n, n), properties = caps)
}

LINOP_NAMES <- ls(asNamespace("linop"), all.names = TRUE)

profile_solver <- function(method, n) {
  spec <- METHODS[[method]]
  A <- if (spec$kind == "rect") build_op("rect", n) else {
    f <- if (spec$kind == "spd") make_spd(n) else make_nonsym(n)
    linop(f, dim = c(n, n), properties = spec$caps)
  }
  b <- with_preserved_seed(5L, stats::rnorm(nrow(A)))

  run <- function() solve(A, b, method = method, maxit = STEPS, tol = 0)

  t0 <- Sys.time(); run()
  one <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  reps <- max(1L, min(20L, ceiling(3 / max(one, 1e-3))))

  pf <- tempfile(fileext = ".Rprof")
  Rprof(pf, interval = 0.002, gc.profiling = FALSE)
  for (i in seq_len(reps)) run()
  Rprof(NULL)
  s <- summaryRprof(pf); unlink(pf)

  self <- s$by.self; total <- s$sampling.time
  grab <- function(nm) {
    row <- paste0("\"", nm, "\"")
    if (row %in% rownames(self)) self[row, "self.time"] else 0
  }
  ## `rep` appears on no path in linop except scale_cols(), so its self time is
  ## that expansion and nothing else.
  rep_t <- grab("rep")
  conj_t <- grab("Conj")
  linop_t <- sum(self$self.time[gsub('^"|"$', "", rownames(self)) %in% LINOP_NAMES])

  say("")
  say("--- %-9s n = %.0e   %.2fs sampled over %d reps", method, n, total, reps)
  top <- head(self[order(-self$self.time), , drop = FALSE], 10L)
  for (i in seq_len(nrow(top))) {
    say("    %-34s self %7.3f  %5.1f%%", rownames(top)[i], top$self.time[i],
        100 * top$self.time[i] / total)
  }

  data.frame(method = method, n = n, steps = STEPS, reps = reps,
             seconds_run = one, sampled_seconds = total,
             rep_self = rep_t, rep_fraction = rep_t / total,
             conj_self = conj_t, conj_fraction = conj_t / total,
             linopR_self = linop_t, linopR_fraction = linop_t / total,
             stringsAsFactors = FALSE)
}

say("======================================================== B. the solvers")
say("")
validate_applies()

srows <- list()
for (n in c(1e5, 1e6)) {
  for (method in names(METHODS)) {
    r <- tryCatch(profile_solver(method, n),
                  error = function(e) { say("!!! %s n = %.0e: %s", method, n, conditionMessage(e)); NULL })
    if (is.null(r)) next
    srows[[length(srows) + 1L]] <- r
    utils::write.csv(do.call(rbind, srows), CSV_S, row.names = FALSE)
    say("  -> rep %.1f%%  Conj %.1f%%  linop R %.1f%%  of %.2fs",
        100 * r$rep_fraction, 100 * r$conj_fraction, 100 * r$linopR_fraction, r$sampled_seconds)
    gc(FALSE)
  }
}
S <- do.call(rbind, srows)
utils::write.csv(S, CSV_S, row.names = FALSE)

say("")
say("================================ summary")
say("%-9s %8s %9s %9s %10s", "method", "n", "rep %", "Conj %", "linop R %")
for (i in seq_len(nrow(S))) {
  say("%-9s %8.0e %8.1f%% %8.1f%% %9.1f%%", S$method[i], S$n[i],
      100 * S$rep_fraction[i], 100 * S$conj_fraction[i], 100 * S$linopR_fraction[i])
}
say("")
say("done: %s  %s", CSV_H, CSV_S)
writeLines(format(Sys.time()), file.path(OUT_DIR, paste0("allocation-sweep", suffix, ".done")))
