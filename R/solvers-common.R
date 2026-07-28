## What every Krylov method in section 7.1 shares, kept here rather than in the
## first file that needed it. The alternative is a copy per solvers-*.R, and the
## preamble below had reached three identical copies before the fourth method
## made the shape of the duplication obvious.

## The ratio of extreme diagonals of the projected problem's triangular factor at
## which the Krylov space stops being extended.
##
## A Krylov space can go on admitting new directions long after those directions
## have stopped meaning anything: on a fixture with kappa(A) = 3.6e25 GMRES's
## recurrence residual falls monotonically past step 18 while the true residual
## climbs from 0.60 to 1.7e3 and ||x|| reaches 1.7e10, because the minimiser over
## a numerically meaningless space is numerically meaningless and the projected
## problem cannot tell.
##
## max|R_ii| / min|R_ii| is a lower bound on cond(R), so a test against this stops
## later than a sharp estimate would and never earlier, which is the right
## direction for a test that ends an iteration. The same bound reaches the
## bidiagonal methods through a second route: the singular values of the Golub-
## Kahan bidiagonal lie inside those of A, so its condition is a lower bound on
## cond(A) as well, and 1/eps is the limit in both cases.
KRYLOV_CONDITION_LIMIT <- 1 / .Machine$double.eps

## The factor a direction has to lose in the first Gram-Schmidt pass before a
## second one is worth taking, from Daniel, Gragg, Kaufman and Stewart. At
## 1/sqrt(2) a direction that kept more than about 71% of its length is left
## alone, which is the usual setting and is where the measurements in
## dev_notes/gmres-and-the-second-pass.md were taken.
##
## Both families need it and for the same reason. GMRES orthogonalises against a
## stored Krylov basis to build the Hessenberg it minimises over; the two
## eigensolvers orthogonalise against a stored basis to keep a short recurrence
## from losing the orthogonality the projection assumes. One criterion, one place.
REORTH_ETA <- 1 / sqrt(2)

## The preamble: the operator is a linop of the shape the method needs, the
## right-hand side conforms, the budget is a positive integer, and the starting
## iterate has the right shape and the right storage mode.
##
## The iterate lives in the domain and the residual in the codomain, which is one
## space for a square method and two for a least-squares one. Writing both cases
## here rather than in each solver is what keeps a rectangular method from
## carrying its own copy of the checks that have nothing to do with its shape.
solver_setup <- function(A, b, x0, maxit, method, square = TRUE) {
  if (!is_linop(A)) stopf("%s() expects a linop", method)
  m <- A$dim[1L]
  n <- A$dim[2L]
  if (square && m != n) {
    stopf(paste0("%s() needs a square operator; this one is %d x %d.\n",
                 "  A rectangular system is a least-squares problem and takes a different method."),
          method, m, n)
  }

  was_vector <- is.null(dim(b))
  B <- as_block(b)
  if (nrow(B) != m) {
    stopf("non-conformable: operator is %d x %d, right-hand side has %d rows",
          m, n, nrow(B))
  }
  k <- ncol(B)
  maxit <- as.integer(maxit %||% min(10 * as.numeric(n), .Machine$integer.max))
  if (is.na(maxit) || maxit < 1L) stopf("maxit must be a positive integer")

  X <- if (is.null(x0)) {
    matrix(0, n, k)
  } else {
    x <- as_block(x0)
    if (!identical(dim(x), c(n, k))) {
      stopf("x0 is %d x %d; it has to be %d x %d", nrow(x), ncol(x), n, k)
    }
    x
  }
  if (A$dtype == "complex" || is.complex(B) || is.complex(X)) {
    storage.mode(X) <- "complex"
  }

  list(m = m, n = n, k = k, B = B, X = X, maxit = maxit, was_vector = was_vector)
}

## What the two preconditioner sides that transform the residual can find at run
## time and no declaration can rule out in advance.
##
## The split form runs on L^-1 A L^-H for M = L L^H, which exists only for a
## hermitian positive definite M, and the quantity that gives it away is the one
## the method was about to take a square root of. PRECOND_REQUIREMENTS asks these
## methods only for `fixed`, so the check happens where it can be made honestly,
## and the message separates a declaration the iteration has contradicted from a
## property nobody claimed.
split_precond_not_hpd <- function(method, quantity, declared) {
  stopf(paste0("%s() reached %s <= 0 on a split-preconditioned solve.\n",
               "  %s\n",
               "  Split preconditioning runs on L^-1 A L^-H for M = L L^H, which exists only for a\n",
               "  hermitian positive definite M; without one the quantity the method minimises is\n",
               "  not a norm. Use side = 'left' or side = 'right', which ask nothing of M."),
        method, quantity,
        if (declared)
          "The preconditioner declares positive_definite = TRUE and the iteration contradicts it."
        else
          "The preconditioner declares no definiteness, and the split form needs it.")
}

left_precond_singular <- function(method) {
  stopf(paste0("%s() found M^-1 r = 0 for a nonzero residual r.\n",
               "  A left preconditioner has to be invertible: it works on M^-1 A, and a singular M\n",
               "  collapses the residual it starts from to nothing."),
        method)
}
