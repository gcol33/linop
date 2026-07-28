## Section 7.2. What the reference eigensolver and the reference singular value
## solver share, kept here rather than in whichever file needed it first.
##
## Both are labelled reference rather than production, and the label is a
## statement about storage and restarting rather than about accuracy. Neither
## restarts implicitly: a round builds a subspace of at most `ncv` vectors, keeps
## all of them, extracts, locks whatever converged and starts the next round from
## what did not. That is O(ncv) vectors of storage and O(ncv^2 n) work in the
## orthogonalisation, which is the price of not having an implicit restart, and it
## is the gap `linop.primme` closes in Phase 3.
##
## What is not reference-grade about them is the certificate, which is the same
## object every solver in the package returns and rests on the same measurements.

## Section 1.1: which takes strings, and RSpectra's short forms are accepted
## because the incumbent's vocabulary is what users already have in their fingers.
##
## Magnitude and algebraic order are separate entries rather than one flag,
## because for an indefinite hermitian operator they are genuinely different
## requests and a single "smallest" would have to guess which was meant.
EIGEN_WHICH <- c(largest = "LM", smallest = "SM",
                 largest_algebraic = "LA", smallest_algebraic = "SA",
                 LM = "LM", SM = "SM", LA = "LA", SA = "SA")

## Singular values are non-negative, so magnitude order and algebraic order
## coincide and only two of the four names mean anything distinct here.
SVD_WHICH <- c(largest = "LM", smallest = "SM", LM = "LM", SM = "SM")

normalise_which <- function(which, table, verb) {
  if (!is_scalar_string(which)) stopf("which must be a single string")
  if (!which %in% names(table)) {
    stopf("%s does not take which = '%s'; it is one of: %s",
          verb, which, paste(unique(names(table)), collapse = ", "))
  }
  unname(table[[which]])
}

## The order the requested end of the spectrum comes in. order() is stable, so
## equal values keep the order the projected problem produced them in.
ritz_order <- function(theta, which) {
  switch(which,
         LM = order(-abs(theta)),
         SM = order(abs(theta)),
         LA = order(-theta),
         SA = order(theta))
}

## Full reorthogonalisation of the block W against the stored orthonormal basis
## Q, by block classical Gram-Schmidt with a second pass taken only where the
## first one cancelled.
##
## Classical rather than modified: the projection is one matrix product against
## the whole basis instead of a loop over its columns, which is where an R
## implementation of a long recurrence spends its time. Classical Gram-Schmidt on
## its own loses orthogonality at the rate the basis is ill conditioned, and what
## restores it is the second pass rather than the sweep order.
##
## The coefficients are discarded here, and that is what keeps a Lanczos
## recurrence with full reorthogonalisation a Lanczos recurrence: they are a
## correction to a three-term recurrence rather than entries of the projected
## problem. Feeding them into the projection instead would make it Arnoldi. GMRES
## needs them, which is why it orthogonalises in its own loop rather than through
## this one, and the only thing the two share is REORTH_ETA.
##
## The mask is exact, so a column that did not need the second pass has an exact
## zero multiple subtracted from it and is left bit for bit as it was.
orth_against <- function(Q, W) {
  if (is.null(Q) || !ncol(Q)) return(W)
  w0 <- col_norms(W)
  W <- W - Q %*% crossprod(Conj(Q), W)
  again <- col_norms(W) < REORTH_ETA * w0
  if (any(again)) {
    W <- W - Q %*% scale_cols(crossprod(Conj(Q), W), as.numeric(again))
  }
  W
}

## A unit starting vector, in the storage mode the operator implies. Random
## unless the caller supplied one, and seeded either way so a report is
## reproducible without the solver moving the caller's stream.
eigen_start_vector <- function(n, dtype, v0, seed) {
  if (!is.null(v0)) {
    v <- as_block(v0)
    if (nrow(v) != n || ncol(v) != 1L) {
      stopf("v0 is %d x %d; it has to be %d x 1", nrow(v), ncol(v), n)
    }
  } else {
    v <- with_preserved_seed(seed, {
      z <- matrix(stats::rnorm(n), n, 1L)
      if (dtype == "complex") z <- z + 1i * matrix(stats::rnorm(n), n, 1L)
      z
    })
  }
  if (dtype == "complex") storage.mode(v) <- "complex"
  nv <- col_norms(v)
  if (nv == 0) stopf("the starting vector is zero")
  v / nv
}

## A replacement direction, for a round that has to start again from a subspace
## the last one has already exhausted. Random, and seeded on the round so two
## rounds do not draw the same vector.
eigen_random_vector <- function(n, dtype, seed, round) {
  with_preserved_seed(seed + round, {
    z <- matrix(stats::rnorm(n), n, 1L)
    if (dtype == "complex") z <- z + 1i * matrix(stats::rnorm(n), n, 1L)
    if (dtype == "complex") storage.mode(z) <- "complex"
    z
  })
}

## The preamble both verbs share: the operator is a linop of the right shape, the
## number of pairs asked for fits inside it, and the subspace is wide enough to
## hold them and the budget is positive.
eigen_setup <- function(A, k, ncv, maxit, verb, dim_limit) {
  if (!is_linop(A)) stopf("%s() expects a linop", verb)
  k <- as.integer(k)
  if (is.na(k) || k < 1L) stopf("k must be a positive integer")
  if (k > dim_limit) {
    stopf("k = %d exceeds the %d pairs this operator has", k, dim_limit)
  }
  ncv <- as.integer(ncv %||% max(2L * k + 1L, k + 20L))
  if (is.na(ncv) || ncv < k) {
    stopf("ncv = %d cannot hold %d pairs", ncv, k)
  }
  ncv <- min(ncv, dim_limit)
  ## A subspace of exactly k leaves the last pair no direction to converge in,
  ## since the k-1 already locked take k-1 of them. Asking for the whole spectrum
  ## is the one case where that is not a mistake: there one round spans the space.
  if (ncv <= k && k < dim_limit) {
    stopf(paste0("ncv must exceed k; a subspace of exactly %d leaves the last pair no\n",
                 "  direction to converge in once the other %d are locked"), k, k - 1L)
  }
  maxit <- as.integer(maxit %||% max(300L, 20L * k))
  if (is.na(maxit) || maxit < 1L) stopf("maxit must be a positive integer")
  list(k = k, ncv = ncv, maxit = maxit)
}

## --------------------------------------------------------------- the results

new_eigen_result <- function(values, vectors, certificate, k, which, method,
                             iterations, maxit, rounds, converged, residuals,
                             dispatch = NULL) {
  certificate$dispatch <- dispatch
  structure(list(values = values, vectors = vectors, certificate = certificate,
                 k = k, which = which, method = method,
                 iterations = iterations, maxit = maxit,
                 restarts = max(0L, rounds - 1L),
                 nconv = sum(converged), converged = converged,
                 residuals = residuals),
            class = "linop_eigen")
}

new_svd_result <- function(d, u, v, certificate, k, which, method,
                           iterations, maxit, rounds, converged, residuals,
                           dispatch = NULL) {
  certificate$dispatch <- dispatch
  structure(list(d = d, u = u, v = v, certificate = certificate,
                 k = k, which = which, method = method,
                 iterations = iterations, maxit = maxit,
                 restarts = max(0L, rounds - 1L),
                 nconv = sum(converged), converged = converged,
                 residuals = residuals),
            class = "linop_svd")
}

## The header both results print, so the two never drift apart in wording.
print_spectral_header <- function(x, label) {
  cat(sprintf("<%s> %d of %d requested %s converged, which = '%s'\n",
              x$method, x$nconv, x$k, label, x$which))
  cat(sprintf("  %d iterations of at most %d, %d restart%s\n",
              x$iterations, x$maxit, x$restarts,
              if (x$restarts == 1L) "" else "s"))
  invisible(NULL)
}

#' @export
print.linop_eigen <- function(x, ...) {
  print_spectral_header(x, "eigenpairs")
  cat("  values: ", paste(format(x$values, digits = 7), collapse = "  "), "\n", sep = "")
  cat(sprintf("  certificate: %s\n", x$certificate$overall))
  invisible(x)
}

#' @export
print.linop_svd <- function(x, ...) {
  print_spectral_header(x, "singular triplets")
  cat("  values: ", paste(format(x$d, digits = 7), collapse = "  "), "\n", sep = "")
  cat(sprintf("  certificate: %s\n", x$certificate$overall))
  invisible(x)
}

#' @export
summary.linop_eigen <- function(object, ...) {
  print(object)
  print(object$certificate)
  invisible(object)
}

#' @export
summary.linop_svd <- function(object, ...) {
  print(object)
  print(object$certificate)
  invisible(object)
}
