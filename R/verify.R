## Section 5.10. The adoption mechanism: a third party writes linop.myclass(),
## runs one function, and knows the ecosystem will work.

#' Check that a claim is valid
#'
#' Called on an operator, asks whether it satisfies the contract. Called on a
#' result together with the operator it came from, asks whether the result is
#' valid. Internally these run different checks; publicly they answer one
#' question.
#'
#' @param x A `linop`, or a solver result.
#' @param ... Passed to methods.
#' @return An object of class `linop_certificate`.
#' @export
verify <- function(x, ...) UseMethod("verify")

#' @rdname verify
#' @param tol Relative tolerance.
#' @param n_probe Number of random probe vectors.
#' @param seed Random seed, so a report is reproducible.
#' @param block Block widths to test.
#' @export
verify.linop <- function(x, tol = 1e-8, n_probe = 10L, seed = 1L, block = c(1L, 3L), ...) {
  A <- x
  ## Gate 1: running the operator suite against a solve object is an error, not
  ## a purity failure. A history-dependent solve returns different output for the
  ## same input by design (section 4.2).
  ##
  ## Ahead of the finite-dimension gate, because this method can be reached
  ## directly with an object that is not a linop at all, and "this is a linsolve"
  ## is the accurate thing to say about it.
  if (is_linsolve(A)) {
    stopf(paste0("this is a linsolve, not a linop, so the operator conformance suite does not apply.\n",
                 "  A solve object may be adaptive and history-dependent by design.\n",
                 "  Use the linsolve suite instead: verify() dispatches to it automatically."))
  }
  ## Gate 2: every check below is a probe, and a probe needs a block.
  require_finite_dim(A, "verify")
  set.seed(seed)
  checks <- list()
  add <- function(name, status, detail = "", source = "computation",
                  guarantee = "identity", confidence = 1) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, status = status, source = source, guarantee = guarantee,
      confidence = confidence, detail = detail, stringsAsFactors = FALSE)
  }

  m <- A$dim[1L]; n <- A$dim[2L]
  cplx <- A$dtype == "complex"
  rnd <- function(nr, nc, force_complex = FALSE) {
    z <- matrix(stats::rnorm(nr * nc), nr, nc)
    if (cplx || force_complex) z <- z + 1i * matrix(stats::rnorm(nr * nc), nr, nc)
    z
  }
  has_adjoint <- !(A$node == "fun" && is.null(A$args$adjoint))
  nrm <- function(v) sqrt(sum(Mod(v)^2))

  ## 1. adjoint consistency -- catches the most common third-party bug, a wrong
  ##    or missing conjugation.
  if (has_adjoint) {
    worst <- 0
    for (i in seq_len(n_probe)) {
      xv <- rnd(n, 1); yv <- rnd(m, 1)
      lhs <- sum(Conj(yv) * linop_apply(A, xv, "N"))
      rhs <- sum(Conj(linop_apply(A, yv, "C")) * xv)
      scale <- max(nrm(xv) * nrm(yv), .Machine$double.eps)
      worst <- max(worst, Mod(lhs - rhs) / scale)
    }
    add("adjoint consistency", if (worst <= tol) "pass" else "fail",
        sprintf("max relative gap %.3e", worst))
  } else {
    add("adjoint consistency", "not_checked", "operator declares no adjoint",
        source = "computation", guarantee = "identity", confidence = NA_real_)
  }

  ## 2. transpose consistency, and the identities tying the view nodes together
  if (has_adjoint) {
    ## adjoint(t(A)) and t(adjoint(A)) are both m x n, so they take n rows.
    Xn <- rnd(n, 2)
    d1 <- max(Mod(linop_apply(make_view(make_view(A, "transpose"), "adjoint"), Xn, "N") -
                  linop_apply(make_view(A, "conjugate"), Xn, "N")))
    d2 <- max(Mod(linop_apply(make_view(make_view(A, "adjoint"), "transpose"), Xn, "N") -
                  linop_apply(make_view(A, "conjugate"), Xn, "N")))
    ## adjoint(Conj(A)) and t(A) are both n x m, so they take m rows.
    Xm <- rnd(m, 2)
    d3 <- max(Mod(linop_apply(make_view(make_view(A, "conjugate"), "adjoint"), Xm, "N") -
                  linop_apply(make_view(A, "transpose"), Xm, "N")))
    worst <- max(d1, d2, d3)
    add("view identities", if (worst <= tol * max(1, worst + 1)) "pass" else "fail",
        sprintf("adjoint(t(A)), t(adjoint(A)), adjoint(Conj(A)); max gap %.3e", worst))

    ## for a real operator the two agree in action while staying distinct nodes
    if (!cplx) {
      g <- max(Mod(linop_apply(A, Xm, "T") - linop_apply(A, Xm, "C")))
      tv <- make_view(A, "transpose"); av <- make_view(A, "adjoint")
      ## Where a view node is actually built it must stay distinct. Leaves that
      ## absorb the view (identity, diag) are simplification rules in their own
      ## right (section 5.8) and legitimately collapse to the same node.
      built <- tv$node %in% names(VIEW_MODE) || av$node %in% names(VIEW_MODE)
      distinct <- !built || !identical(tv$node, av$node)
      add("real transpose equals adjoint",
          if (g <= tol && distinct) "pass" else "fail",
          sprintf("action gap %.3e, nodes distinct: %s%s", g, distinct,
                  if (!built) " (view absorbed by the leaf)" else ""))
    }
  } else {
    add("view identities", "not_checked", "operator declares no adjoint",
        confidence = NA_real_)
  }

  ## 3. linearity -- also the test that catches an adaptive solve wrapped as a linop
  x1 <- rnd(n, 1); x2 <- rnd(n, 1)
  a <- if (cplx) 0.7 + 0.3i else 0.7
  b <- if (cplx) -1.2 + 0.5i else -1.2
  lin <- max(Mod(linop_apply(A, a * x1 + b * x2, "N") -
                 (a * linop_apply(A, x1, "N") + b * linop_apply(A, x2, "N"))))
  ref <- max(nrm(linop_apply(A, x1, "N")), .Machine$double.eps)
  add("linearity", if (lin <= tol * ref) "pass" else "fail",
      sprintf("max gap %.3e", lin))

  ## 4. block consistency
  worst <- 0
  for (bw in block) {
    X <- rnd(n, bw)
    Yb <- linop_apply(A, X, "N")
    for (j in seq_len(bw)) {
      Yj <- linop_apply(A, X[, j, drop = FALSE], "N")
      worst <- max(worst, max(Mod(Yb[, j] - Yj[, 1L])))
    }
  }
  add("block consistency", if (worst <= tol) "pass" else "fail",
      sprintf("widths %s, max gap %.3e", paste(block, collapse = ","), worst))

  ## 5. shape and dtype, for real and complex X
  ok <- TRUE; note <- character()
  Xr <- matrix(stats::rnorm(n * 2), n, 2)
  Yr <- linop_apply(A, Xr, "N")
  if (!identical(dim(Yr), c(m, 2L))) { ok <- FALSE; note <- c(note, "real block shape") }
  Xi <- matrix(stats::rnorm(n * 2), n, 2)
  Yz <- linop_apply(A, Xr + 1i * Xi, "N")
  if (!identical(dim(Yz), c(m, 2L))) { ok <- FALSE; note <- c(note, "complex block shape") }
  if (!is.complex(Yz)) { ok <- FALSE; note <- c(note, "complex input did not promote") }
  add("shape and dtype", if (ok) "pass" else "fail",
      if (ok) "real and complex blocks" else paste(note, collapse = "; "))

  ## 6. complex linearity. Section 5.5 requires that a real operator applied to a
  ##    complex block promotes the result rather than dropping the imaginary
  ##    part. linop_apply() upcasts a real result to complex, which satisfies
  ##    check 5 on its own, so an apply that discarded the imaginary part of its
  ##    input passes every check above. A linear operator satisfies
  ##    A(x + iy) = A(x) + i A(y), and that is what fails when it does.
  gapC <- max(Mod(Yz - (Yr + 1i * linop_apply(A, Xi, "N"))))
  refC <- max(nrm(Yz), .Machine$double.eps)
  add("complex linearity", if (gapC <= tol * refC) "pass" else "fail",
      sprintf("A(x + iy) against A(x) + i A(y); max gap %.3e", gapC))

  ## 7. declared capabilities probed. A contradiction is an error. Agreement adds
  ##    a probe record beside the declaration and never upgrades it further.
  probe_rows <- probe_capabilities(A, tol = tol, n_probe = n_probe, rnd = rnd)
  if (nrow(probe_rows)) {
    contradicted <- probe_rows$verdict == "contradicted"
    add("declared capabilities",
        if (any(contradicted)) "fail" else "pass",
        paste(sprintf("%s: %s", probe_rows$capability, probe_rows$verdict), collapse = "; "),
        source = "probe", guarantee = "heuristic", confidence = NA_real_)
  } else {
    add("declared capabilities", "not_checked", "nothing declared to probe",
        source = "probe", guarantee = "heuristic", confidence = NA_real_)
  }

  ## 8. purity
  X <- rnd(n, 2); Xcopy <- X
  Y1 <- linop_apply(A, X, "N"); Y2 <- linop_apply(A, X, "N")
  pure <- identical(Y1, Y2) && identical(X, Xcopy)
  add("purity", if (pure) "pass" else "fail",
      if (pure) "two applies agree, input unmodified"
      else "output differs between applies, or the input was modified")

  ## 9. gemm agreement with the synthesised form
  X <- rnd(n, 2); Y0 <- rnd(m, 2)
  g <- linop_gemm(A, X, Y0, alpha = 2, beta = -3, mode = "N")
  ref2 <- 2 * linop_apply(A, X, "N") - 3 * Y0
  gd <- max(Mod(g - ref2))
  add("gemm agreement", if (gd <= tol) "pass" else "fail", sprintf("max gap %.3e", gd))

  ## 10. materialisation agreement with column-by-column application
  if (prod(as.numeric(A$dim)) <= 1e6) {
    M <- linop_materialize(A)
    E <- diag(1, n); if (cplx) storage.mode(E) <- "complex"
    Mcol <- linop_apply(A, E, "N")
    md <- max(Mod(as.matrix(M) - Mcol))
    add("materialisation agreement", if (md <= tol) "pass" else "fail",
        sprintf("max gap %.3e", md))
  } else {
    add("materialisation agreement", "not_checked", "operator too large to materialise",
        confidence = NA_real_)
  }

  build_certificate(do.call(rbind, checks), subject = "operator", probes = probe_rows)
}

probe_capabilities <- function(A, tol, n_probe, rnd) {
  n <- A$dim[2L]
  rows <- list()
  declared <- Filter(function(nm) !is.na(capv(A, nm)), CAPABILITY_NAMES)
  square <- A$dim[1L] == A$dim[2L]
  has_adjoint <- !(A$node == "fun" && is.null(A$args$adjoint))

  test <- function(nm, gap) {
    claimed <- capv(A, nm)
    consistent <- if (isTRUE(claimed)) gap <= tol else gap > tol
    rows[[length(rows) + 1L]] <<- data.frame(
      capability = nm, claimed = claimed, gap = gap,
      verdict = if (consistent) "consistent" else "contradicted",
      stringsAsFactors = FALSE)
  }

  for (nm in declared) {
    if (nm == "hermitian" && square && has_adjoint) {
      X <- rnd(n, 2)
      test(nm, max(Mod(linop_apply(A, X, "N") - linop_apply(A, X, "C"))))
    } else if (nm == "symmetric" && square && has_adjoint) {
      X <- rnd(n, 2)
      test(nm, max(Mod(linop_apply(A, X, "N") - linop_apply(A, X, "T"))))
    } else if (nm == "positive_definite" && square) {
      worst <- Inf
      for (i in seq_len(n_probe)) {
        v <- rnd(n, 1)
        q <- Re(sum(Conj(v) * linop_apply(A, v, "N")))
        worst <- min(worst, q / max(sum(Mod(v)^2), .Machine$double.eps))
      }
      ## a positive Rayleigh quotient is not evidence of definiteness; this only
      ## detects a contradiction
      rows[[length(rows) + 1L]] <- data.frame(
        capability = nm, claimed = capv(A, nm), gap = worst,
        verdict = if (isTRUE(capv(A, nm)) && worst <= 0) "contradicted" else "consistent",
        stringsAsFactors = FALSE)
    } else if (nm == "diagonal" && square) {
      X <- rnd(n, 1)
      Y <- linop_apply(A, X, "N")
      d <- linop_apply(A, matrix(1, n, 1), "N")
      test(nm, max(Mod(Y - d * X)))
    }
  }
  if (!length(rows)) {
    return(data.frame(capability = character(), claimed = logical(),
                      gap = numeric(), verdict = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}
