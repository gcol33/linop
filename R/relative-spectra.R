## For a coordinate section J = P H P, P H^2 P = J^2 + E^* E, where
## E contains the two couplings across the cut. Squaring J alone discards
## these rows and changes the quadratic relative spectrum.
jacobi_gap_methods <- function() {
  list(second_order = jacobi_second_order, dissipative = jacobi_dissipative)
}

gap_number <- function(x, name, lower = 0, whole = FALSE) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
      !is.finite(x) || x < lower || (whole && x != floor(x))) {
    stopf("%s must be a finite %snumber >= %g", name,
          if (whole) "whole " else "real ", lower)
  }
  x
}

## A symmetric eigendecomposition is checked by its reconstruction and Gram
## defects. The polar factor of U is orthogonal, and
## ||U - polar(U)|| <= e / (1 + sqrt(1-e)), e = ||U^T U - I|| < 1.
## Weyl then encloses each ordered eigenvalue of the Bloch matrix.
jacobi_band_bounds <- function(P) {
  intervals <- lapply(c(1, -1), function(phase) {
    H <- jacobi_bloch(P, phase)
    e <- eigen(H, symmetric = TRUE)
    d <- e$values
    U <- e$vectors
    m <- nrow(H)
    hn <- max(rowSums(abs(H)))
    arithmetic <- (32 * m^2 + 64) * .Machine$double.eps * hn
    gram <- max(rowSums(abs(crossprod(U) - diag(m)))) + (32 * m^2 + 64) * .Machine$double.eps
    if (!is.finite(gram) || gram >= 1) stopf("the Bloch basis could not be validated")
    defect <- max(rowSums(abs(H - tcrossprod(scale_cols(U, d), U))))
    polar <- gram / (1 + sqrt(1 - gram))
    err <- defect + max(abs(d)) * polar * (sqrt(1 + gram) + 1) + arithmetic
    cbind(rev(d - err), rev(d + err))
  })
  cbind(pmin(intervals[[1L]][, 1L], intervals[[2L]][, 1L]),
        pmax(intervals[[1L]][, 2L], intervals[[2L]][, 2L]))
}

jacobi_second_order <- function(A, barrier) {
  J <- section_materialize(A)
  m <- nrow(J)
  edge2 <- jacobi_section_edge(A, diag(m))^2
  ## Eliminating y from J x + y = z x and -E^* E x + J y = z y
  ## gives ((J-z)^2 + E^* E) x = 0.
  C <- rbind(cbind(J, diag(m)), cbind(-diag(edge2, m), J))
  e <- eigen(C)
  take <- which(Im(e$values) >= 0)
  z <- e$values[take]
  X <- e$vectors[seq_len(m), take, drop = FALSE]
  X <- scale_cols(X, 1 / col_norms(X))
  ## A conjugate pair can split into two real numerical roots as its imaginary
  ## part reaches rounding. Collinear real representatives describe one vector.
  keep <- rep(TRUE, length(z))
  real <- integer(0)
  for (i in which(Im(z) == 0)) {
    overlap <- Mod(cross_adjoint(X[, real, drop = FALSE], X[, i, drop = FALSE]))
    if (any(overlap >= 1 - 32 * m * .Machine$double.eps)) keep[i] <- FALSE
    else real <- c(real, i)
  }
  z <- z[keep]
  X <- X[, keep, drop = FALSE]
  JX <- linop_apply(A, X, "N")
  QX <- linop_apply(A, JX, "N") - scale_cols(JX, 2 * z) +
        scale_cols(X, z^2) + edge2 * X
  list(roots = z, vectors = X, residual = col_norms(QX),
       scale = (jacobi_norm_bound(A$args$A) + Mod(z))^2,
       raw_radius = abs(Im(z)), separation = rep(TRUE, length(z)),
       ratio = NULL, mass = NULL, barrier_n = NULL)
}

jacobi_dissipative <- function(A, barrier) {
  n <- A$args$n
  P <- A$args$A
  radius <- barrier$n %||% (P$args$radius + (n - P$args$radius) %/% 2L)
  if (radius < P$args$radius || radius >= n) {
    stopf("barrier$n must cover the perturbation radius %d and be below section width %d",
          P$args$radius, n)
  }
  p <- as.numeric(abs(seq.int(-n, n)) <= radius)
  J <- section_materialize(A)
  diag(J) <- diag(J) + 1i * barrier$gamma * p
  e <- eigen(J)
  X <- scale_cols(e$vectors, 1 / col_norms(e$vectors))
  z <- e$values
  residual <- linop_apply(A, X, "N") + 1i * barrier$gamma * p * X - scale_cols(X, z)
  mass <- colSums(Mod(X)^2 * p) / colSums(Mod(X)^2)
  ratio <- Im(z) / barrier$gamma
  list(roots = z, vectors = X, residual = col_norms(residual),
       scale = jacobi_norm_bound(P) + barrier$gamma + Mod(z),
       raw_radius = NULL, separation = ratio > 1 / 2,
       ratio = ratio, mass = mass, barrier_n = radius)
}

## On a computed quadratic vector x, ||(H-Re(z))x||/||x|| equals |Im(z)|
## in exact arithmetic. Measuring it directly avoids cancellation in the
## quadratic coefficients when two roots approach the real axis.
## The same measurement validates a dissipative candidate on H itself.
jacobi_gap_section <- function(A, k, args) {
  method <- args$method
  if (k > A$dim[1L]) stopf("k exceeds the section dimension %d", A$dim[1L])
  extracted <- jacobi_gap_methods()[[method]](A, args$barrier)
  z <- extracted$roots
  in_gap <- Re(z) > args$gap[1L] & Re(z) < args$gap[2L]
  indices <- which(in_gap)
  target <- args$sigma + if (method == "dissipative") 1i * args$barrier$gamma else 0i
  indices <- indices[order(Mod(z[indices] - target))]
  indices <- utils::head(indices, k)
  roots <- z[indices]
  X <- extracted$vectors[, indices, drop = FALSE]
  q <- Re(roots)
  P <- A$args$A
  hn <- jacobi_norm_bound(P)
  m <- A$dim[1L]
  want <- section_target(P, args$tol)
  r <- cert_rows()
  measured <- band_residuals(A, q, X, function(X) jacobi_section_edge(A, X))
  ## This operation count bounds complex tridiagonal products, subtraction,
  ## squared moduli, norm summation and division by the vector norm. eps is
  ## twice the unit roundoff. The additive term covers gradual underflow.
  neps <- (32 * (m + 2) + 64) * .Machine$double.eps
  if (neps >= 1 / 2) stopf("the section exceeds the roundoff validation budget")
  flr <- neps / (1 - neps) * (hn + abs(q)) +
         args$floor_const * .Machine$double.eps * hn +
         sqrt((m + 2) * .Machine$double.xmin)
  bound <- measured$total + flr
  enclosures <- cbind(lower = q - bound, upper = q + bound)
  isolated <- enclosures[, 1L] > args$gap[1L] & enclosures[, 2L] < args$gap[2L]
  separation <- extracted$separation[indices]
  residual <- extracted$residual[indices] / measured$norm
  scale <- extracted$scale[indices]
  validated <- is.finite(bound) & is.finite(residual) & is.finite(scale) & measured$norm > 0
  if (!all(validated)) stopf("the computed spectral candidates could not be validated")
  met <- bound <= want & isolated & separation
  nonempty <- length(q) > 0L
  level <- function(value, target, floor = 0) {
    if (nonempty) cert_level(max(value), target, max(c(0, floor))) else "not_checked"
  }
  r$add("arithmetic floor", if (nonempty) "pass" else "not_checked",
        sprintf("complex band arithmetic and norm bounds at %d sites; largest radius correction %.3e",
                m, max(c(0, flr))), source = "theorem",
        guarantee = "deterministic_bound", confidence = 1)
  r$add("finite residual", level(measured$residual, want, flr),
        sprintf("worst ||(H_n - Re(z))x||/||x|| = %.3e", max(c(0, measured$residual))))
  r$add("truncation bound", level(measured$truncation, want, flr),
        sprintf("two cut rows of H, worst norm %.3e at n = %d",
                max(c(0, measured$truncation)), A$args$n),
        source = "theorem", guarantee = "deterministic_bound", confidence = 1)
  relative_enclosures <- NULL
  if (method == "second_order") {
    ## Q(z) - (Q(z)x)x^*/||x||^2 is singular and differs from Q(z) by
    ## ||Q(z)x||/||x||. Boulton--Strauss Theorem 2 gives radius sqrt(b^2+delta).
    defect <- residual + neps / (1 - neps) * scale +
              sqrt((m + 2) * .Machine$double.xmin)
    radius <- sqrt(Im(roots)^2 + defect) + flr
    relative_enclosures <- cbind(lower = q - radius, upper = q + radius)
    r$add("second-order inclusion", if (nonempty) "pass" else "not_checked",
          sprintf("sqrt(Im(z)^2 + delta) with checked quadratic defect; largest radius %.3e and scaled residual %.3e; each interval intersects sigma(H)",
                  max(c(0, radius)), max(c(0, residual / scale))),
          evidence = evidence("theorem", "deterministic_bound", 1,
                              depends_on = list(cape(P, "hermitian"))))
  } else {
    r$add("perturbed residual", level(residual / scale, args$tol, neps),
          sprintf("worst scaled perturbed residual = %.3e", max(c(0, residual / scale))))
  }
  if (method == "dissipative") {
    r$add("barrier separation",
          if (!nonempty) "not_checked" else if (all(separation)) "pass" else "fail",
          sprintf("Im(z)/gamma = %s; gamma = %.6g, projection |j| <= %d; finite-width classification",
                  paste(format(extracted$ratio[indices], digits = 8), collapse = ", "),
                  args$barrier$gamma, extracted$barrier_n),
          source = "theorem", guarantee = "estimate", confidence = NA_real_)
  }
  r$add("isolation",
        if (!nonempty) "not_checked" else if (all(isolated)) "pass" else "fail",
        sprintf("each whole interval must lie in the validated essential gap (%.12g, %.12g)",
                args$gap[1L], args$gap[2L]),
        source = "theorem", guarantee = "deterministic_bound", confidence = 1)
  r$add("enclosure", level(bound, want, flr),
        if (nonempty) sprintf("[%.12g, %.12g]; each interval intersects sigma(H)%s; largest radius %.3e against %.3e",
                              min(enclosures[, 1L]), max(enclosures[, 2L]),
                              if (all(isolated)) " in discrete eigenvalues" else "", max(bound), want)
        else "no candidate at this width; no assertion that the gap is empty",
        evidence = evidence("theorem", "deterministic_bound", 1,
                            depends_on = list(cape(P, "hermitian"))))
  r$add("target identity", "not_checked",
        "existence in an isolated interval does not count eigenvalues or establish the nearest one",
        guarantee = "estimate", confidence = NA_real_)
  r$add("gap completeness", "not_checked",
        "a finite relative spectrum or barrier calculation cannot prove that the rest of the gap is empty",
        guarantee = "estimate", confidence = NA_real_)
  cert_add_convergence(r, met, k, 1L, 1L, "section accuracy or gap isolation is insufficient",
                       "section accuracy or gap isolation is insufficient")
  v <- list(residual = measured$residual, truncation = measured$truncation,
            bound = bound, target = want, floor = flr, n = A$args$n,
            roots = roots, enclosures = enclosures, gap = args$gap,
            relative_enclosures = relative_enclosures,
            bands = args$bands, converged = met,
            pencil_residual = residual / scale,
            raw_radius = extracted$raw_radius[indices],
            barrier_ratio = extracted$ratio[indices],
            projection_mass = extracted$mass[indices], barrier_n = extracted$barrier_n)
  cert <- build_certificate(r$collect(), "finite section", values = v,
                            evidence = r$collect_evidence())
  fit <- new_eigen_result(q, X, cert, k, "nearest shift in essential gap", method,
                          1L, 1L, 1L, met, measured$total,
                          dispatch = list(requested = method, chosen = method,
                                          reason = "dense coordinate-section extraction"))
  fit$roots <- roots
  fit$enclosures <- enclosures
  fit$all_roots <- extracted$roots
  fit
}

jacobi_gap_hooks <- function() {
  hook <- function(A, k, control, args) {
    do.call(jacobi_gap_eigs, c(list(A = A, k = k, section = control), args))
  }
  m <- names(jacobi_gap_methods())
  stats::setNames(rep(list(hook), length(m)), m)
}

jacobi_gap_eigs <- function(A, k, which, sigma, B, tol, maxit, ncv, v0, seed,
                            method, inner, section, barrier, floor_const, norm_control) {
  if (!is.null(B) || !is.null(maxit) || !is.null(ncv) || !is.null(v0) ||
      length(inner) || length(norm_control)) {
    stopf("method = '%s' uses dense section extraction; B, maxit, ncv, v0, inner and norm_control do not apply",
          method)
  }
  if (!identical(which, "largest")) stopf("a gap method takes sigma; which does not apply")
  gap_number(k, "k", 1, TRUE)
  gap_number(tol, "tol", .Machine$double.xmin)
  gap_number(floor_const, "floor_const")
  if (is.null(sigma)) stopf("a gap method requires a real sigma inside an essential-spectrum gap")
  gap_number(sigma, "sigma", -Inf)
  P <- if (A$node == "section") A$args$A else A
  if (!all(is.finite(c(P$args$a, P$args$b, P$args$a_inf, P$args$b_inf)))) {
    stopf("gap methods require finite Jacobi coefficients")
  }
  bands <- jacobi_band_bounds(P)
  gaps <- if (nrow(bands) > 1L) cbind(utils::head(bands[, 2L], -1L), utils::tail(bands[, 1L], -1L))
          else matrix(numeric(0), 0L, 2L)
  hit <- which(sigma > gaps[, 1L] & sigma < gaps[, 2L])
  if (length(hit) != 1L) stopf("sigma must lie strictly inside a validated gap between essential bands")
  if (!is.list(barrier)) stopf("barrier must be a list")
  if (length(barrier) && (is.null(names(barrier)) || any(!nzchar(names(barrier))))) {
    stopf("barrier arguments must be named")
  }
  if (method == "second_order" && length(barrier)) stopf("barrier applies to method = 'dissipative'")
  bad <- setdiff(names(barrier), c("gamma", "n"))
  if (length(bad)) stopf("barrier takes gamma and n; it does not take: %s", paste(bad, collapse = ", "))
  barrier <- utils::modifyList(list(gamma = 0.01, n = NULL), barrier)
  gap_number(barrier$gamma, "barrier$gamma", .Machine$double.xmin)
  if (!is.null(barrier$n)) gap_number(barrier$n, "barrier$n", 0, TRUE)
  args <- list(method = method, sigma = sigma, tol = tol, floor_const = floor_const,
               barrier = barrier, gap = gaps[hit, ], bands = bands)
  if (A$node == "section") {
    if (length(section)) stopf("section width is already fixed by finite_section()")
    return(jacobi_gap_section(A, k, args))
  }
  ctl <- section_width_control(section, SECTION_N_START, P$args$radius)
  n <- ctl$n
  if (!is.null(barrier$n) && barrier$n >= n) stopf("barrier$n must be below the first section width %d", n)
  width_search(P, n, ctl$n_max, k, args, accept = section_target(P, tol) / 2,
               predict = function(v, history) 2L * utils::tail(history$n, 1L),
               run = jacobi_gap_section,
               measure = function(v) if (length(v$bound)) max(v$bound) else Inf,
               stop = function(v) if (!length(v$bound)) "no candidate at this width; gap emptiness is unproved")
}
