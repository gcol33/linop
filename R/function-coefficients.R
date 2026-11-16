standard_chop <- function(coeffs, tol = .Machine$double.eps) {
  n <- length(coeffs)
  if (tol >= 1) return(1L)
  if (n < 17L) return(n)
  envelope <- rev(cummax(rev(Mod(coeffs))))
  if (envelope[1L] == 0) return(1L)
  envelope <- envelope / envelope[1L]
  for (j in 2:n) {
    j2 <- floor(1.25 * j + 5.5)
    if (j2 > n) return(n)
    e1 <- envelope[j]
    if (e1 == 0 || envelope[j2] / e1 > 3 * (1 - log(e1) / log(tol))) break
  }
  plateau <- j - 1L
  if (envelope[plateau] == 0) return(plateau)
  j3 <- sum(envelope >= tol^(7 / 6))
  if (j3 < j2) {
    j2 <- j3 + 1L
    envelope[j2] <- tol^(7 / 6)
  }
  max(which.min(log10(envelope[seq_len(j2)]) +
                 seq(0, -log10(tol) / 3, length.out = j2)) - 1L, 1L)
}

function_sampler <- function(f, name) {
  scalar <- FALSE
  calls <- 0L
  checked <- function(y, n) {
    if (!(is.numeric(y) || is.complex(y)) || !is.null(dim(y)) || length(y) != n) {
      stopf("%s function must return a numeric vector of length %d", name, n)
    }
    if (any(!is.finite(y))) stopf("%s function returned non-finite values", name)
    y
  }
  one <- function(x) {
    calls <<- calls + 1L
    y <- tryCatch(f(x), error = function(e)
      stopf("%s function failed: %s", name, conditionMessage(e)))
    checked(y, 1L)
  }
  list(eval = function(x) {
    if (scalar) return(vapply(x, one, complex(1)))
    calls <<- calls + 1L
    y <- tryCatch(f(x), error = identity)
    if (inherits(y, "error")) {
      scalar <<- TRUE
      return(vapply(x, one, complex(1)))
    }
    checked(y, length(x))
  }, calls = function() calls, scalar = function() scalar)
}

cheb_transform <- function(y, grid) {
  n <- length(y)
  c <- stats::fft(c(y, rev(y[2:(n - 1L)])))[seq_len(n)] / (n - 1L)
  c[c(1L, n)] <- c[c(1L, n)] / 2
  if (!is.complex(y) || all(Im(y) == 0)) Re(c) else c
}

fourier_transform <- function(y, grid) {
  n <- length(y)
  c <- stats::fft(y) / n
  d <- (n - 1L) %/% 2L
  k <- seq.int(-d, d)
  c[k %% n + 1L]
}

fourier_chop <- function(c, tol) {
  n <- length(c); d <- (n - 1L) %/% 2L
  paired <- c(Mod(c[d + 1L]), rep(Mod(c[d + 1L + seq_len(d)]) +
                                 Mod(c[d + 1L - seq_len(d)]), each = 2L))
  cut <- standard_chop(paired, tol)
  list(cut = 2L * (cut %/% 2L) + 1L, happy = cut < n)
}

basis_recurrence <- function(kind, n) {
  k <- seq_len(n)
  if (kind == "hermite") list(a = sqrt(k / 2), b = numeric(n))
  else list(a = -k, b = 2 * (k - 1L) + 1)
}

## Implicit QL diagonalizes a Jacobi matrix while rotating its first coordinate.
gauss_jacobi <- function(kind, n) {
  rec <- basis_recurrence(kind, n)
  d <- rec$b; e <- rec$a; e[n] <- 0
  z <- c(1, numeric(n - 1L))
  for (l in seq_len(n)) {
    iter <- 0L
    repeat {
      m <- l
      while (m < n && abs(e[m]) > .Machine$double.eps * (abs(d[m]) + abs(d[m + 1L]))) m <- m + 1L
      if (m == l) break
      iter <- iter + 1L
      if (iter > 100L) stopf("%s Gauss transform: Jacobi QL did not converge", kind)
      g <- (d[l + 1L] - d[l]) / (2 * e[l])
      r <- sqrt(g^2 + 1)
      g <- d[m] - d[l] + e[l] / (g + if (g >= 0) r else -r)
      s <- 1; c <- 1; p <- 0
      for (i in seq.int(m - 1L, l)) {
        f <- s * e[i]; b <- c * e[i]
        if (abs(g) <= abs(f)) {
          c <- g / f; r <- sqrt(c^2 + 1)
          e[i + 1L] <- f * r; s <- 1 / r; c <- c * s
        } else {
          s <- f / g; r <- sqrt(s^2 + 1)
          e[i + 1L] <- g * r; c <- 1 / r; s <- s * c
        }
        g <- d[i + 1L] - p
        r <- (d[i] - g) * s + 2 * c * b
        p <- s * r; d[i + 1L] <- g + p; g <- c * r - b
        f <- z[i + 1L]; z[i + 1L] <- s * z[i] + c * f
        z[i] <- c * z[i] - s * f
      }
      d[l] <- d[l] - p; e[l] <- g; e[m] <- 0
    }
  }
  ord <- order(d)
  list(x = d[ord], w = z[ord]^2, kind = kind)
}

gauss_transform <- function(y, grid) {
  n <- length(y); rec <- basis_recurrence(grid$kind, n)
  prev <- numeric(n); cur <- sqrt(grid$w)
  wy <- cur * y
  out <- rep(0 * y[1L], n)
  for (j in seq_len(n)) {
    out[j] <- sum(cur * wy)
    nextp <- ((grid$x - rec$b[j]) * cur - if (j > 1L) rec$a[j - 1L] * prev else 0) / rec$a[j]
    prev <- cur; cur <- nextp
  }
  out
}

orthogonal_values <- function(c, x, kind) {
  rec <- basis_recurrence(kind, length(c))
  prev <- numeric(length(x)); cur <- rep(1, length(x)); out <- c[1L] * cur
  for (j in seq_len(length(c) - 1L)) {
    nextp <- ((x - rec$b[j]) * cur - if (j > 1L) rec$a[j - 1L] * prev else 0) / rec$a[j]
    out <- out + c[j + 1L] * nextp
    prev <- cur; cur <- nextp
  }
  out
}

coefficient_registry <- function() {
  ordinary_chop <- function(c, tol) {
    cut <- standard_chop(c, tol)
    list(cut = cut, happy = cut < length(c))
  }
  list(
    chebyshev = list(grid = function(n) list(x = sin(pi * seq(0.5, -0.5, length.out = n))),
      transform = cheb_transform, chop = ordinary_chop,
      keep = function(c, cut) c[seq_len(cut)],
      eval = function(c, x) vapply(x, function(t) sum(c * cos((seq_along(c) - 1L) * acos(t))), complex(1)),
      probes = c(-0.3579, 0.1234), start = 17L, next_n = function(n) 2L * n - 1L),
    fourier = list(grid = function(n) list(x = 2 * pi * (0:(n - 1L)) / n),
      transform = fourier_transform, chop = fourier_chop,
      keep = function(c, cut) { d <- length(c) %/% 2L; h <- cut %/% 2L; c[d + 1L + (-h:h)] },
      eval = function(c, x) vapply(x, function(t) sum(c * exp(1i * (seq_along(c) - 1L - length(c) %/% 2L) * t)), complex(1)),
      probes = c(0.1234, 2.3579), start = 17L, next_n = function(n) 2L * n - 1L),
    hermite = list(grid = function(n) gauss_jacobi("hermite", n),
      transform = gauss_transform, chop = ordinary_chop,
      keep = function(c, cut) c[seq_len(cut)],
      eval = function(c, x) orthogonal_values(c, x, "hermite"),
      probes = c(-0.3579, 0.1234), start = 17L, next_n = function(n) 2L * n - 1L),
    laguerre = list(grid = function(n) gauss_jacobi("laguerre", n),
      transform = gauss_transform, chop = ordinary_chop,
      keep = function(c, cut) c[seq_len(cut)],
      eval = function(c, x) orthogonal_values(c, x, "laguerre"),
      probes = c(0.1234, 2.3579), start = 17L, next_n = function(n) 2L * n - 1L))
}

function_coefficients <- function(f, kind, name, tol = .Machine$double.eps, nmax = 65536L) {
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0 || tol >= 1) stopf("function coefficient tol must be in (0, 1)")
  if (length(nmax) != 1L || !is.finite(nmax) || nmax < 17 || nmax != trunc(nmax)) stopf("function coefficient nmax must be a whole number at least 17")
  spec <- coefficient_registry()[[kind]]
  sample <- function_sampler(f, name)
  probes <- sample$eval(spec$probes)
  n <- spec$start; vscale <- max(Mod(probes)); grids <- integer()
  repeat {
    grid <- spec$grid(n)
    y <- sample$eval(grid$x)
    vscale <- max(vscale, Mod(y))
    c <- spec$transform(y, grid)
    if (any(!is.finite(c))) stopf("%s function: %s transform produced non-finite coefficients", name, kind)
    chop <- spec$chop(c, tol)
    grids <- c(grids, n)
    if (chop$happy) {
      kept <- spec$keep(c, chop$cut)
      if (kind != "fourier" && any(Im(kept) != 0)) {
        stopf("%s function must be real in %s", name, kind)
      }
      err <- max(Mod(spec$eval(kept, spec$probes) - probes))
      if (err <= max(100 * .Machine$double.eps, 10 * tol) * max(vscale, .Machine$double.xmin) * length(kept)) {
        if (all(Im(kept) == 0)) kept <- Re(kept)
        discarded <- if (kind == "fourier") {
          h <- length(kept) %/% 2L; d <- length(c) %/% 2L
          c[abs(seq_along(c) - d - 1L) > h]
        } else c[-seq_along(kept)]
        return(list(coefficients = kept, name = name, kind = kind,
                    chop = chop, grids = grids, calls = sample$calls(),
                    scalar = sample$scalar(), vscale = vscale,
                    tail = max(c(0, Mod(discarded))) / max(vscale, .Machine$double.xmin),
                    sample_error = err / max(vscale, .Machine$double.xmin), tol = tol))
      }
    }
    if (n >= nmax) stopf("%s function unresolved in %s at nmax = %d", name, kind, nmax)
    n <- min(spec$next_n(n), nmax)
    if (kind == "fourier" && n %% 2L == 0L) n <- n - 1L
    if (grids[length(grids)] == n) stopf("%s function unresolved in %s at nmax = %d", name, kind, nmax)
  }
}

function_records <- function(A) {
  if (!is_linop(A)) return(list())
  c(A$args$function_coefficients %||% list(),
    unlist(lapply(node_children(A), function_records), recursive = FALSE))
}

basis_polynomial_matrix <- function(v, X, kind, representation, out = matrix(0, nrow(X), ncol(X))) {
  M <- nrow(X)
  prev <- matrix(0, M, M); cur <- diag(1, M)
  rec <- basis_recurrence(kind, length(v))
  for (j in seq_along(v)) {
    if (v[j] != 0) out <- out + v[j] * cur
    nextp <- if (representation == "monomial") cur %*% X else
      ((X - diag(rec$b[j], M)) %*% cur - if (j > 1L) rec$a[j - 1L] * prev else 0) / rec$a[j]
    prev <- cur; cur <- nextp
  }
  out
}

monomial_to_orthogonal <- function(v, kind) {
  n <- length(v)
  X <- basis_matrices(kind, max(n, 2L))$X
  cur <- c(1, numeric(nrow(X) - 1L)); out <- 0 * cur
  for (j in seq_along(v)) {
    out <- out + v[j] * cur
    cur <- as.vector(X %*% cur)
  }
  out[seq_len(n)]
}

space_multiplier <- function(A, orthogonal = FALSE) {
  v <- A$args$v
  if (orthogonal && !identical(A$args$multiplier_basis, "orthogonal")) {
    v <- monomial_to_orthogonal(v, A$args$kind)
  }
  v
}

carry_function_records <- function(A, records) {
  if (length(records)) A$args$function_coefficients <- records
  A
}

