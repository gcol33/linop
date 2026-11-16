coefficient_shells <- function(c, kind) {
  if (kind != "fourier") return(Mod(c))
  d <- length(c) %/% 2L
  c(Mod(c[d + 1L]), pmax(Mod(c[d + 1L - seq_len(d)]),
                         Mod(c[d + 1L + seq_len(d)])))
}

coefficient_discarded <- function(full, kept, kind) {
  if (kind != "fourier") return(full[-seq_along(kept)])
  full[abs(seq_along(full) - length(full) %/% 2L - 1L) > length(kept) %/% 2L]
}

coefficient_contour <- function(record) {
  if (!is.null(record$continuation$value)) return(record$continuation$value)
  result <- function(value) { record$continuation$value <- value; value }
  unavailable <- function(reason) result(list(available = FALSE, reason = reason))
  kind <- record$kind
  if (!kind %in% c("chebyshev", "fourier")) {
    return(unavailable("ellipse and strip estimates apply to the interval and circle"))
  }
  s <- coefficient_shells(record$full_coefficients, kind)
  active <- which(s > max(s) * 1e-10) - 1L
  alpha <- log(16)
  if (length(active) >= 4L) {
    k <- utils::tail(active, ceiling(length(active) / 2))
    centered <- k - mean(k)
    slope <- -sum(centered * log(s[k + 1L])) / sum(centered^2)
    if (!is.finite(slope) || slope <= 0) return(unavailable("geometric coefficient decay was not detected"))
    alpha <- min(log(4), slope * .8)
  }
  if (alpha * (length(s) - 1L) < 1) {
    return(unavailable("coefficient decay is too flat to estimate an analytic domain"))
  }
  m <- 2^ceiling(log2(max(256, 128 / alpha)))
  theta <- 2 * pi * (0:(m - 1L)) / m
  rho <- exp(alpha)
  z <- if (kind == "chebyshev") (rho * exp(1i * theta) + exp(-1i * theta) / rho) / 2
       else c(theta + 1i * alpha, theta - 1i * alpha)
  y <- tryCatch(record$fun(z), error = identity)
  if (inherits(y, "error") || !(is.numeric(y) || is.complex(y)) ||
      !is.null(dim(y)) || length(y) != length(z) || any(!is.finite(y))) {
    return(unavailable("function does not accept the complex contour in one vector call"))
  }
  M <- 2 * max(Mod(y))
  k <- seq_len(min(8L, length(s) - 1L, floor(log(1e4)/alpha)))
  if (kind == "chebyshev") {
    c <- stats::fft(y) / m
    plus <- 2 * c[k + 1L] * exp(-alpha * k)
    minus <- 2 * c[m - k + 1L] * exp(alpha * k)
    expected <- record$full_coefficients[k + 1L]
  } else {
    plus <- stats::fft(y[seq_len(m)]) / m
    minus <- stats::fft(y[m + seq_len(m)]) / m
    plus <- plus[k + 1L] * exp(alpha * k)
    minus <- minus[k + 1L] * exp(-alpha * k)
    d <- length(record$full_coefficients) %/% 2L
    expected <- record$full_coefficients[d + 1L + k]
  }
  if (max(Mod(c(plus - expected, minus - expected))) > 1e-7 * max(M, record$vscale, .Machine$double.xmin)) {
    return(unavailable("complex samples are inconsistent with analytic continuation of the real samples"))
  }
  degree <- length(s) - 1L
  remainder <- if (kind == "chebyshev") 4 * M * exp(-alpha * degree) / expm1(alpha)
               else 4 * M * exp(-alpha * degree) / (-expm1(-alpha))
  result(list(available = TRUE, alpha = alpha, rho = rho, M = M, samples = length(z),
              interpolation = remainder))
}

coefficient_tail_bound <- function(record) {
  contour <- coefficient_contour(record)
  tail <- sum(Mod(coefficient_discarded(record$full_coefficients, record$coefficients, record$kind)))
  roundoff <- 8 * .Machine$double.eps * (1 + log2(utils::tail(record$grids, 1L))) *
    max(record$vscale, sum(Mod(record$full_coefficients)))
  contour$discarded <- tail
  contour$roundoff <- roundoff
  contour$bound <- if (contour$available) contour$interpolation + tail + roundoff else NA_real_
  contour
}

scale_function_records <- function(records, a) {
  lapply(records, function(x) {
    if (!startsWith(x$name, "boundary")) x$factor <- (x$factor %||% 1) * a
    x
  })
}

function_request_operator <- function(A, tol, floor_const = SOLVE_FLOOR_CONST) {
  records <- A$args$function_coefficients
  if (!length(records) || !A$node %in% c("fourier", "chebyshev")) return(A)
  potential <- which(vapply(records, function(x)
    x$name %in% c("multiplier", "derivative[[1]]"), logical(1)))
  if (!length(potential)) return(A)
  scale <- if (A$node == "fourier") A$args$vnorm + Mod(utils::tail(A$args$coef, 1L)) else 1
  budget <- tol * scale / (10 * length(potential))
  mult <- if (A$node == "fourier") A$args$vhat else A$args$terms[[1L]]
  changed <- FALSE
  add <- if (A$node == "fourier") add_centred else add_padded
  for (i in potential) {
    z <- records[[i]]
    factor <- z$factor %||% 1
    if (Mod(factor) == 0) next
    relative <- min(.1, max(z$tol, budget / (Mod(factor) * max(z$vscale, .Machine$double.xmin))))
    spec <- coefficient_registry()[[z$kind]]
    chop <- spec$chop(z$full_coefficients, relative)
    if (!chop$happy) next
    kept <- spec$keep(z$full_coefficients, chop$cut)
    if (length(kept) >= length(z$coefficients)) next
    candidate <- z
    candidate$coefficients <- kept
    bound <- coefficient_tail_bound(candidate)
    if (!bound$available || Mod(factor) * bound$bound > budget) next
    mult <- add(mult, factor * add(kept, -z$coefficients))
    candidate$tol <- relative
    candidate$chop <- chop
    candidate$tail <- max(c(0, Mod(coefficient_discarded(z$full_coefficients, kept, z$kind)))) /
      max(z$vscale, .Machine$double.xmin)
    candidate$sample_error <- z$sample_error + sum(Mod(add(kept, -z$coefficients))) /
      max(z$vscale, .Machine$double.xmin)
    records[[i]] <- candidate
    changed <- TRUE
  }
  if (!changed) return(A)
  mult <- if (all(Im(mult) == 0)) Re(mult) else mult
  if (A$node == "fourier") {
    while (length(mult) > 1L && mult[1L] == 0 && utils::tail(mult, 1L) == 0) mult <- mult[-c(1L, length(mult))]
    out <- new_space_coefficients("fourier", A$args$coef, mult, NULL)
    n <- max(BAND_N_START, out$args$d)
    norm <- max(Mod(fourier_symbol(out, seq.int(-n, n)))) + out$args$vnorm
    floor <- floor_const * .Machine$double.eps * norm
    if (floor > tol * fourier_scale(out, 0) / 4) return(A)
  } else {
    terms <- A$args$terms
    terms[[1L]] <- trim_poly(mult)
    out <- new_chebyshev(terms, 0, A$args$boundary)
  }
  carry_function_records(out, records)
}
