gauss_taylor_coefficients <- function(kind, n, x, value = 0, derivative = 1,
                                      order = if (kind == "laguerre") 64L else 48L) {
  # Hermite: u'' + (2n+1-x^2)u = 0, with u = exp(-x^2/2) H_n.
  # Laguerre: v_ss + ((n+1/2)x-x^2/4)v = 0, with s = log(x), v = exp(-x/2) L_n.
  a <- numeric(order + 1L); a[1L] <- value; a[2L] <- derivative
  at <- function(k) if (k < 0L) 0 else a[k + 1L]
  if (kind == "laguerre") {
    j <- 0:order
    q <- ((n+.5)*x-x^2*2^j/4)/exp(lfactorial(j))
  }
  for (k in 0:(order - 2L)) {
    a[k + 3L] <- if (kind == "hermite") {
      ((x^2 - 2*n - 1)*at(k) + 2*x*at(k-1L) + at(k-2L))/((k+2)*(k+1))
    } else {
      -sum(q[seq_len(k+1L)]*rev(a[seq_len(k+1L)]))/((k+2)*(k+1))
    }
  }
  a
}
gauss_taylor_value <- function(a, t) {
  y <- utils::tail(a, 1L); dy <- 0
  for (j in (length(a)-1L):1L) { dy <- dy*t+y; y <- y*t+a[j] }
  c(value = y, derivative = dy)
}
gauss_taylor_root <- function(a, step, x) {
  converged <- FALSE
  for (iter in 1:12) {
    y <- gauss_taylor_value(a, step)
    change <- y[1L]/y[2L]
    step <- step-change
    if (abs(change) <= 2*.Machine$double.eps*max(1, abs(x+step))) {
      converged <- TRUE
      break
    }
  }
  if (!converged || !is.finite(step) || step <= 0) stop("Gauss Taylor root did not converge")
  list(x = x+step, derivative = gauss_taylor_value(a, step)[2L])
}
gauss_jacobi <- function(kind, n) {
  if (n == 1L) return(list(x = if (kind == "hermite") 0 else 1, w = 1, kind = kind))
  if (kind == "hermite") {
    positive <- numeric(n %/% 2L); logs <- positive
    odd <- n %% 2L == 1L
    x <- 0; logderiv <- 0
    for (j in seq_along(positive)) {
      a <- gauss_taylor_coefficients(kind, n, x, if (j == 1L && !odd) 1 else 0,
                               if (j == 1L && !odd) 0 else 1)
      guess <- (if (j == 1L && !odd) pi/2 else pi)/sqrt(2*n+1-x^2)
      root <- gauss_taylor_root(a, guess, x)
      x <- root$x; logderiv <- logderiv+log(abs(root$derivative))
      positive[j] <- x; logs[j] <- -x^2-2*logderiv
    }
    x <- c(-rev(positive), if (odd) 0, positive)
    logs <- c(rev(logs), if (odd) 0, logs)
  } else {
    x <- numeric(n); logs <- x
    # j_{0,1} is the first positive zero of the Bessel function J_0.
    first <- 2.404825557695773^2/(4*n+2)
    for (iter in 1:8) {
      term <- value <- 1; derivative <- 0
      for (j in seq_len(min(n, 80L))) {
        term <- -term*(n-j+1)*first/j^2
        value <- value+term; derivative <- derivative+j*term/first
        if (abs(term) < .Machine$double.eps^2) break
      }
      first <- first-value/derivative
    }
    x[1L] <- first
    logs[1L] <- log(first)-first
    logderiv <- 0
    for (j in 2:n) {
      a <- gauss_taylor_coefficients(kind, n, x[j-1L])
      Q <- (n+.5)*x[j-1L]-x[j-1L]^2/4
      root <- gauss_taylor_root(a, 2*log1p(pi/(2*sqrt(Q))), log(x[j-1L]))
      x[j] <- exp(root$x); logderiv <- logderiv+log(abs(root$derivative))
      logs[j] <- log(x[j])-x[j]-2*logderiv
    }
  }
  w <- exp(logs-max(logs)); w <- w/sum(w)
  list(x = x, w = w, kind = kind)
}

