analytic_function_fixtures <- function(kind) {
  fs <- if (kind == "chebyshev") list(exp = exp, sin = function(x) sin(pi*x),
    runge = function(x) 1/(1+25*x^2), log = function(x) log(1.1-x),
    polynomial = function(x) 1+2*x+3*x^2)
    else list(exp_cos = function(x) exp(cos(x)), sin = sin,
      rational = function(x) 1/(2-cos(x)), complex_exp = function(x) exp(1i*x))
  set.seed(3801)
  for (seed in 1:24) {
    r <- runif(3, -.75, .75); w <- runif(3, -.5, 1)
    fs[[paste0("seed", seed)]] <- local({
      rates <- r; weights <- w; basis <- kind
      function(x) Reduce(`+`, lapply(seq_along(rates), function(j) {
        t <- if (basis == "fourier") cos(x) else x
        weights[j] * (1-rates[j]^2)/(1-2*rates[j]*t+rates[j]^2)
      }))
    })
    attr(fs[[paste0("seed", seed)]], "rates") <- r
    attr(fs[[paste0("seed", seed)]], "weights") <- w
  }
  fs
}

analytic_error_upper <- function(kind, fixture, f, coefficients) {
  L <- max(if (fixture %in% c("exp", "sin", "exp_cos", "complex_exp")) 64L else 512L,
           length(coefficients))
  k <- 0:L
  bessel_tail <- function(a) {
    degree <- L+1
    log_bound <- log(2) + degree*log(a/2)-lfactorial(degree) +
      a^2/(4*(degree+1))-log1p(-a/(2*(degree+1)))
    max(exp(log_bound), .Machine$double.xmin)
  }
  rates <- attr(f, "rates")
  if (!is.null(rates)) {
    weights <- attr(f, "weights")
    exact <- vapply(k, function(j) sum(weights*rates^j), numeric(1))
    if (kind == "chebyshev") exact[-1L] <- 2*exact[-1L]
    remainder <- 2*sum(abs(weights)*abs(rates)^(L+1)/(1-abs(rates)))
  } else if (kind == "chebyshev") {
    if (fixture == "exp") {
      exact <- 2*besselI(1, k); exact[1L] <- exact[1L]/2
      remainder <- bessel_tail(1)
    } else if (fixture == "sin") {
      exact <- ifelse(k %% 2L == 1L, 2*(-1)^((k-1L) %/% 2L)*besselJ(pi, k), 0)
      remainder <- bessel_tail(pi)
    } else if (fixture == "runge") {
      r <- -(sqrt(26)-1)/(sqrt(26)+1)
      exact <- ifelse(k %% 2L == 0L, 2*r^(k %/% 2L)/sqrt(26), 0)
      exact[1L] <- 1/sqrt(26)
      remainder <- 2*abs(r)^(L %/% 2L+1)/(sqrt(26)*(1-abs(r)))
    } else if (fixture == "log") {
      r <- 1.1-sqrt(1.1^2-1)
      exact <- c(-log(2*r), -2*r^k[-1L]/k[-1L])
      remainder <- 2*r^(L+1)/((L+1)*(1-r))
    } else {
      exact <- c(2.5, 2, 1.5, numeric(L-2L)); remainder <- 0
    }
  } else {
    if (fixture == "exp_cos") {
      exact <- besselI(1, k); remainder <- bessel_tail(1)
    } else if (fixture == "rational") {
      r <- 2-sqrt(3)
      exact <- r^k/sqrt(3)
      remainder <- 2*r^(L+1)/(sqrt(3)*(1-r))
    } else {
      exact <- rep(0i, L+1L)
      exact[2L] <- if (fixture == "sin") -.5i else 1
      remainder <- 0
    }
  }
  if (kind == "fourier") {
    negative <- if (fixture == "complex_exp") rep(0i, L) else Conj(rev(exact[-1L]))
    full <- c(negative, exact)
    error <- linop:::add_centred(full, -coefficients)
  } else error <- exact-linop:::pad_to(coefficients, L+1L)
  sum(Mod(error))+remainder
}
