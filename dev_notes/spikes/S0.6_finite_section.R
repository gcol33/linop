## S0.6 -- finite-section bound feasibility for the discrete Schrodinger operator
##   (H u)_j = u_{j-1} + u_{j+1} + V_j u_j    on ell^2(Z),  V real, finite support
##
## Verifies numerically the three closed forms derived in dev_notes/S0.6-*.md:
##   (a) decay rate      rho(lambda) = |lambda|/2 - sqrt((lambda/2)^2 - 1)
##   (b) a posteriori    dist(lambda_n, spec H) <= ||r|| / ||u||, r supported on 2 sites
##   (c) Kato-Temple     quadratic improvement using the known band edge +-2
##
## Exact ground truth used: single-site potential V = v*delta_0 has exactly one
## eigenvalue outside [-2,2], at lambda = sign(v)*sqrt(v^2+4), with eigenvector
## u_j = rho^|j|, rho = (-|v| + sqrt(v^2+4))/2.

options(digits = 12)
cat("S0.6 -- finite-section bounds, discrete Schrodinger on ell^2(Z)\n")
cat("R:", R.version.string, "\n\n")

## --- closed forms ----------------------------------------------------------
rho_of  <- function(lambda) { a <- abs(lambda) / 2; a - sqrt(a^2 - 1) }   # in (0,1) for |lambda|>2
alpha_of <- function(lambda) acosh(abs(lambda) / 2)                       # rho = exp(-alpha)

## --- finite section H_n on indices -n..n (dimension N = 2n+1) --------------
make_Hn <- function(n, Vfun) {
  N <- 2L * n + 1L
  idx <- seq.int(-n, n)
  H <- matrix(0, N, N)
  diag(H) <- Vfun(idx)
  for (i in seq_len(N - 1L)) { H[i, i + 1L] <- 1; H[i + 1L, i] <- 1 }
  list(H = H, idx = idx)
}

## --- a posteriori residual bound -------------------------------------------
## Extend the finite eigenvector by zero. (H - lam) u~ is nonzero only at
## j = +-(n+1), with values u_n and u_{-n}. So ||r|| = sqrt(u_n^2 + u_{-n}^2).
resid_bound <- function(vec, n) {
  N <- length(vec)
  sqrt(vec[1]^2 + vec[N]^2) / sqrt(sum(vec^2))
}

## ===========================================================================
## 1. exact ground truth, single-site potential
## ===========================================================================
cat("=== 1. single-site potential V = v*delta_0: exact lambda = sqrt(v^2+4) ===\n")
for (v in c(0.5, 1, 2, 4)) {
  lam_exact <- sqrt(v^2 + 4)
  rho_exact <- (-v + sqrt(v^2 + 4)) / 2
  cat(sprintf("\n v = %-4g  lambda_exact = %.12f   rho_exact = %.12f   alpha = %.6f\n",
              v, lam_exact, rho_exact, alpha_of(lam_exact)))
  cat(sprintf("   check rho from lambda: %.12f  (match: %s)\n",
              rho_of(lam_exact), isTRUE(all.equal(rho_of(lam_exact), rho_exact))))
  Vf <- function(j) ifelse(j == 0L, v, 0)
  cat("     n    lambda_n           err=|lam_n-lam|   ||r||/||u||   err<=bound   bound/rho^n\n")
  for (n in c(5, 10, 20, 40, 80)) {
    hh <- make_Hn(n, Vf)
    e  <- eigen(hh$H, symmetric = TRUE)
    k  <- which.max(e$values)
    lam_n <- e$values[k]; vec <- e$vectors[, k]
    err <- abs(lam_n - lam_exact)
    bnd <- resid_bound(vec, n)
    cat(sprintf("  %4d  %.12f  %.3e     %.3e     %-5s      %.4f\n",
                n, lam_n, err, bnd, err <= bnd * (1 + 1e-9), bnd / rho_exact^n))
  }
}

## ===========================================================================
## 2. decay rate of the computed eigenvector matches rho(lambda)
## ===========================================================================
cat("\n=== 2. measured decay of the eigenvector vs closed-form rho ===\n")
v <- 1; Vf <- function(j) ifelse(j == 0L, v, 0)
n <- 60; hh <- make_Hn(n, Vf)
e <- eigen(hh$H, symmetric = TRUE); k <- which.max(e$values)
vec <- abs(e$vectors[, k]); lam_n <- e$values[k]
ratios <- vec[(n + 2):(n + 25)] / vec[(n + 1):(n + 24)]   # u_{j+1}/u_j for j = 0..23
cat(sprintf("  closed-form rho(lambda_n) = %.12f\n", rho_of(lam_n)))
cat(sprintf("  measured ratios u_{j+1}/u_j, j=1..10: %s\n",
            paste(sprintf("%.9f", ratios[2:11]), collapse = " ")))
cat(sprintf("  max |measured - closed form| over j=1..20: %.3e\n",
            max(abs(ratios[2:21] - rho_of(lam_n)))))

## ===========================================================================
## 3. Kato-Temple, using the known essential-spectrum edge +2
## ===========================================================================
cat("\n=== 3. Kato-Temple bracket (band edge a = 2 is known from the symbol) ===\n")
cat("   lambda in [q - eta^2/(b-q),  q + eta^2/(q-a)] with q the Rayleigh quotient\n")
lam_exact <- sqrt(v^2 + 4)
cat("     n   eta=||r||/||u||   simple |err|<=eta    KT half-width      KT contains truth\n")
for (n in c(5, 10, 20, 40)) {
  hh <- make_Hn(n, Vf); e <- eigen(hh$H, symmetric = TRUE)
  k <- which.max(e$values); vec <- e$vectors[, k]; vec <- vec / sqrt(sum(vec^2))
  q <- as.numeric(t(vec) %*% hh$H %*% vec)
  eta <- resid_bound(vec, n)
  a <- 2                                   # top of essential spectrum
  kt <- eta^2 / (q - a)                    # upper-side Kato-Temple width
  contains <- (lam_exact >= q - kt - 1e-14) && (lam_exact <= q + kt + 1e-14)
  cat(sprintf("  %4d   %.6e     %.6e      %.6e     %s\n", n, eta, eta, kt, contains))
}

## ===========================================================================
## 4. the refusal case: V = 0, no eigenvalues outside [-2,2]
## ===========================================================================
cat("\n=== 4. refusal case V = 0: spectrum is purely continuous on [-2,2] ===\n")
for (n in c(10, 40, 100)) {
  hh <- make_Hn(n, function(j) rep(0, length(j)))
  ev <- eigen(hh$H, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf("  n=%3d  N=%4d  range [%.9f, %.9f]  any |lam|>2? %s   max|lam| - 2 = %.3e\n",
              n, 2 * n + 1, min(ev), max(ev), any(abs(ev) > 2), max(abs(ev)) - 2))
}
cat("  -> every finite-section eigenvalue lies strictly inside (-2,2); they are\n")
cat("     discretisations of continuous spectrum, not eigenvalues. eigs() must refuse.\n")

## ===========================================================================
## 5. multi-site potential, no closed form: bound must still be computable
## ===========================================================================
cat("\n=== 5. three-site potential V = (1.5, -0.5, 2.0) at j = -1,0,1 ===\n")
Vf3 <- function(j) ifelse(j == -1L, 1.5, ifelse(j == 0L, -0.5, ifelse(j == 1L, 2.0, 0)))
ref <- eigen(make_Hn(400, Vf3)$H, symmetric = TRUE, only.values = TRUE)$values
lam_ref <- max(ref)
cat(sprintf("  reference (n=400) top eigenvalue: %.12f   outside band: %s\n",
            lam_ref, abs(lam_ref) > 2))
cat("     n    lambda_n           |lam_n - lam_ref|   eta          err<=eta\n")
for (n in c(5, 10, 20, 40, 80)) {
  hh <- make_Hn(n, Vf3); e <- eigen(hh$H, symmetric = TRUE)
  k <- which.max(e$values); lam_n <- e$values[k]; vec <- e$vectors[, k]
  err <- abs(lam_n - lam_ref); eta <- resid_bound(vec, n)
  cat(sprintf("  %4d  %.12f  %.6e      %.6e   %s\n", n, lam_n, err, eta,
              err <= eta * (1 + 1e-9)))
}

## ===========================================================================
## 6. what must the discretisation object carry?
## ===========================================================================
cat("\n=== 6. inputs needed for each bound ===\n")
cat("  rho(lambda)      <- lambda only                          [closed form]\n")
cat("  eta = ||r||/||u|| <- boundary entries u_{-n}, u_n of the computed vector\n")
cat("  Kato-Temple      <- eta, Rayleigh quotient q, band edge a=2 (from the symbol)\n")
cat("  band edge        <- limiting off-diagonal 1 and diagonal 0 => [-2,2] by Weyl\n")
cat("  support radius R <- needed only to assert n > R so the tail is free\n")
