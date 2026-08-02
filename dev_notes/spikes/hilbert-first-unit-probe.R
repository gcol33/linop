## The first Hilbert unit, built against v0.1 core through the exported surface only.
##
## The open question blocking the package decision is "does the first unit need any
## new core export?", and it is not answerable by reading NAMESPACE: the answer is
## whichever names the unit turns out to reach for. So this builds the unit -- the
## S0.6 class, finite section of a discrete Schrodinger operator on ell^2(Z) -- as a
## provider would, using only what a separate package could call, and records where
## the surface runs out.
##
## What it is not: an implementation. The certificate assembled here is a list of
## numbers, not a linop_certificate, because that is one of the findings.
##
## Run from the repository root:
##   "C:/Program Files/R/R-4.6.0/bin/x64/Rscript.exe" dev_notes/spikes/hilbert-first-unit-probe.R

devtools::load_all(".", quiet = TRUE)
options(digits = 12)

out_dir <- file.path("dev_notes", "spikes", "results")
if (!dir.exists(out_dir)) stop("run from the repository root")
txt <- file(file.path(out_dir, "hilbert-first-unit.txt"), open = "wt")
say <- function(...) {
  l <- paste0(...)
  cat(l, "\n", sep = "")
  cat(l, "\n", sep = "", file = txt)
}

## ==========================================================================
## The unit. (H u)_j = u_{j-1} + u_{j+1} + V_j u_j on ell^2(Z), V real with
## support in [-R, R]. The finite section keeps |j| <= n and is tridiagonal
## with unit off-diagonals, so it is hermitian by construction and the apply
## is three shifted adds.
## ==========================================================================

finite_section <- function(V, n) {
  R <- (length(V) - 1L) / 2L
  if (R != as.integer(R)) stop("V must have odd length, giving V_j for j = -R..R")
  if (n <= R) stop("n must exceed the support radius so the tail is free")
  N <- 2L * n + 1L
  d <- numeric(N)
  d[(n - R + 1L):(n + R + 1L)] <- V

  apply_fs <- function(X) {
    X <- as.matrix(X)
    rbind(X[-1L, , drop = FALSE], 0) + rbind(0, X[-nrow(X), , drop = FALSE]) + d * X
  }

  H <- linop(apply_fs, adjoint = apply_fs, dim = c(N, N),
             properties = c("hermitian", "symmetric"))
  set_provenance(H, "linop.hilbert",
                 structure(list(V = V, R = R, n = n, band_edge = 2,
                                norm_bound = 2 + max(abs(V))),
                           class = "finite_section"))
}

## The provider's methods, registered against the payload class.
registerS3method("provenance_summary", "finite_section", function(p, ...) {
  sprintf("finite section of a discrete Schrodinger operator, n = %d, supp V = [-%d, %d]",
          p$payload$n, p$payload$R, p$payload$R)
}, envir = globalenv())
registerS3method("provenance_refine", "finite_section", function(p, n_new, ...) {
  finite_section(p$payload$V, n_new)
}, envir = globalenv())
registerS3method("provenance_lift", "finite_section", function(p, x, ...) {
  ## the finite vector as an element of ell^2(Z): index it
  data.frame(j = seq.int(-p$payload$n, p$payload$n), value = as.numeric(x))
}, envir = globalenv())
registerS3method("provenance_original_residual", "finite_section",
                 function(p, result, ...) {
  ## (H - lambda) u~ is nonzero only at j = +-(n+1), where it takes the values
  ## u_n and u_{-n}. Two entries of the computed vector and nothing else.
  u <- as.numeric(result$vectors[, 1L])
  sqrt(u[1L]^2 + u[length(u)]^2) / sqrt(sum(u^2))
}, envir = globalenv())

## --------------------------------------------------------- closed forms ---
rho_of <- function(lambda) { a <- abs(lambda) / 2; a - sqrt(a^2 - 1) }
lambda_delta <- function(v) sign(v) * sqrt(v^2 + 4)

## ==========================================================================
## 1. What core gives a matrix-free provider that knows its own structure
## ==========================================================================
say("=== 1. the capability a provider can state, and the one it can prove ===")

H <- finite_section(c(1), 20)
say("  matrix-free leaf, properties = c('hermitian'):")
say("    value  : ", format(cap(H, "hermitian")$value))
say("    source : ", cap(H, "hermitian")$evidence$source)
say("    guarantee: ", cap(H, "hermitian")$evidence$guarantee)

Hd <- linop(as.matrix(H))
say("  the same operator as a dense leaf, with no properties given:")
say("    value  : ", format(cap(Hd, "hermitian")$value))
say("    source : ", cap(Hd, "hermitian")$evidence$source)
say("    guarantee: ", cap(Hd, "hermitian")$evidence$guarantee)
say("")
say("  So the route that materialises the operator gets the better evidence, and")
say("  the matrix-free route -- the one the layer exists for -- gets a bare")
say("  declaration. normalise_properties() stamps ev_declared() unconditionally;")
say("  ev_construction() is internal and there is no properties= form reaching it.")
say("")
vcert <- verify(H)
say("  verify(H) on the matrix-free leaf: ", vcert$overall)
say("    checks: ", paste(vcert$checks$check, collapse = ", "))
say("")

## ==========================================================================
## 2. The provenance envelope, obtained the way a caller obtains one
## ==========================================================================
say("=== 2. provenance round trip through the caller's path ===")
p <- provenance(H)
say("  summary   : ", provenance_summary(p))
say("  lift      : ", nrow(provenance_lift(p, rep(1, 41))), " indexed entries, j from ",
    provenance_lift(p, rep(1, 41))$j[1L])
say("  refine    : dim of refine(p, 30) = ",
    paste(dim(provenance_refine(p, 30)), collapse = " x "))
say("  after t(H) %*% H: ", provenance_summary(provenance(t(H) %*% H)))
say("")

## ==========================================================================
## 3. eigs() against the closed form, with the three-part certificate's numbers
## ==========================================================================
say("=== 3. eigs() on V = 1 * delta_0, exact lambda = sqrt(5) ===")
v <- 1
lam_exact <- lambda_delta(v)
rho_exact <- rho_of(lam_exact)
eps <- .Machine$double.eps
say(sprintf("  lambda_exact = %.12f   rho_exact = %.12f", lam_exact, rho_exact))
say("     n   lambda_n         true err     eta          floor       err<=eta+floor  eta/rho^n  eigs")

rows <- list()
for (n in c(5, 10, 20, 40, 80)) {
  Hn <- finite_section(c(v), n)
  fit <- eigs(Hn, k = 1L, which = "largest", tol = 1e-12, ncv = min(2L * n + 1L, 40L))
  q <- fit$values[1L]
  eta <- provenance_original_residual(provenance(Hn), fit)
  floor_v <- 4 * provenance(Hn)$payload$norm_bound * eps
  err <- abs(q - lam_exact)
  ok <- err <= eta + floor_v
  say(sprintf("  %4d  %.12f  %.3e   %.3e   %.3e   %-5s        %.4f     %s",
              n, q, err, eta, floor_v, ok, eta / rho_exact^n,
              fit$certificate$overall))
  rows[[length(rows) + 1L]] <- data.frame(n = n, lambda = q, err = err, eta = eta,
                                          floor = floor_v, holds = ok,
                                          ratio = eta / rho_exact^n,
                                          eigs = fit$certificate$overall,
                                          nconv = fit$nconv)
}
say("")
say("  Without the floor the n = 80 row is the S0.6 failure: eta has decayed below")
say("  the arithmetic plateau, so a converged answer certifies as fail. The floor")
say("  is the provider's own, 2 + max|V| in closed form, no norm estimate needed.")
say("")

## ==========================================================================
## 4. Kato-Temple against the known band edge
## ==========================================================================
say("=== 4. Kato-Temple, band edge a = 2 known from the symbol ===")
say("     n   eta (linear)   KT half-width   contains sqrt(5)")
for (n in c(5, 10, 20, 40)) {
  Hn <- finite_section(c(v), n)
  fit <- eigs(Hn, k = 1L, which = "largest", tol = 1e-12, ncv = min(2L * n + 1L, 40L))
  q <- fit$values[1L]
  eta <- provenance_original_residual(provenance(Hn), fit)
  a <- provenance(Hn)$payload$band_edge
  kt <- eta^2 / (q - a)
  contains <- (lam_exact >= q - kt - 1e-14) && (lam_exact <= q + kt + 1e-14)
  say(sprintf("  %4d   %.6e   %.6e    %s", n, eta, kt, contains))
}
say("")

## ==========================================================================
## 5. The refusal, on q - a > 0 rather than on a separate check
## ==========================================================================
say("=== 5. refusal at V = 0: the isolation condition and the refusal are one ===")
say("     n    q (top Ritz)      q - a         issue?  eigs nconv  eigs cert")
for (n in c(10, 40, 100)) {
  Hn <- finite_section(c(0), n)
  fit <- eigs(Hn, k = 1L, which = "largest", tol = 1e-12, ncv = min(2L * n + 1L, 40L))
  q <- fit$values[1L]
  gap <- q - 2
  say(sprintf("  %4d   %.12f   %.6e   %-5s   %d of 1      %s",
              n, q, gap, gap > 0, fit$nconv, fit$certificate$overall))
}
say("  Every top Ritz value sits strictly inside (-2, 2): these are discretisations")
say("  of continuous spectrum and H_free has no eigenvalues at all. q - a > 0 fails")
say("  on its own terms, so nothing separate has to detect the case.")
say("")
say("  The eigensolver is a second signal and not a substitute for the first. Where")
say("  it converges, the value it returns is still inside the band and the layer")
say("  still has to refuse; where it stalls, the stall is about the clustering near")
say("  the band edge rather than about the absence of an eigenvalue. Only q - a > 0")
say("  distinguishes 'no eigenvalue here' from 'a hard eigenvalue here'.")
say("")

## ==========================================================================
## 6. What the unit reached for, and whether a separate package could reach it
## ==========================================================================
say("=== 6. the exported surface, measured ===")
wanted <- c("linop", "set_provenance", "provenance", "strip_provenance",
            "provenance_lift", "provenance_refine", "provenance_original_residual",
            "provenance_summary", "eigs", "verify", "capability", "evidence",
            "cap", "cap_value", "requirement", "evidence_satisfies",
            "build_certificate", "cert_rows", "solve_certificate",
            "arithmetic_floor", "norm2", "ev_construction", "linop_apply",
            "SOLVE_FLOOR_CONST")
exported <- getNamespaceExports("linop")
tbl <- data.frame(name = wanted, exported = wanted %in% exported,
                  stringsAsFactors = FALSE)
tbl$used <- tbl$name %in% c("linop", "set_provenance", "provenance", "eigs",
                            "provenance_lift", "provenance_refine",
                            "provenance_original_residual", "provenance_summary",
                            "cap", "verify")
for (i in seq_len(nrow(tbl))) {
  say(sprintf("  %-30s exported: %-5s  used by the unit: %s",
              tbl$name[i], tbl$exported[i], tbl$used[i]))
}
say("")
missing_used <- tbl$name[tbl$used & !tbl$exported]
say("  reached and not exported: ",
    if (length(missing_used)) paste(missing_used, collapse = ", ") else "none")
say("  needed to build a linop_certificate: ",
    paste(tbl$name[tbl$name %in% c("build_certificate", "cert_rows",
                                   "arithmetic_floor", "norm2") &
                   !tbl$exported], collapse = ", "))

utils::write.csv(do.call(rbind, rows),
                 file.path(out_dir, "hilbert-first-unit.csv"), row.names = FALSE)
utils::write.csv(tbl, file.path(out_dir, "hilbert-first-unit-surface.csv"),
                 row.names = FALSE)
close(txt)
