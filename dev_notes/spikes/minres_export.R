## Write the fixtures and this implementation's iterates for the SciPy
## cross-check, and read the comparison back. The external half of Gate 2's
## reference line for MINRES; the dense-definition half is minres_reference.R.
##
##   Rscript dev_notes/spikes/minres_export.R write
##   python  dev_notes/spikes/minres_scipy.py
##   Rscript dev_notes/spikes/minres_export.R compare
##
## Real only. SciPy's minres casts a complex operator to real with a
## ComplexWarning and solves a different problem, so the complex fixtures are the
## dense reference's alone.

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

dir <- file.path(tempdir(), "minres-xcheck")
dir <- Sys.getenv("MINRES_XDIR", dir)
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
cat("dir:", dir, "\n")

indef_spectrum <- function(n, kappa) {
  rep(c(-1, 1), length.out = n) * exp(seq(log(1), log(1 / kappa), length.out = n))
}

## Four conditionings, an exact breakdown and a near-breakdown, so the comparison
## covers the regime where the recurrence is still reproducible, the regime where
## it is not, and the two cases the gate line names.
cases <- list(
  list(name = "kappa1e2",  kind = "indef", kappa = 1e2,  n = 40, seed = 1),
  list(name = "kappa1e4",  kind = "indef", kappa = 1e4,  n = 40, seed = 2),
  list(name = "kappa1e6",  kind = "indef", kappa = 1e6,  n = 40, seed = 3),
  list(name = "kappa1e10", kind = "indef", kappa = 1e10, n = 40, seed = 4),
  list(name = "breakdown", kind = "invariant", n = 60, d = 3L, seed = 5),
  list(name = "nearbreak", kind = "perturbed", n = 60, delta = 1e-10, seed = 6))
steps <- c(1L, 2L, 3L, 4L, 6L, 8L, 12L, 16L)

build <- function(cs) {
  if (identical(cs$kind, "indef")) {
    M <- spd_prescribed(cs$n, indef_spectrum(cs$n, cs$kappa), seed = cs$seed)
    set.seed(5000L + cs$seed)
    return(list(M = M, b = matrix(stats::rnorm(cs$n), cs$n, 1L)))
  }
  ## Both remaining kinds are the indefinite shifted Laplacian, whose
  ## eigenvectors are closed form, so the invariant subspace is constructed.
  M <- as.matrix(shifted_laplacian_1d(cs$n, 1.5))
  V <- laplacian_1d_eigenvectors(cs$n)
  if (identical(cs$kind, "invariant")) {
    set.seed(5000L + cs$seed)
    b <- V[, seq_len(cs$d) * 5L, drop = FALSE] %*% stats::rnorm(cs$d)
  } else {
    b <- V[, 7L] + cs$delta * V[, 23L]
  }
  list(M = M, b = matrix(b, cs$n, 1L))
}

## The published recurrence, transcribed straight from the Paige and Saunders
## form with the preconditioner set to the identity and no shift. This is the
## reference the test suite carries, so what SciPy is for is checking the
## transcription: a reference written here and never compared to anything is a
## second copy of the same beliefs.
reference_minres <- function(M, b, steps) {
  n <- length(b)
  eps <- .Machine$double.eps
  x <- rep(0, n)
  r1 <- b
  y <- r1
  beta1 <- sqrt(sum(r1 * y))
  if (beta1 <= 0) return(x)
  oldb <- 0; beta <- beta1; dbar <- 0; epsln <- 0
  phibar <- beta1
  cs <- -1; sn <- 0
  w <- rep(0, n); w2 <- rep(0, n)
  r2 <- r1
  for (itn in seq_len(steps)) {
    if (beta <= 0) break
    v <- y / beta
    y <- as.vector(M %*% v)
    if (itn >= 2) y <- y - (beta / oldb) * r1
    alfa <- sum(v * y)
    y <- y - (alfa / beta) * r2
    r1 <- r2; r2 <- y
    oldb <- beta
    beta <- sqrt(max(sum(r2 * y), 0))
    oldeps <- epsln
    delta <- cs * dbar + sn * alfa
    gbar <- sn * dbar - cs * alfa
    epsln <- sn * beta
    dbar <- -cs * beta
    gamma <- max(sqrt(gbar^2 + beta^2), eps)
    cs <- gbar / gamma
    sn <- beta / gamma
    phi <- cs * phibar
    phibar <- sn * phibar
    w1 <- w2; w2 <- w
    w <- (v - oldeps * w1 - delta * w2) / gamma
    x <- x + phi * w
  }
  x
}

action <- commandArgs(trailingOnly = TRUE)[1L]

if (identical(action, "write")) {
  for (cs in cases) {
    d <- build(cs)
    utils::write.table(d$M, file.path(dir, paste0(cs$name, "-A.txt")),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(d$b, file.path(dir, paste0(cs$name, "-b.txt")),
                       row.names = FALSE, col.names = FALSE)
    for (k in steps) {
      x <- linop:::minres_solve(linop(d$M, properties = c(hermitian = TRUE)), d$b,
                                tol = 0, maxit = k)$x
      utils::write.table(x, file.path(dir, sprintf("%s-mine-%d.txt", cs$name, k)),
                         row.names = FALSE, col.names = FALSE)
    }
  }
  writeLines(as.character(steps), file.path(dir, "steps.txt"))
  writeLines(vapply(cases, function(cs) cs$name, character(1)),
             file.path(dir, "cases.txt"))
  cat("written\n")
}

if (identical(action, "compare")) {
  cat("  mine   this implementation against SciPy 1.17.1\n")
  cat("  trans  the transcription the test suite carries, against SciPy\n\n")
  cat(sprintf("%-10s %6s  %11s  %11s  %11s  %9s\n",
              "case", "steps", "mine", "trans", "||x||", "scipy its"))
  for (cs in cases) {
    d <- build(cs)
    for (k in steps) {
      fm <- file.path(dir, sprintf("%s-mine-%d.txt", cs$name, k))
      fs <- file.path(dir, sprintf("%s-scipy-%d.txt", cs$name, k))
      if (!file.exists(fs)) next
      mine <- as.matrix(utils::read.table(fm))
      theirs <- as.matrix(utils::read.table(fs))
      trans <- matrix(reference_minres(d$M, as.vector(d$b), k), ncol = 1L)
      its <- scan(file.path(dir, sprintf("%s-scipyit-%d.txt", cs$name, k)),
                  quiet = TRUE)
      nx <- max(Mod(theirs))
      rel <- function(u) if (nx > 0) max(Mod(u - theirs)) / nx else max(Mod(u - theirs))
      cat(sprintf("%-10s %6d  %11.3e  %11.3e  %11.3e  %9d\n", cs$name, k,
                  rel(mine), rel(trans), nx, as.integer(its)))
    }
  }
}
