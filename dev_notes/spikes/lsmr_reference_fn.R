## Fong and Saunders 2011, single column, no damping, no stopping test and no
## preconditioner. Shared by the spikes so they compare against one transcription.
reference_lsmr <- function(M, b, steps) {
  n <- ncol(M)
  Mh <- Conj(t(M))
  x <- rep(0, n)
  if (is.complex(M) || is.complex(b)) storage.mode(x) <- "complex"

  beta <- sqrt(sum(Mod(b)^2)); u <- b / beta
  v <- as.vector(Mh %*% u); alpha <- sqrt(sum(Mod(v)^2)); v <- v / alpha

  zetabar <- alpha * beta
  alphabar <- alpha
  rho <- 1; rhobar <- 1; cbar <- 1; sbar <- 0; zeta <- 0
  h <- v; hbar <- rep(0, n)
  if (is.complex(x)) { storage.mode(h) <- "complex"; storage.mode(hbar) <- "complex" }

  betadd <- beta; betad <- 0; rhodold <- 1; tautildeold <- 0; thetatilde <- 0
  normr <- numeric(steps); normar <- numeric(steps)

  for (k in seq_len(steps)) {
    ut <- as.vector(M %*% v) - alpha * u
    beta <- sqrt(sum(Mod(ut)^2))
    u <- if (beta > 0) ut / beta else ut
    vt <- as.vector(Mh %*% u) - beta * v
    alpha_new <- sqrt(sum(Mod(vt)^2))

    rhoold <- rho
    rho <- sqrt(alphabar^2 + beta^2)
    cs <- alphabar / rho
    sn <- beta / rho
    thetanew <- sn * alpha_new
    alphabar <- cs * alpha_new

    rhobarold <- rhobar
    zetaold <- zeta
    thetabar <- sbar * rho
    rhotemp <- cbar * rho
    rhobar <- sqrt(rhotemp^2 + thetanew^2)
    cbar <- rhotemp / rhobar
    sbar <- thetanew / rhobar
    zeta <- cbar * zetabar
    zetabar <- -sbar * zetabar

    hbar <- h - (thetabar * rho / (rhoold * rhobarold)) * hbar
    x <- x + (zeta / (rho * rhobar)) * hbar
    v <- if (alpha_new > 0) vt / alpha_new else vt
    h <- v - (thetanew / rho) * h
    alpha <- alpha_new

    betahat <- cs * betadd
    betadd <- -sn * betadd
    thetatildeold <- thetatilde
    rhotildeold <- sqrt(rhodold^2 + thetabar^2)
    ctildeold <- rhodold / rhotildeold
    stildeold <- thetabar / rhotildeold
    thetatilde <- stildeold * rhobar
    rhodold <- ctildeold * rhobar
    betad <- -stildeold * betad + ctildeold * betahat
    tautildeold <- (zetaold - thetatildeold * tautildeold) / rhotildeold
    taud <- (zeta - thetatilde * tautildeold) / rhodold

    normr[k] <- sqrt((betad - taud)^2 + betadd^2)
    normar[k] <- abs(zetabar)
  }
  list(x = matrix(x, n, 1L), normr = normr, normar = normar)
}
