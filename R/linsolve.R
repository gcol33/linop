## Section 4.2. A linsolve computes an approximation to A^-1 b and may be
## adaptive. It is NOT a linop: an object is a linop if and only if its action on
## x is determined without reference to x (section 4.1), and CG chooses its
## polynomial by minimising over the Krylov space of b.
##
## Internal in v0.1: the public surface exposes solve() only (section 1.1).
## Designing the class now does not oblige exposing it now.

LINSOLVE_FIDELITY    <- c("exact", "linear_approximation", "variable_inexact")
LINSOLVE_DETERMINACY <- c("fixed", "input_dependent", "history_dependent")
LINSOLVE_RANDOMNESS  <- c("deterministic", "stochastic")

new_linsolve <- function(apply_inverse, dim, fidelity, determinacy,
                         randomness = "deterministic",
                         adjoint = NULL, error_bound = FALSE, method = "unknown") {
  if (!fidelity %in% LINSOLVE_FIDELITY) {
    stopf("fidelity must be one of: %s", paste(LINSOLVE_FIDELITY, collapse = ", "))
  }
  if (!determinacy %in% LINSOLVE_DETERMINACY) {
    stopf("determinacy must be one of: %s", paste(LINSOLVE_DETERMINACY, collapse = ", "))
  }
  if (!randomness %in% LINSOLVE_RANDOMNESS) {
    stopf("randomness must be one of: %s", paste(LINSOLVE_RANDOMNESS, collapse = ", "))
  }
  structure(list(apply_inverse = apply_inverse, dim = as.integer(dim),
                 contract = list(fidelity = fidelity, determinacy = determinacy,
                                 randomness = randomness,
                                 adjoint = if (is.null(adjoint)) "absent" else "available",
                                 error_bound = if (error_bound) "available" else "absent"),
                 adjoint = adjoint, method = method),
            class = "linsolve")
}

is_linsolve <- function(x) inherits(x, "linsolve")

## Backends and eigensolvers declare which contracts they accept. Ordinary
## Arnoldi must not receive a history-dependent solve merely because both objects
## implement solve().
accepts_contract <- function(accepted, contract) {
  if ("exact" %in% accepted && contract$fidelity == "exact") return(TRUE)
  if ("fixed_linear" %in% accepted &&
      contract$fidelity == "linear_approximation" &&
      contract$determinacy == "fixed") return(TRUE)
  if ("variable_inexact" %in% accepted &&
      contract$fidelity == "variable_inexact" &&
      contract$error_bound == "available") return(TRUE)
  FALSE
}

check_solve_contract <- function(S, accepted, consumer) {
  if (!is_linsolve(S)) stopf("expected a linsolve object")
  if (!accepts_contract(accepted, S$contract)) {
    stopf(paste0("%s does not accept this solve object.\n",
                 "  declared: fidelity = %s, determinacy = %s, error_bound = %s\n",
                 "  accepted: %s"),
          consumer, S$contract$fidelity, S$contract$determinacy,
          S$contract$error_bound, paste(accepted, collapse = ", "))
  }
  invisible(TRUE)
}

#' @export
print.linsolve <- function(x, ...) {
  cat(sprintf("<linsolve> %d x %d, method '%s'\n", x$dim[1L], x$dim[2L], x$method))
  cat(sprintf("  fidelity:    %s\n", x$contract$fidelity))
  cat(sprintf("  determinacy: %s\n", x$contract$determinacy))
  cat(sprintf("  randomness:  %s\n", x$contract$randomness))
  cat(sprintf("  adjoint:     %s\n", x$contract$adjoint))
  cat(sprintf("  error_bound: %s\n", x$contract$error_bound))
  invisible(x)
}

## A history-dependent solve returns different output for the same input by
## design, so it fails the purity check. linsolve objects get their own suite,
## which tests the declared contract rather than assuming purity, and running the
## linop suite against one is itself an error.

#' @rdname verify
#' @export
verify.linsolve <- function(x, tol = 1e-8, n_probe = 5L, seed = 1L, ...) {
  S <- x
  set.seed(seed)
  n <- S$dim[2L]
  checks <- list()
  add <- function(name, status, detail = "", guarantee = "identity", confidence = 1) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, status = status, source = "computation", guarantee = guarantee,
      confidence = confidence, detail = detail, stringsAsFactors = FALSE)
  }

  b <- matrix(stats::rnorm(n), n, 1)
  z1 <- S$apply_inverse(b); z2 <- S$apply_inverse(b)

  ## Repeatability is required only of a contract that claims it.
  if (S$contract$determinacy == "fixed" && S$contract$randomness == "deterministic") {
    add("repeatability", if (isTRUE(all.equal(z1, z2))) "pass" else "fail",
        "declared fixed and deterministic")
  } else {
    add("repeatability", "not_checked",
        sprintf("determinacy '%s' does not promise it", S$contract$determinacy),
        confidence = NA_real_)
  }

  ## Additivity is required only of a linear_approximation contract.
  if (S$contract$fidelity %in% c("exact", "linear_approximation")) {
    b1 <- matrix(stats::rnorm(n), n, 1); b2 <- matrix(stats::rnorm(n), n, 1)
    g <- max(Mod(S$apply_inverse(b1 + b2) - (S$apply_inverse(b1) + S$apply_inverse(b2))))
    add("additivity", if (g <= tol * max(1, max(Mod(z1)))) "pass" else "fail",
        sprintf("max gap %.3e", g))
  } else {
    add("additivity", "not_checked",
        "variable_inexact solves are not linear, by declaration", confidence = NA_real_)
  }

  add("shape", if (identical(dim(as_block(z1)), c(n, 1L))) "pass" else "fail", "")

  build_certificate(do.call(rbind, checks), subject = "linsolve")
}
