## Section 3's backend protocol, internal. Which engine computes a spectral
## result is chosen here, and the choice is reported in the certificate's
## dispatch, every candidate with its verdict.
##
## The one-package decision keeps this internal: an engine is a value of
## method =, never a second exported verb, and the registry is how eigs() and
## svds() find out which engines exist rather than a surface anyone else
## registers into. The protocol's fields are section 3's, plus the two it left
## implicit, the package an external engine lives in and the function that runs
## it.
##
## Selection is deterministic and is never a benchmark:
##
##   1. method = a name: that backend, or an error giving its reason.
##   2. method = "auto": the eligible backend of highest priority, where
##      options(linop.backend_priority = c(name = p)) overrides the declared
##      priority.
##
## An external engine is eligible only when its package is installed and
## reports the contract version this file was written against. A mismatch is a
## rejection with a reason, not an attempt.

BACKEND_CONTRACT <- "1.1"

.backends <- new.env(parent = emptyenv())

## supports(A, request) returns TRUE or a string saying why not, never FALSE,
## so every exclusion is explainable (section 3). capabilities names what the
## method requires of the operator, and targets the values of which it can
## select; both are checked here for every engine, so a requirement is stated once
## in the registry entry rather than in each verb. run takes the operator (for
## eigs, also the operator the iteration runs on), ctl, the list of k, sel,
## tol, maxit, ncv, v0, seed and floor_const, and the engine; it returns the
## pairs in the shape lanczos_core() and gkl_core() return them.
##
## evidence_minimum is a requirement(), or a function returning one, the evidence for each capability has to
## satisfy before auto acts on it. It filters auto and does not gate a method
## named by the caller, who is then asserting their own declaration: the split
## solve() makes for CG and MINRES, for the same reason.
register_backend <- function(name, family, run, priority, provides_label,
                             engine_version = NA_character_, package = NULL,
                             info_fn = NULL,
                             contract_version = BACKEND_CONTRACT,
                             dtypes = c("double", "complex"),
                             requires = "apply", targets = character(0),
                             capabilities = character(0),
                             evidence_minimum = NULL,
                             default_ncv = function(A, k) NULL,
                             supports = function(A, request) TRUE) {
  key <- paste(family, name, sep = ":")
  assign(key, list(name = name, family = family, run = run,
                   priority = as.numeric(priority), label = provides_label,
                   engine_version = engine_version, package = package,
                   info_fn = info_fn,
                   contract_version = contract_version, dtypes = dtypes,
                   requires = requires, targets = targets,
                   capabilities = capabilities,
                   evidence_minimum = evidence_minimum,
                   default_ncv = default_ncv, supports = supports),
         envir = .backends)
  invisible(key)
}

backends_for <- function(family) {
  keys <- ls(.backends, pattern = paste0("^", family, ":"))
  lapply(keys, get, envir = .backends)
}

## The exported functions of an external engine, or a string saying why there
## are none. The package name is data in the registry entry. An engine package
## written for linop reports the contract it speaks through info_fn; a library
## written for its own users has no such function, and its version is the record.
backend_engine <- function(b) {
  if (is.null(b$package)) return(NULL)
  if (!requireNamespace(b$package, quietly = TRUE)) {
    return(sprintf("package %s is not installed", b$package))
  }
  if (is.null(b$info_fn)) {
    ns <- asNamespace(b$package)
    return(list(info = list(version = as.character(utils::packageVersion(b$package))),
                fn = function(name) getExportedValue(ns, name)))
  }
  info <- getExportedValue(b$package, b$info_fn)()
  if (!identical(info$contract_version, b$contract_version)) {
    return(sprintf("%s speaks contract %s; this linop speaks %s", b$package,
                   info$contract_version, b$contract_version))
  }
  ns <- asNamespace(b$package)
  list(info = info, fn = function(name) getExportedValue(ns, name))
}

backend_priority <- function(b) {
  over <- getOption("linop.backend_priority")
  if (!is.null(over) && b$name %in% names(over)) as.numeric(over[[b$name]])
  else b$priority
}

## The verdict on every backend of a family, in the order auto considers them.
backend_verdicts <- function(family, A, request) {
  bs <- backends_for(family)
  pr <- vapply(bs, backend_priority, numeric(1))
  bs <- bs[order(-pr)]
  lapply(bs, function(b) {
    why <- backend_refusal(b, A, request)
    list(backend = b, priority = backend_priority(b), eligible = is.null(why),
         reason = why %||% "eligible")
  })
}

## Why a backend cannot take this request, or NULL. The checks run cheapest
## first, and the one that loads an external package runs last.
backend_refusal <- function(b, A, request) {
  if (!A$dtype %in% b$dtypes) return(sprintf("does not take a %s operator", A$dtype))
  for (cap in b$capabilities) {
    why <- capability_refusal(A, cap)
    if (!is.null(why)) return(why)
    if (isTRUE(request$auto)) {
      why <- evidence_refusal(A, cap, b$evidence_minimum)
      if (!is.null(why)) return(why)
    }
  }
  if (length(b$targets) && !is.null(request$which) && !request$which %in% b$targets) {
    return(sprintf("does not select which = '%s'", request$which))
  }
  eng <- backend_engine(b)
  if (is.character(eng)) return(eng)
  s <- b$supports(A, request)
  if (isTRUE(s)) NULL else s
}

select_backend <- function(family, A, request, method, verb) {
  request$auto <- identical(method, "auto")
  verdicts <- backend_verdicts(family, A, request)
  names_all <- vapply(verdicts, function(v) v$backend$name, character(1))
  if (!is_scalar_string(method) || !method %in% c("auto", names_all)) {
    stopf("method is one of %s",
          paste0("'", c("auto", sort(names_all)), "'", collapse = ", "))
  }
  table <- data.frame(
    backend = names_all,
    priority = vapply(verdicts, function(v) v$priority, numeric(1)),
    status = vapply(verdicts, function(v) if (v$eligible) "eligible" else "rejected",
                    character(1)),
    reason = vapply(verdicts, function(v) v$reason, character(1)),
    stringsAsFactors = FALSE)
  if (identical(method, "auto")) {
    ok <- which(table$status == "eligible")
    if (!length(ok)) {
      stopf("%s() has no eligible engine:\n%s", verb,
            paste0("  ", table$backend, ": ", table$reason, collapse = "\n"))
    }
    pick <- ok[1L]
    why <- if (pick == 1L) "highest priority" else
      sprintf("highest priority eligible; %s", paste0(
        table$backend[seq_len(pick - 1L)], " ", table$reason[seq_len(pick - 1L)],
        collapse = "; "))
  } else {
    pick <- match(method, table$backend)
    if (table$status[pick] != "eligible") {
      stopf("%s(method = '%s') cannot run: %s", verb, method, table$reason[pick])
    }
    why <- "requested by name"
  }
  table$status[pick] <- "selected"
  b <- verdicts[[pick]]$backend
  list(backend = b, engine = backend_engine(b), table = table, why = why)
}

## The dispatch record a spectral result carries.
backend_dispatch <- function(method, sel, label) {
  list(requested = method, chosen = sel$backend$name,
       reason = sprintf("%s (%s)", label, sel$why),
       candidates = sel$table)
}

## ------------------------------------------------------------ the registry --

## The hermitian engines refuse a complex shift: A - sigma I is then not
## hermitian, and Arnoldi takes it instead.
hermitian_shift_supports <- function(request) {
  if (is.complex(request$sigma)) {
    "a complex shift leaves A - sigma I non-hermitian"
  } else TRUE
}

## A pencil reaches the hermitian recurrence only in the B inner product, which
## is an inner product only for B hermitian positive definite.
## Under auto the evidence for both has to meet the minimum CG's does, since the
## inner product is exactly as good as B's definiteness.
definite_pencil_supports <- function(request) {
  B <- request$B
  if (is.null(B)) return(TRUE)
  why <- capability_refusal(B, "hermitian") %||% capability_refusal(B, "positive_definite")
  if (is.null(why) && isTRUE(request$auto)) {
    why <- evidence_refusal(B, "hermitian", CG_PD_REQUIREMENT) %||%
           evidence_refusal(B, "positive_definite", CG_PD_REQUIREMENT)
  }
  if (is.null(why)) TRUE
  else paste0("the B inner product ", sub("the operator", "B", why, fixed = TRUE))
}

## Why auto will not act on a capability the operator does hold: its evidence is
## below the minimum. NULL when it meets it or there is no minimum.
evidence_refusal <- function(A, name, req) {
  if (is.function(req)) req <- req()
  if (is.null(req) || evidence_satisfies(cape(A, name), req)) return(NULL)
  sprintf(paste0("auto acts on %s only at computed, constructed or proven evidence, ",
                 "and the operator's is %s; name the method to rely on it"),
          name, cape(A, name)$source)
}

## The minimum auto applies to a hermitian engine: MINRES's, the same reading of
## the evidence fields solve() uses to choose a hermitian method. A function,
## because this file is collated before the one defining it.
eigs_hermitian_requirement <- function() MINRES_HERMITIAN_REQUIREMENT

register_backend("lanczos", "eigs", priority = 0,
                 provides_label = "Lanczos with full reorthogonalisation",
                 targets = c("LM", "SM", "LA", "SA"), capabilities = "hermitian",
                 evidence_minimum = eigs_hermitian_requirement,
                 default_ncv = function(A, k) default_ncv(A, k),
                 supports = function(A, request) {
                   s <- hermitian_shift_supports(request)
                   if (isTRUE(s)) definite_pencil_supports(request) else s
                 },
                 run = function(Op, A, ctl, engine) {
                   run <- lanczos_core(Op, A, ctl$k, ctl$sel, ctl$tol, ctl$maxit, ctl$ncv,
                                       ctl$v0, ctl$seed, ctl$floor_const, B = ctl$B)
                   if (!is.null(ctl$B)) return(run)
                   deflation_sweep(function(k, maxit, seed, prelocked)
                     lanczos_core(Op, A, k, ctl$sel, ctl$tol, maxit, ctl$ncv, NULL, seed,
                                  ctl$floor_const, prelocked = prelocked),
                     run, ctl$k, ctl$sel, ctl$tol, ctl$maxit, ctl$seed)
                 })

## Below the reference Lanczos, so that a hermitian operator keeps the method
## whose projected problem is real and symmetric; above nothing, because nothing
## else runs on an operator that establishes no capability.
register_backend("arnoldi", "eigs", priority = -10,
                 provides_label = "Arnoldi with Krylov-Schur restarting",
                 targets = c("LM", "SM", "LA", "SA", "LI", "SI"),
                 default_ncv = function(A, k) default_ncv(A, k),
                 run = function(Op, A, ctl, engine) {
                   herm <- isTRUE(capv(A, "hermitian")) &&
                     isTRUE(definite_pencil_supports(list(B = ctl$B)))
                   run <- arnoldi_core(Op, A, ctl$k, ctl$sel, ctl$tol, ctl$maxit, ctl$ncv,
                                       ctl$v0, ctl$seed, ctl$floor_const,
                                       hermitian = herm, B = ctl$B)
                   ## The sweep for missed copies, where copies exist as separate
                   ## eigenvectors: a hermitian operator. A repeated eigenvalue of
                   ## a non-normal one is generically defective, and a candidate
                   ## there could split a conjugate pair.
                   if (!herm || !is.null(ctl$B)) return(run)
                   deflation_sweep(function(k, maxit, seed, prelocked)
                     arnoldi_core(Op, A, k, ctl$sel, ctl$tol, maxit, ctl$ncv, NULL, seed,
                                  ctl$floor_const, hermitian = herm, prelocked = prelocked),
                     run, ctl$k, ctl$sel, ctl$tol, ctl$maxit, ctl$seed)
                 })

register_backend("golub-kahan", "svds", priority = 0,
                 provides_label = "Golub-Kahan-Lanczos with full reorthogonalisation",
                 requires = c("apply", "adjoint"), targets = c("LM", "SM"),
                 default_ncv = function(A, k) default_ncv(A, k),
                 run = function(A, ctl, engine) gkl_svds(A, ctl))

register_backend("primme", "eigs", priority = 100,
                 provides_label = "PRIMME with implicit restarting",
                 engine_version = "3.2.3", package = "linop.primme",
                 info_fn = "primme_engine", targets = c("LM", "SM", "LA", "SA"),
                 capabilities = "hermitian", evidence_minimum = eigs_hermitian_requirement,
                 supports = function(A, request) {
                   if (!is.null(request$B)) return("it does not take B")
                   s <- hermitian_shift_supports(request)
                   if (isTRUE(s)) primme_supports(request) else s
                 },
                 run = function(Op, A, ctl, engine) primme_eigs_run(Op, A, ctl, engine))

register_backend("primme", "svds", priority = 100,
                 provides_label = "PRIMME with implicit restarting",
                 engine_version = "3.2.3", package = "linop.primme",
                 info_fn = "primme_engine", requires = c("apply", "adjoint"),
                 targets = c("LM", "SM"),
                 supports = function(A, request) {
                   if (isTRUE(request$v0)) {
                     "it starts from a left and a right vector, and v0 is only the right one"
                   } else primme_supports(request)
                 },
                 run = function(A, ctl, engine) primme_svds_run(A, ctl, engine))

## Reachable by name only: below Arnoldi, which takes every operator it takes.
register_backend("rspectra", "eigs", priority = -20,
                 provides_label = "RSpectra (Spectra's implicitly restarted Arnoldi and Lanczos)",
                 package = "RSpectra", dtypes = "double",
                 targets = c("LM", "SM", "LA", "SA", "LI", "SI"),
                 default_ncv = function(A, k) default_ncv(A, k),
                 supports = rspectra_eigs_supports,
                 run = function(Op, A, ctl, engine) rspectra_eigs_run(Op, A, ctl, engine))

register_backend("rspectra", "svds", priority = -20,
                 provides_label = "RSpectra (Spectra's implicitly restarted Lanczos)",
                 package = "RSpectra", dtypes = "double",
                 requires = c("apply", "adjoint"), targets = "LM",
                 default_ncv = function(A, k) default_ncv(A, k),
                 supports = rspectra_svds_supports,
                 run = function(A, ctl, engine) rspectra_svds_run(A, ctl, engine))
