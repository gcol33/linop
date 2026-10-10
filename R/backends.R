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
## so every exclusion is explainable (section 3). run takes the operator (for
## eigs, also the operator the iteration runs on), ctl, the list of k, sel,
## tol, maxit, ncv, v0, seed and floor_const, and the engine; it returns the
## pairs in the shape lanczos_core() and gkl_core() return them.
register_backend <- function(name, family, run, priority, provides_label,
                             engine_version = NA_character_, package = NULL,
                             info_fn = NULL,
                             contract_version = BACKEND_CONTRACT,
                             dtypes = c("double", "complex"),
                             requires = "apply", targets = character(0),
                             default_ncv = function(A, k) NULL,
                             supports = function(A, request) TRUE) {
  key <- paste(family, name, sep = ":")
  assign(key, list(name = name, family = family, run = run,
                   priority = as.numeric(priority), label = provides_label,
                   engine_version = engine_version, package = package,
                   info_fn = info_fn,
                   contract_version = contract_version, dtypes = dtypes,
                   requires = requires, targets = targets,
                   default_ncv = default_ncv, supports = supports),
         envir = .backends)
  invisible(key)
}

backends_for <- function(family) {
  keys <- ls(.backends, pattern = paste0("^", family, ":"))
  lapply(keys, get, envir = .backends)
}

## The exported functions of an external engine, or a string saying why there
## are none. The package name is data in the registry entry.
backend_engine <- function(b) {
  if (is.null(b$package)) return(NULL)
  if (!requireNamespace(b$package, quietly = TRUE)) {
    return(sprintf("package %s is not installed", b$package))
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
    why <- if (!A$dtype %in% b$dtypes) {
      sprintf("does not take a %s operator", A$dtype)
    } else {
      eng <- backend_engine(b)
      if (is.character(eng)) eng
      else {
        s <- b$supports(A, request)
        if (isTRUE(s)) NULL else s
      }
    }
    list(backend = b, priority = backend_priority(b), eligible = is.null(why),
         reason = why %||% "eligible")
  })
}

select_backend <- function(family, A, request, method, verb) {
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

register_backend("lanczos", "eigs_hermitian", priority = 0,
                 provides_label = "Lanczos with full reorthogonalisation",
                 targets = c("LM", "SM", "LA", "SA"),
                 default_ncv = function(A, k) default_ncv(A, k),
                 run = function(Op, A, ctl, engine)
                   lanczos_core(Op, A, ctl$k, ctl$sel, ctl$tol, ctl$maxit, ctl$ncv,
                                ctl$v0, ctl$seed, ctl$floor_const))

register_backend("golub-kahan", "svds", priority = 0,
                 provides_label = "Golub-Kahan-Lanczos with full reorthogonalisation",
                 requires = c("apply", "adjoint"), targets = c("LM", "SM"),
                 default_ncv = function(A, k) default_ncv(A, k),
                 run = function(A, ctl, engine)
                   gkl_core(A, ctl$k, ctl$sel, ctl$tol, ctl$maxit, ctl$ncv,
                            ctl$v0, ctl$seed, ctl$floor_const))

register_backend("primme", "eigs_hermitian", priority = 100,
                 provides_label = "PRIMME with implicit restarting",
                 engine_version = "3.2.3", package = "linop.primme",
                 info_fn = "primme_engine", targets = c("LM", "SM", "LA", "SA"),
                 supports = function(A, request) primme_supports(request),
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
