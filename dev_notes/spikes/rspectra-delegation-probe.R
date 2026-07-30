## What RSpectra can carry as a delegation target, measured rather than read off its
## signature. Plan section 7.2 lists RSpectra delegation for v0.1 and
## dev_notes/eigs-svds-and-the-third-certificate.md defers it to Phase 3 with the
## backend registry, describing its interface as "single vector and real". This probe
## establishes what happens at each of those two limits when they are crossed, since a
## delegating wrapper has to refuse what it cannot pass on and an interface that warns
## and then answers is the case that needs a number.
##
## Run from the repository root:
##   "C:/Program Files/R/R-4.6.0/bin/x64/Rscript.exe" dev_notes/spikes/rspectra-delegation-probe.R
##
## Writes results/rspectra-delegation.csv and results/rspectra-delegation.txt.

if (!requireNamespace("RSpectra", quietly = TRUE)) {
  stop("RSpectra is not installed; nothing to probe")
}

out_dir <- file.path("dev_notes", "spikes", "results")
if (!dir.exists(out_dir)) {
  stop("run from the repository root; dev_notes/spikes/results not found")
}

txt <- file(file.path(out_dir, "rspectra-delegation.txt"), open = "wt")
say <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = txt)
}

rows <- list()
record <- function(case, route, dtype, outcome, detail) {
  rows[[length(rows) + 1L]] <<- data.frame(
    case = case, route = route, dtype = dtype,
    outcome = outcome, detail = detail, stringsAsFactors = FALSE)
}

say("RSpectra ", as.character(packageVersion("RSpectra")),
    " | R ", paste(R.version$major, R.version$minor, sep = "."),
    " | ", R.version$platform)
say("")

n <- 60
set.seed(1)

## ---------------------------------------------------------------- 1. real, non-symmetric
## The gap eigs() names at R/eigs.R:111. RSpectra is ARPACK, so this is the route that
## would close it without any compiled code entering linop.
Ans <- matrix(rnorm(n * n), n, n)
got <- RSpectra::eigs(function(x, args) as.numeric(args %*% x),
                      k = 3, n = n, args = Ans)$values
ref <- eigen(Ans, only.values = TRUE)$values
ref <- ref[order(-Mod(ref))][seq_len(3)]
d1 <- max(Mod(sort(Mod(got)) - sort(Mod(ref))))
say("1. real non-symmetric, matrix-free")
say("   max |diff| vs dense eigen(): ", format(d1, digits = 4))
record("non-symmetric", "matrix-free", "double", "works",
       sprintf("max |diff| vs dense = %.3e", d1))

## ---------------------------------------------------------------- 2. block apply
## linop's tier-1 apply is block. If RSpectra ever asked for more than one column the
## delegation could keep that; this records whether it does.
shapes <- character(0)
Asym <- crossprod(Ans)
invisible(RSpectra::eigs_sym(
  function(x, args) {
    shapes <<- c(shapes, paste(dim(as.matrix(x)), collapse = "x"))
    as.numeric(args %*% x)
  }, k = 3, n = n, args = Asym))
say("")
say("2. callback input shapes over a whole run: ", paste(unique(shapes), collapse = ", "),
    "  (", length(shapes), " calls)")
record("block apply", "matrix-free", "double",
       if (identical(unique(shapes), paste0(n, "x1"))) "single column only" else "block seen",
       paste0("shapes: ", paste(unique(shapes), collapse = ", ")))

## ---------------------------------------------------------------- 3. complex, dense
set.seed(7)
m <- 40
Ac <- matrix(complex(real = rnorm(m * m), imaginary = rnorm(m * m)), m, m)
Ac <- Ac + t(Conj(Ac))
stopifnot(max(Mod(Ac - t(Conj(Ac)))) == 0)

dense_err <- tryCatch({ RSpectra::eigs(Ac, k = 3); NA_character_ },
                      error = function(e) conditionMessage(e))
say("")
say("3. complex dense: ", if (is.na(dense_err)) "accepted" else paste0("ERROR: ", dense_err))
record("complex", "dense", "complex",
       if (is.na(dense_err)) "accepted" else "refused",
       if (is.na(dense_err)) "" else dense_err)

## ---------------------------------------------------------------- 4. complex, matrix-free
## The case that matters. The callback returns a complex vector and RSpectra coerces it.
nwarn <- 0L
vals <- withCallingHandlers(
  RSpectra::eigs(function(x, args) as.vector(args %*% x), k = 3, n = m, args = Ac)$values,
  warning = function(w) { nwarn <<- nwarn + 1L; invokeRestart("muffleWarning") })

truth <- eigen(Ac, only.values = TRUE)$values
truth <- Re(truth[order(-Mod(truth))][seq_len(3)])
re_only <- eigen(Re(Ac), only.values = TRUE)$values
re_only <- Re(re_only[order(-Mod(re_only))][seq_len(3)])

d_truth <- max(abs(sort(Re(vals)) - sort(truth)))
d_re    <- max(abs(sort(Re(vals)) - sort(re_only)))

say("")
say("4. complex hermitian, matrix-free")
say("   returned            : ", paste(format(Re(vals), digits = 8), collapse = "  "))
say("   eigenvalues of A    : ", paste(format(truth, digits = 8), collapse = "  "))
say("   eigenvalues of Re(A): ", paste(format(re_only, digits = 8), collapse = "  "))
say("   warnings raised     : ", nwarn)
say("   max |diff| vs A     : ", format(d_truth, digits = 4))
say("   max |diff| vs Re(A) : ", format(d_re, digits = 4))
record("complex", "matrix-free", "complex", "answers the wrong operator",
       sprintf("vs A = %.3e, vs Re(A) = %.3e, warnings = %d", d_truth, d_re, nwarn))

## ---------------------------------------------------------------- 5. svds, matrix-free
## svds needs both actions. linop supplies them as modes N and C of one operator, so the
## question is whether RSpectra takes the second as a separate argument.
say("")
say("5. svds formals: ", paste(names(formals(RSpectra::svds)), collapse = ", "))
B <- matrix(rnorm(n * (n - 20)), n, n - 20)
svd_res <- tryCatch(
  RSpectra::svds(function(x, args) as.numeric(args %*% x),
                 k = 3, nu = 3, nv = 3, dim = dim(B), args = B,
                 Atrans = function(x, args) as.numeric(crossprod(args, x))),
  error = function(e) conditionMessage(e))
if (is.character(svd_res)) {
  say("   matrix-free svds: ERROR: ", svd_res)
  record("svds", "matrix-free", "double", "refused", svd_res)
} else {
  d5 <- max(abs(sort(svd_res$d) - sort(svd(B)$d[seq_len(3)])))
  say("   matrix-free svds: works, max |diff| vs svd() = ", format(d5, digits = 4))
  record("svds", "matrix-free", "double", "works",
         sprintf("max |diff| vs svd() = %.3e", d5))
}

## ---------------------------------------------------------------- 6. targets and sigma
say("")
say("6. eigs formals:     ", paste(names(formals(RSpectra::eigs)), collapse = ", "))
say("   eigs_sym formals: ", paste(names(formals(RSpectra::eigs_sym)), collapse = ", "))
sm <- tryCatch({
  v <- RSpectra::eigs_sym(function(x, args) as.numeric(args %*% x),
                          k = 3, which = "SM", n = n, args = Asym)$values
  sprintf("SM works, smallest = %.6g", min(v))
}, error = function(e) paste0("ERROR: ", conditionMessage(e)))
say("   which = 'SM', matrix-free: ", sm)
record("smallest", "matrix-free", "double",
       if (grepl("^ERROR", sm)) "refused" else "works", sm)

sig_dense <- tryCatch({
  v <- RSpectra::eigs_sym(Asym, k = 3, sigma = 0)$values
  sprintf("works, smallest = %.6g", min(v))
}, error = function(e) paste0("ERROR: ", conditionMessage(e)))
sig_fn <- tryCatch({
  RSpectra::eigs_sym(function(x, args) as.numeric(args %*% x),
                     k = 3, sigma = 0, n = n, args = Asym)
  "accepted"
}, error = function(e) paste0("ERROR: ", conditionMessage(e)))
say("   sigma on a function: ", sig_fn)
say("   sigma on a matrix:   ", sig_dense)
record("shift-invert", "matrix-free", "double",
       if (identical(sig_fn, "accepted")) "accepted" else "refused",
       paste0("function: ", sig_fn, " | dense: ", sig_dense))

## ----------------------------------------------------------------
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "rspectra-delegation.csv"), row.names = FALSE)
say("")
say("wrote ", file.path(out_dir, "rspectra-delegation.csv"))
close(txt)
