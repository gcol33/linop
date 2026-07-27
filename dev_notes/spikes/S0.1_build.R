## S0.1 -- dispatch spike. Builds three minimal packages (S3, S4, S7) each
## defining a matrix-like class, installs them to a temp library, and probes
## dispatch for every generic in the linop public budget.

root <- "C:/Users/GILLES~1/AppData/Local/Temp/claude/C--Users-Gilles-Colling-Documents-dev-linop/60941760-700a-476f-a7e6-656dcb3a2335/scratchpad/s01"
unlink(root, recursive = TRUE)
dir.create(root, recursive = TRUE, showWarnings = FALSE)
lib <- file.path(root, "lib")
dir.create(lib)

write_pkg <- function(name, desc_extra, namespace, code) {
  d <- file.path(root, name)
  dir.create(file.path(d, "R"), recursive = TRUE)
  writeLines(c(
    paste0("Package: ", name),
    "Version: 0.0.1",
    paste0("Title: Dispatch probe ", name),
    "Author: spike",
    "Maintainer: spike <spike@example.com>",
    "Description: Minimal dispatch probe.",
    "License: MIT + file LICENSE",
    desc_extra
  ), file.path(d, "DESCRIPTION"))
  writeLines(c("YEAR: 2026", "COPYRIGHT HOLDER: spike"), file.path(d, "LICENSE"))
  writeLines(namespace, file.path(d, "NAMESPACE"))
  writeLines(code, file.path(d, "R", "cls.R"))
  d
}

## ---------------------------------------------------------------- S3 -------
## Two variants in one package: `opA` uses individual methods, `opG` uses the
## matrixOps *group* generic, to learn whether one group method covers all three
## matrix-multiply generics.
s3 <- write_pkg(
  "dispS3", "Depends: R (>= 4.4.0)",
  c('export(opA, opG)',
    'S3method("%*%", opA)', 'S3method(crossprod, opA)', 'S3method(tcrossprod, opA)',
    'S3method(t, opA)', 'S3method(Conj, opA)', 'S3method(solve, opA)',
    'S3method(dim, opA)', 'S3method(Ops, opA)', 'S3method(as.matrix, opA)',
    'S3method(print, opA)',
    'S3method(matrixOps, opG)', 'S3method(dim, opG)', 'S3method(Complex, opG)'),
  c(
    'opA <- function(m) structure(list(m = m), class = "opA")',
    'opG <- function(m) structure(list(m = m), class = "opG")',
    'un <- function(x) if (inherits(x, c("opA","opG"))) x$m else x',
    '',
    'dim.opA <- function(x) dim(x$m)',
    't.opA <- function(x) opA(t(x$m))',
    'Conj.opA <- function(z) opA(Conj(z$m))',
    'solve.opA <- function(a, b, ...) if (missing(b)) opA(solve(a$m)) else solve(a$m, un(b))',
    'as.matrix.opA <- function(x, ...) x$m',
    'print.opA <- function(x, ...) { cat("<opA>\\n"); invisible(x) }',
    '"%*%.opA" <- function(x, y) structure(list(tag = "opA-matmul", v = un(x) %*% un(y)), class = "probe")',
    'crossprod.opA <- function(x, y = NULL) structure(list(tag = "opA-crossprod", v = if (is.null(y)) crossprod(un(x)) else crossprod(un(x), un(y))), class = "probe")',
    'tcrossprod.opA <- function(x, y = NULL) structure(list(tag = "opA-tcrossprod", v = if (is.null(y)) tcrossprod(un(x)) else tcrossprod(un(x), un(y))), class = "probe")',
    'Ops.opA <- function(e1, e2) structure(list(tag = paste0("opA-Ops-", .Generic)), class = "probe")',
    '',
    'dim.opG <- function(x) dim(x$m)',
    'matrixOps.opG <- function(x, y) structure(list(tag = paste0("opG-matrixOps-", .Generic)), class = "probe")',
    'Complex.opG <- function(z) structure(list(tag = paste0("opG-Complex-", .Generic)), class = "probe")'
  ))

## ---------------------------------------------------------------- S4 -------
s4 <- write_pkg(
  "dispS4", "Depends: R (>= 4.4.0)\nImports: methods",
  c('export(op4)', 'exportClasses(op4cls)', 'exportMethods("%*%", crossprod, tcrossprod, t, Conj, solve, dim)',
    'importFrom(methods, setClass, setMethod, setGeneric, new, representation)'),
  c(
    'setClass("op4cls", representation(m = "matrix"))',
    'op4 <- function(m) new("op4cls", m = m)',
    'setMethod("dim", "op4cls", function(x) dim(x@m))',
    'setMethod("t", "op4cls", function(x) op4(t(x@m)))',
    'setMethod("Conj", "op4cls", function(z) op4(Conj(z@m)))',
    'setMethod("solve", signature("op4cls", "ANY"), function(a, b, ...) structure(list(tag = "S4-solve"), class = "probe"))',
    'setMethod("%*%", signature("op4cls", "ANY"), function(x, y) structure(list(tag = "S4-matmul"), class = "probe"))',
    'setMethod("%*%", signature("ANY", "op4cls"), function(x, y) structure(list(tag = "S4-matmul-rhs"), class = "probe"))',
    'setMethod("crossprod", signature("op4cls", "ANY"), function(x, y = NULL) structure(list(tag = "S4-crossprod"), class = "probe"))',
    'setMethod("tcrossprod", signature("op4cls", "ANY"), function(x, y = NULL) structure(list(tag = "S4-tcrossprod"), class = "probe"))'
  ))

## ---------------------------------------------------------------- S7 -------
s7 <- write_pkg(
  "dispS7", "Depends: R (>= 4.4.0)\nImports: S7",
  c('export(op7)', 'import(S7)'),
  c(
    'op7 <- S7::new_class("op7", properties = list(m = S7::class_any))',
    'S7::method(`%*%`, list(op7, S7::class_any)) <- function(x, y) structure(list(tag = "S7-matmul"), class = "probe")',
    'S7::method(`%*%`, list(S7::class_any, op7)) <- function(x, y) structure(list(tag = "S7-matmul-rhs"), class = "probe")',
    'S7::method(crossprod, list(op7, S7::class_any)) <- function(x, y = NULL) structure(list(tag = "S7-crossprod"), class = "probe")',
    'S7::method(tcrossprod, list(op7, S7::class_any)) <- function(x, y = NULL) structure(list(tag = "S7-tcrossprod"), class = "probe")',
    'S7::method(t, op7) <- function(x) structure(list(tag = "S7-t"), class = "probe")',
    'S7::method(Conj, op7) <- function(z) structure(list(tag = "S7-Conj"), class = "probe")',
    'S7::method(solve, list(op7, S7::class_any)) <- function(a, b, ...) structure(list(tag = "S7-solve"), class = "probe")',
    'S7::method(dim, op7) <- function(x) dim(S7::prop(x, "m"))',
    '.onLoad <- function(libname, pkgname) S7::methods_register()'
  ))

for (p in c(s3, s4, s7)) {
  cat("\n#### installing", basename(p), "\n")
  out <- tryCatch(
    utils::install.packages(p, lib = lib, repos = NULL, type = "source",
                            INSTALL_opts = "--no-docs --no-byte-compile"),
    error = function(e) cat("INSTALL ERROR:", conditionMessage(e), "\n"))
}

cat("\n#### installed:\n")
print(list.files(lib))
cat("\nLIB=", lib, "\n", sep = "")
writeLines(lib, file.path(root, "libpath.txt"))
