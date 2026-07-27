lib <- "C:/Users/GILLES~1/AppData/Local/Temp/claude/C--Users-Gilles-Colling-Documents-dev-linop/60941760-700a-476f-a7e6-656dcb3a2335/scratchpad/s01/lib"
.libPaths(c(lib, .libPaths()))
library(dispS3); library(dispS4)

tag <- function(e) {
  v <- tryCatch(eval(e, envir = globalenv()), error = function(err) paste("ERR:", substr(conditionMessage(err), 1, 55)))
  if (inherits(v, "probe")) v$tag
  else if (is.character(v)) v
  else if (inherits(v, c("opA", "opG", "op4cls"))) paste0("<", class(v)[1], ">")
  else paste(class(v)[1], ":", paste(utils::head(as.vector(v), 4), collapse = ","))
}

M <- matrix(c(1, 2, 3, 4), 2); Z <- matrix(complex(real = 1:4, imaginary = 4:1), 2)
A <- opA(M); G <- opG(M); S <- op4(M)

cat("################ S3, individual methods (class opA) ################\n")
for (e in c("A %*% M", "M %*% A", "A %*% A", "crossprod(A)", "crossprod(A, M)",
            "tcrossprod(A)", "t(A)", "Conj(A)", "solve(A)", "solve(A, c(1,2))",
            "dim(A)", "nrow(A)", "ncol(A)", "A + A", "2 * A", "as.matrix(A)"))
  cat(sprintf("  %-18s -> %s\n", e, tag(parse(text = e)[[1]])))

cat("\n################ S3, matrixOps GROUP method (class opG) ################\n")
for (e in c("G %*% M", "M %*% G", "crossprod(G)", "crossprod(G, M)", "tcrossprod(G)",
            "dim(G)", "nrow(G)", "ncol(G)", "Conj(G)"))
  cat(sprintf("  %-18s -> %s\n", e, tag(parse(text = e)[[1]])))

cat("\n################ S4 (control, class op4cls) ################\n")
for (e in c("S %*% M", "M %*% S", "crossprod(S)", "tcrossprod(S)", "t(S)",
            "Conj(S)", "solve(S)", "dim(S)", "nrow(S)", "ncol(S)"))
  cat(sprintf("  %-18s -> %s\n", e, tag(parse(text = e)[[1]])))

cat("\n################ are nrow/ncol generic at all? ################\n")
cat("  nrow is primitive/generic? ", is.primitive(base::nrow), " body: ",
    paste(deparse(body(base::nrow)), collapse = " "), "\n", sep = "")
cat("  ncol body: ", paste(deparse(body(base::ncol)), collapse = " "), "\n", sep = "")
cat("  -> nrow/ncol need NO method; they derive from dim(). Confirmed above.\n")

cat("\n################ matrixOps group membership, this R ################\n")
cat("  R version:", as.character(getRversion()), "\n")
cat("  methods::getGroupMembers('matrixOps') -> ",
    paste(tryCatch(methods::getGroupMembers("matrixOps"), error = function(e) "ERR"), collapse = ", "), "\n")
cat("  .S3PrimitiveGenerics contains crossprod? ", "crossprod" %in% .S3PrimitiveGenerics, "\n")
cat("  .S3PrimitiveGenerics contains %*%?       ", "%*%" %in% .S3PrimitiveGenerics, "\n")

cat("\n################ complex: t() vs Conj(t()) must differ ################\n")
cat("  identical(t(Z), Conj(t(Z))) =", identical(t(Z), Conj(t(Z))), " (must be FALSE)\n")
cat("  crossprod(Z) hermitian? ", isTRUE(all.equal(crossprod(Z), Conj(t(crossprod(Z))))), "\n")
cat("  t(Z)%*%Z symmetric non-herm? ",
    isTRUE(all.equal(t(Z) %*% Z, t(t(Z) %*% Z))) && !isTRUE(all.equal(t(Z) %*% Z, Conj(t(t(Z) %*% Z)))), "\n")
