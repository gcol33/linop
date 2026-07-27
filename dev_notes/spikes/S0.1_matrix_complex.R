suppressMessages(library(Matrix))
cat("Matrix version:", as.character(packageVersion("Matrix")), "\n\n")
Z <- matrix(complex(real = c(1, 2, 3, 4), imaginary = c(1, -1, 2, 0.5)), 2, 2)

cat("Matrix(Z) ->", tryCatch({ m <- Matrix(Z); paste(class(m), collapse=",") },
                             error = function(e) paste("ERR:", conditionMessage(e))), "\n")
cat("as(Z,'dgCMatrix') ->", tryCatch(paste(class(as(Z, "dgCMatrix")), collapse=","),
                             error = function(e) paste("ERR:", substr(conditionMessage(e),1,60))), "\n")
for (cl in c("zMatrix", "zgCMatrix", "zgeMatrix", "zsparseMatrix")) {
  cat(sprintf("  isVirtualClass/exists %-14s: %s\n", cl,
      tryCatch(methods::existsMethod("show", cl) || methods::isVirtualClass(cl) ||
               !is.null(methods::getClassDef(cl, package = "Matrix")),
               error = function(e) "no")))
}
cat("\nsparseMatrix() with complex x ->",
    tryCatch(paste(class(sparseMatrix(i = c(1,2), j = c(1,2), x = c(1+1i, 2-1i))), collapse = ","),
             error = function(e) paste("ERR:", substr(conditionMessage(e), 1, 70))), "\n")

cat("\n--- real sparse works, for the adapter ---\n")
S <- sparseMatrix(i = c(1,2,2), j = c(1,1,2), x = c(1, 2, 3))
cat("class:", paste(class(S), collapse=","), "| crossprod(S) class:", paste(class(crossprod(S)), collapse=","), "\n")
cat("dsCMatrix (symmetric) recognised:", is(crossprod(S), "symmetricMatrix"), "\n")

cat("\n--- which Matrix classes will the adapter meet? ---\n")
for (nm in c("dgCMatrix","dsCMatrix","dtCMatrix","dgeMatrix","dsyMatrix","ddiMatrix","dgRMatrix","dgTMatrix"))
  cat(sprintf("  %-11s exists: %s\n", nm, !is.null(methods::getClassDef(nm, package = "Matrix"))))
