Z <- matrix(complex(real = c(1, 2, 3, 4), imaginary = c(1, -1, 2, 0.5)), 2, 2)
cat("Z =\n"); print(Z)

cp   <- crossprod(Z)
tZZ  <- t(Z) %*% Z
hZZ  <- Conj(t(Z)) %*% Z

cat("\ncrossprod(Z)      =\n"); print(cp)
cat("\nt(Z) %*% Z  (A^T A)=\n"); print(tZZ)
cat("\nConj(t(Z))%*%Z (A^H A)=\n"); print(hZZ)

cat("\n### which does base crossprod equal? ###\n")
cat("  crossprod(Z) == t(Z)%*%Z      (A^T A):", isTRUE(all.equal(cp, tZZ)), "\n")
cat("  crossprod(Z) == Conj(t(Z))%*%Z(A^H A):", isTRUE(all.equal(cp, hZZ)), "\n")

herm <- function(M) isTRUE(all.equal(M, Conj(t(M))))
symm <- function(M) isTRUE(all.equal(M, t(M)))
cat("\n### properties ###\n")
cat("  crossprod(Z)  hermitian:", herm(cp),  " symmetric:", symm(cp),  "\n")
cat("  A^T A         hermitian:", herm(tZZ), " symmetric:", symm(tZZ), "\n")
cat("  A^H A         hermitian:", herm(hZZ), " symmetric:", symm(hZZ), "\n")

cat("\n### tcrossprod ###\n")
tcp <- tcrossprod(Z)
cat("  tcrossprod(Z) == Z%*%t(Z)      :", isTRUE(all.equal(tcp, Z %*% t(Z))), "\n")
cat("  tcrossprod(Z) == Z%*%Conj(t(Z)):", isTRUE(all.equal(tcp, Z %*% Conj(t(Z)))), "\n")

cat("\n### Matrix package, for comparison ###\n")
suppressMessages(library(Matrix))
MZ <- Matrix(Z)
mcp <- Matrix::crossprod(MZ)
cat("  Matrix::crossprod(Z) == A^T A:", isTRUE(all.equal(as.matrix(mcp), tZZ, check.attributes = FALSE)), "\n")
cat("  Matrix::crossprod(Z) == A^H A:", isTRUE(all.equal(as.matrix(mcp), hZZ, check.attributes = FALSE)), "\n")
cat("  class:", class(mcp), "\n")

cat("\n### and what is Matrix's Hermitian-producing verb? ###\n")
cat("  Matrix::t is transpose:", isTRUE(all.equal(as.matrix(Matrix::t(MZ)), t(Z), check.attributes = FALSE)), "\n")
