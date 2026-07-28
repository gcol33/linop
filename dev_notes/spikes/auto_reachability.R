## Which operators can reach which branch of method = "auto", given the evidence
## minima of plan section 7.1.
##
##   Rscript dev_notes/spikes/auto_reachability.R

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-linop.R")

n <- 12
Mspd <- spd_prescribed(n, seq(1, 4, length.out = n))
Mhpd <- hpd_prescribed(n, seq(1, 4, length.out = n))

ops <- list(
  "linop_eye" = linop_eye(n),
  "linop_scaling(+d)" = linop_scaling(seq(1, 4, length.out = n)),
  "linop_scaling(mixed)" = linop_scaling(c(-1, seq(1, 4, length.out = n - 1L))),
  "linop(dense spd), bare" = linop(Mspd),
  "linop(dense spd), declared" = linop(Mspd, properties = c(hermitian = TRUE,
                                                            positive_definite = TRUE)),
  "linop(dense hpd), bare" = linop(Mhpd),
  "laplacian_1d (fun leaf)" = laplacian_1d(n),
  "shifted_laplacian_1d" = shifted_laplacian_1d(n, 0.9),
  "convdiff_1d" = convdiff_1d(n, 0.3),
  "adjoint(X) %*% X" = adjoint(linop(rmat(2 * n, n, seed = 1))) %*%
                        linop(rmat(2 * n, n, seed = 1)),
  "kms as_spd_linop" = as_spd_linop(kms_matrix(n, 0.4)),
  "sum of two declared" = linop(Mspd, properties = c(hermitian = TRUE,
                                                     positive_definite = TRUE)) +
                          linop(Mspd, properties = c(hermitian = TRUE,
                                                     positive_definite = TRUE)))

src <- function(A, nm) {
  e <- linop:::cape(A, nm)
  if (is.null(e)) "-" else e$source
}
val <- function(A, nm) format(capv(A, nm))

cat(sprintf("%-28s %-10s %-16s %-10s %-16s %-8s %s\n",
            "operator", "herm", "herm source", "pd", "pd source", "diag", "auto picks"))
cat(strrep("-", 118), "\n", sep = "")
for (nm in names(ops)) {
  A <- ops[[nm]]
  pick <- linop:::auto_method(A)$method
  cat(sprintf("%-28s %-10s %-16s %-10s %-16s %-8s %s\n", nm,
              val(A, "hermitian"), src(A, "hermitian"),
              val(A, "positive_definite"), src(A, "positive_definite"),
              val(A, "diagonal"), pick))
}

cat("\n=== does any operator in the package reach the cg branch of auto ===\n")
reach <- vapply(ops, function(A) linop:::auto_method(A)$method, character(1))
cat(sprintf("  direct   %s\n", paste(names(reach)[reach == "direct"], collapse = ", ")))
cat(sprintf("  cg       %s\n", paste(names(reach)[reach == "cg"], collapse = ", ")))
cat(sprintf("  minres   %s\n", paste(names(reach)[reach == "minres"], collapse = ", ")))
cat(sprintf("  gmres    %s\n", paste(names(reach)[reach == "gmres"], collapse = ", ")))

cat("\n=== what arithmetic does to the certificate attribute ===\n")
A <- linop(Mspd, properties = c(hermitian = TRUE, positive_definite = TRUE))
set.seed(1)
b <- stats::rnorm(n)
x <- solve(A, b, method = "cg", tol = 1e-11)
keeps <- function(lab, e) cat(sprintf("  %-22s certificate %s\n", lab,
                                      if (is.null(attr(e, "certificate"))) "dropped" else "kept"))
keeps("x", x)
keeps("x + 0", x + 0)
keeps("2 * x", 2 * x)
keeps("x + x", x + x)
keeps("x[1:3]", x[1:3])
keeps("as.numeric(x)", as.numeric(x))
keeps("c(x)", c(x))
keeps("sum(x)", sum(x))
keeps("crossprod-ish t(x)", t(x))
X <- solve(A, cbind(b, b), method = "cg", tol = 1e-11)
keeps("matrix X", X)
keeps("X + 0", X + 0)
keeps("X[, 1]", X[, 1])
