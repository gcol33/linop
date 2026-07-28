## Smoke run for the section 7.2 reference solvers, against closed forms.
devtools::load_all(".", quiet = TRUE)

lap <- function(n) {
  linop(function(X) rbind(X[-1, , drop = FALSE], 0) +
                    rbind(0, X[-n, , drop = FALSE]) - 2 * X,
        dim = c(n, n), properties = c(hermitian = TRUE))
}
lap_eig <- function(n) -4 * sin(seq_len(n) * pi / (2 * (n + 1)))^2

n <- 60
A <- lap(n)
truth <- sort(lap_eig(n))

cat("--- eigs, smallest algebraic ---\n")
fit <- eigs(A, k = 4, which = "smallest_algebraic")
print(fit)
cat("gap to truth: ", max(abs(fit$values - truth[1:4])), "\n")

cat("\n--- eigs, largest algebraic ---\n")
fit2 <- eigs(A, k = 4, which = "largest_algebraic")
cat("values: ", fit2$values, "\n")
cat("gap to truth: ", max(abs(sort(fit2$values) - sort(truth[(n - 3):n]))), "\n")

cat("\n--- eigs, largest magnitude ---\n")
fit3 <- eigs(A, k = 3, which = "largest")
cat("values: ", fit3$values, "\n")
cat("gap to truth: ", max(abs(sort(fit3$values) - sort(truth[1:3]))), "\n")

cat("\n--- certificate ---\n")
print(fit$certificate)

cat("\n--- dense symmetric, exact norm ---\n")
set.seed(1)
M <- matrix(rnorm(25), 5, 5); M <- M + t(M)
D <- linop(M)
fd <- eigs(D, k = 2, which = "largest")
ev <- eigen(M)$values
cat("values: ", fd$values, " truth: ", ev[order(-abs(ev))][1:2], "\n")
print(fd$certificate)

cat("\n--- complex hermitian ---\n")
set.seed(2)
Z <- matrix(rnorm(36) + 1i * rnorm(36), 6, 6); Z <- Z + Conj(t(Z))
CZ <- linop(Z)
fz <- eigs(CZ, k = 2, which = "largest")
ez <- Re(eigen(Z)$values)
cat("values: ", fz$values, " truth: ", ez[order(-abs(ez))][1:2], "\n")
cat("imaginary part of vectors nonzero: ", is.complex(fz$vectors), "\n")

cat("\n--- shift-invert ---\n")
fs <- eigs(A, k = 3, sigma = -1)
cat("values: ", fs$values, "\n")
cat("nearest to -1 in truth: ", truth[order(abs(truth + 1))][1:3], "\n")
print(fs$certificate)

cat("\n--- svds, dense rectangular ---\n")
set.seed(3)
R <- matrix(rnorm(40), 8, 5)
SR <- linop(R)
fv <- svds(SR, k = 3)
cat("d: ", fv$d, " truth: ", svd(R)$d[1:3], "\n")
print(fv$certificate)

cat("\n--- svds, matrix free ---\n")
nd <- 50
Dm <- linop(function(X) X[-1, , drop = FALSE] - X[-nd, , drop = FALSE],
            adjoint = function(X) rbind(0, X) - rbind(X, 0),
            dim = c(nd - 1L, nd))
dsv <- 2 * sin(seq.int(0, nd - 1L) * pi / (2 * nd))
fm <- svds(Dm, k = 3)
cat("d: ", fm$d, " truth: ", rev(sort(dsv))[1:3], "\n")
cat("U orthonormal: ", max(abs(crossprod(fm$u) - diag(3))), "\n")
cat("V orthonormal: ", max(abs(crossprod(fm$v) - diag(3))), "\n")

cat("\n--- svds, complex ---\n")
set.seed(4)
Cm <- matrix(rnorm(30) + 1i * rnorm(30), 6, 5)
fc <- svds(linop(Cm), k = 2)
cat("d: ", fc$d, " truth: ", svd(Cm)$d[1:2], "\n")

cat("\nall smoke runs finished\n")
