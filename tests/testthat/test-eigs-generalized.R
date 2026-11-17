## eigs(A, B = ) for the generalized problem A x = lambda B x: Lanczos in the B
## inner product for a definite pencil, Arnoldi on B^-1 A otherwise, shift-invert
## on (A - sigma B)^-1 B, and the certificate for a pencil.

## The linear-element pencil of the 1-D Laplacian: stiffness K = tridiag(-1, 2, -1)
## and mass M = tridiag(1, 4, 1) / 6 share the sine eigenvectors, so the pencil's
## eigenvalues are 6 (1 - cos t_j) / (2 + cos t_j), t_j = j pi / (n + 1).
fem_mass <- function(n) {
  M <- diag(4, n)
  M[cbind(1:(n - 1), 2:n)] <- 1
  M[cbind(2:n, 1:(n - 1))] <- 1
  M / 6
}
fem_eigenvalues <- function(n) {
  t <- seq_len(n) * pi / (n + 1)
  sort(6 * (1 - cos(t)) / (2 + cos(t)))
}

pencil_residual <- function(A, B, X, values) {
  R <- linop:::linop_apply(A, X, "N") -
    linop:::linop_apply(B, X, "N") * rep(values, each = nrow(X))
  max(sqrt(colSums(Mod(R)^2)) / sqrt(colSums(Mod(X)^2)))
}

nearest_error <- function(values, truth) {
  max(vapply(values, function(v) min(Mod(truth - v)), numeric(1)))
}

## ----------------------------------------------------------- recovery --

test_that("a definite pencil runs Lanczos in the B inner product, at both ends", {
  n <- 100
  K <- laplacian_1d(n)
  M <- linop(fem_mass(n), properties = c(positive_definite = TRUE))
  truth <- fem_eigenvalues(n)
  top <- eigs(K, k = 3, which = "largest_algebraic", B = M, method = "lanczos")
  expect_equal(top$method, "lanczos")
  expect_equal(top$values, rev(truth)[1:3], tolerance = 1e-12)
  expect_lt(pencil_residual(K, M, top$vectors, top$values), 1e-9)
  ## B-orthonormal, not orthonormal
  G <- crossprod(top$vectors, fem_mass(n) %*% top$vectors)
  expect_lt(max(abs(G - diag(3))), 1e-12)
  expect_match(top$certificate$dispatch$reason, "hermitian pencil with B positive definite")
  expect_match(top$certificate$dispatch$reason, "on B\\^-1 A")

  low <- eigs(K, k = 3, which = "smallest_algebraic", B = M, method = "lanczos")
  expect_equal(low$nconv, 3L)
  expect_equal(low$values, truth[1:3], tolerance = 1e-12)
  expect_equal(cert_status(low$certificate, "residual"), "pass")
})

test_that("B = I gives the standard problem's values", {
  n <- 60
  A <- laplacian_1d(n)
  ## linop_eye() establishes definiteness; a dense identity computes symmetry only
  I <- linop_eye(n)
  std <- eigs(A, k = 3, which = "largest_algebraic")
  gen <- eigs(A, k = 3, which = "largest_algebraic", B = I, method = "lanczos")
  expect_equal(gen$method, "lanczos")
  expect_equal(gen$values, std$values, tolerance = 1e-12)
})

test_that("random definite pencils, real and complex, against eigen() over 20 seeds", {
  for (s in 1:20) {
    n <- 40
    dt <- if (s %% 2) "double" else "complex"
    Ad <- hpd_prescribed(n, seq(-3, 5, length.out = n), seed = s)
    if (dt == "double") Ad <- Re(Ad)
    Bd <- hpd_prescribed(n, exp(seq(log(10), 0, length.out = n)), seed = 100L + s)
    if (dt == "double") Bd <- Re(Bd)
    Ad <- (Ad + Conj(t(Ad))) / 2; Bd <- (Bd + Conj(t(Bd))) / 2
    ref <- sort(Re(eigen(solve(Bd, Ad), only.values = TRUE)$values), decreasing = TRUE)
    A <- linop(Ad)
    B <- linop(Bd, properties = c(positive_definite = TRUE))
    fit <- eigs(A, k = 3, which = "largest_algebraic", B = B, seed = s, method = "lanczos")
    info <- sprintf("seed %d %s", s, dt)
    expect_equal(fit$method, "lanczos", info = info)
    expect_equal(fit$values, ref[1:3], tolerance = 1e-10, info = info)
    expect_false(fit$certificate$overall == "fail", info = info)
    expect_equal(cert_status(fit$certificate, "orthogonality"), "pass", info = info)
    expect_equal(cert_status(fit$certificate, "forward error"), "pass", info = info)
  }
})

test_that("a non-hermitian pencil runs Arnoldi on B^-1 A, against eigen() over 20 seeds", {
  for (s in 1:20) {
    n <- 30
    Ad <- rmat(n, n, seed = s)
    Bd <- diag(n) + 0.2 * rmat(n, n, seed = 50L + s)
    ev <- eigen(solve(Bd, Ad), only.values = TRUE)$values
    fit <- eigs(linop(Ad), k = 4, which = "largest", B = linop(Bd), seed = s)
    info <- sprintf("seed %d", s)
    expect_equal(fit$method, "arnoldi", info = info)
    expect_lt(nearest_error(fit$values, ev), 1e-9, label = info)
    expect_equal(sort(Mod(fit$values), decreasing = TRUE),
                 sort(Mod(ev), decreasing = TRUE)[1:4], tolerance = 1e-9, info = info)
    expect_equal(cert_status(fit$certificate, "backward error"), "pass", info = info)
    expect_equal(cert_status(fit$certificate, "forward error"), "not_checked", info = info)
  }
})

test_that("a B that declares nothing takes Arnoldi, and the dispatch says why", {
  n <- 60
  Md <- fem_mass(n)
  M <- linop(function(X) Md %*% X, dim = c(n, n))
  ## A computed hermitian, so the refusal is about B
  fit <- eigs(linop(as.matrix(laplacian_1d(n))), k = 3, which = "largest_real", B = M)
  expect_equal(fit$method, "arnoldi")
  expect_equal(fit$values, rev(fem_eigenvalues(n))[1:3], tolerance = 1e-10)
  cand <- fit$certificate$dispatch$candidates
  expect_match(cand$reason[cand$backend == "lanczos"], "the B inner product requires B to be hermitian")
  expect_false(grepl("\n", cand$reason[cand$backend == "lanczos"]))
  expect_error(eigs(laplacian_1d(n), k = 3, B = M, method = "lanczos"),
               "the B inner product requires B")
  expect_equal(cert_status(fit$certificate, "orthogonality"), "not_checked")
})

test_that("shift-invert on (A - sigma B)^-1 B, definite and not", {
  n <- 100
  K <- laplacian_1d(n)
  M <- linop(fem_mass(n), properties = c(positive_definite = TRUE))
  truth <- fem_eigenvalues(n)
  near <- truth[order(abs(truth - 0.5))][1:3]
  fit <- eigs(K, k = 3, sigma = 0.5, B = M, method = "lanczos")
  expect_equal(fit$method, "lanczos shift-invert")
  expect_match(fit$certificate$dispatch$reason, "\\(A - 0.5 B\\)\\^-1 B")
  expect_equal(sort(fit$values), sort(near), tolerance = 1e-12)

  Ad <- rmat(30, 30, seed = 7L); Bd <- diag(30) + 0.2 * rmat(30, 30, seed = 8L)
  ev <- eigen(solve(Bd, Ad), only.values = TRUE)$values
  sigma <- 0.3 + 0.2i
  fz <- eigs(linop(Ad), k = 2, sigma = sigma, B = linop(Bd))
  expect_equal(fz$method, "arnoldi shift-invert")
  expect_lt(nearest_error(fz$values, ev[order(Mod(ev - sigma))][1:2]), 1e-9)
})

## ---------------------------------------------------------- certificate --

test_that("the pencil's backward error perturbs both operators", {
  n <- 60
  K <- laplacian_1d(n)
  Md <- fem_mass(n)
  M <- linop(Md, properties = c(positive_definite = TRUE))
  fit <- eigs(K, k = 2, which = "largest_algebraic", B = M)
  cert <- fit$certificate
  X <- fit$vectors
  r <- sqrt(colSums((as.matrix(K) %*% X - (Md %*% X) * rep(fit$values, each = n))^2))
  x <- sqrt(colSums(X^2))
  nA <- cert$values$norm; nB <- cert$values$norm_B
  expect_equal(cert$values$backward_error, r / ((nA + abs(fit$values) * nB) * x),
               tolerance = 1e-10)
  expect_match(cert$checks$detail[cert$checks$check == "backward error"],
               "\\(\\|\\|A\\|\\| \\+ \\|theta\\| \\|\\|B\\|\\|\\)")
  expect_match(cert$checks$detail[cert$checks$check == "residual"], "theta B x")
})

test_that("the forward bound rests on the declarations of both operators", {
  n <- 40
  K <- laplacian_1d(n)
  Md <- fem_mass(n)
  declared <- linop(function(X) Md %*% X, dim = c(n, n),
                    properties = c(hermitian = TRUE, positive_definite = TRUE))
  fit <- eigs(K, k = 2, which = "largest_algebraic", B = declared)
  ev <- fit$certificate$evidence[["forward error"]]
  expect_equal(ev$guarantee, "estimate")
  sources <- vapply(ev$depends_on, function(d) d$source, character(1))
  expect_true("user_declaration" %in% sources)
  expect_false(evidence_satisfies(ev, requirement(sources = c("theorem", "construction",
                                                             "computation"))))
  expect_lt(fit$certificate$values$forward_bound[1], 1e-8)
})

test_that("a hermitian pencil with B not definite has no real-spectrum guarantee", {
  ## A = diag(1, -1), B = [[0, 1], [1, 0]] padded: hermitian both, B indefinite,
  ## and the pencil's eigenvalues include +-i.
  n <- 6
  Ad <- diag(c(1, -1, 3, 4, 5, 6))
  Bd <- diag(n); Bd[1:2, 1:2] <- matrix(c(0, 1, 1, 0), 2)
  ev <- eigen(solve(Bd, Ad), only.values = TRUE)$values
  fit <- eigs(linop(Ad), k = 2, which = "largest_imaginary", B = linop(Bd), ncv = 6L)
  expect_equal(fit$method, "arnoldi")
  expect_equal(sort(Im(fit$values)), c(-1, 1), tolerance = 1e-10)
  expect_lt(nearest_error(fit$values, ev), 1e-10)
})

## ------------------------------------------------------- preconditioned --

## A mass matrix whose density rho_j spans `spread`, so kappa(B) is near it and
## ||B^-1|| is far from 1/||B||: the case a Jacobi preconditioner on the inner
## solve is for. The pencil's eigenvalues come from the symmetric reduction
## L^-1 K L^-T with B = L L^T.
graded_mass <- function(n, spread) {
  rho <- spread^seq(0, 1, length.out = n)
  M <- diag(4 * rho)
  off <- sqrt(rho[-1] * rho[-n])
  M[cbind(1:(n - 1), 2:n)] <- off
  M[cbind(2:n, 1:(n - 1))] <- off
  M / 6
}
graded_truth <- function(Kd, Md) {
  Li <- solve(t(chol(Md)))
  sort(eigen(Li %*% Kd %*% t(Li), symmetric = TRUE, only.values = TRUE)$values)
}
counted_mass <- function(Md, env) {
  linop(function(X) { env$B <- env$B + ncol(X); Md %*% X },
        adjoint = function(X) Md %*% X, dim = dim(Md),
        properties = c(hermitian = TRUE, positive_definite = TRUE))
}
counted_jacobi <- function(d, env) {
  preconditioner(function(X) { env$P <- env$P + 1L; X / d },
                 hermitian = TRUE, positive_definite = TRUE)
}

test_that("a preconditioner reaches the solve with B through inner", {
  n <- 100
  Md <- graded_mass(n, 1e4)
  truth <- rev(graded_truth(as.matrix(laplacian_1d(n)), Md))[1:3]
  plain <- new.env(); plain$B <- 0L
  pre <- new.env(); pre$B <- 0L; pre$P <- 0L
  f0 <- eigs(laplacian_1d(n), k = 3, which = "largest_algebraic",
             B = counted_mass(Md, plain), method = "lanczos",
             inner = list(method = "cg"))
  f1 <- eigs(laplacian_1d(n), k = 3, which = "largest_algebraic",
             B = counted_mass(Md, pre), method = "lanczos",
             inner = list(method = "cg", preconditioner = counted_jacobi(diag(Md), pre)))
  expect_gt(pre$P, 0L)
  for (f in list(f0, f1)) {
    expect_equal(f$nconv, 3L)
    expect_equal(f$values, truth, tolerance = 1e-12)
    expect_equal(cert_status(f$certificate, "residual"), "pass")
    expect_equal(cert_status(f$certificate, "forward error"), "pass")
  }
  ## the same answer at a fifth of the applies of B, at the same budget
  expect_identical(f0$maxit, f1$maxit)
  expect_lt(pre$B, plain$B / 3)
})

test_that("a preconditioner reaches the shift-invert solve with A - sigma B", {
  n <- 100
  Kd <- as.matrix(laplacian_1d(n))
  Md <- graded_mass(n, 1e4)
  truth <- graded_truth(Kd, Md)
  sigma <- (truth[5] + truth[6]) / 2
  env <- new.env(); env$B <- 0L; env$P <- 0L
  pc <- counted_jacobi(abs(diag(Kd - sigma * Md)), env)
  fit <- eigs(laplacian_1d(n), k = 2, sigma = sigma, B = counted_mass(Md, env),
              method = "lanczos", inner = list(method = "minres", preconditioner = pc))
  expect_gt(env$P, 0L)
  expect_equal(fit$nconv, 2L)
  expect_equal(sort(fit$values), truth[5:6], tolerance = 1e-10)
  expect_false(fit$certificate$overall == "fail")
})

test_that("the forward row allows for ||B^-1|| beyond 1/||B||, and never past it", {
  n <- 60
  Kd <- as.matrix(laplacian_1d(n))
  Md <- graded_mass(n, 1e6)
  truth <- graded_truth(Kd, Md)
  sigma <- (truth[5] + truth[6]) / 2
  M <- linop(Md, properties = c(positive_definite = TRUE))
  fit <- eigs(laplacian_1d(n), k = 2, sigma = sigma, B = M, method = "lanczos")
  cert <- fit$certificate
  cv <- cert$values
  expect_equal(sort(fit$values), truth[5:6], tolerance = 1e-12)
  ## read with ||B^-1|| as 1/||B||, this converged pair fails the row; the
  ## bound is at rounding level, so by a factor that varies between 7 and 35
  expect_gt(max(cv$forward_bound / (1e-10 * (cv$norm / cv$norm_B + abs(fit$values)))), 2)
  expect_equal(cert_status(cert, "forward error"), "pass")
  detail <- cert$checks$detail[cert$checks$check == "forward error"]
  binv <- as.numeric(sub(".*\\|\\|B\\^-1\\|\\| >= ([^;]+);.*", "\\1", detail))
  true_binv <- 1 / min(eigen(Md, symmetric = TRUE, only.values = TRUE)$values)
  ## a lower bound on ||B^-1||, and here a long way above 1/||B||
  expect_lte(binv, true_binv * (1 + 1e-6))
  expect_gt(binv, 100 / cert$values$norm_B)
})

## ------------------------------------------------------------ refusals --

test_that("what eigs() refuses about B", {
  n <- 10
  A <- laplacian_1d(n)
  M <- linop(fem_mass(n), properties = c(positive_definite = TRUE))
  expect_error(eigs(A, k = 2, B = linop(diag(n + 1))), "the two the same shape")
  expect_error(eigs(A, k = 2, B = "B"), "B must be a linop or a matrix")
  expect_error(eigs(A, k = 2, B = M, which = "largest_imaginary"),
               "pencil is established")
  expect_error(eigs(linop_jacobi(diagonal = 1), k = 1, B = M), "finite dimension")
  ## a plain matrix is accepted and computed symmetric, but definiteness is not
  ## computed, so it goes to Arnoldi
  fit <- eigs(A, k = 2, which = "largest_real", B = fem_mass(n))
  expect_equal(fit$method, "arnoldi")
  expect_equal(fit$values, rev(fem_eigenvalues(n))[1:2], tolerance = 1e-10)
})

test_that("the external engines refuse B by name", {
  skip_if_not_installed("linop.primme")
  M <- linop(fem_mass(10), properties = c(positive_definite = TRUE))
  expect_error(eigs(laplacian_1d(10), k = 2, B = M, method = "primme"), "does not take B")
})

test_that("auto takes Lanczos for a pencil only when both declarations meet the minimum", {
  n <- 40
  K <- linop(as.matrix(laplacian_1d(n)))
  declared <- linop(fem_mass(n), properties = c(positive_definite = TRUE))
  fit <- eigs(K, k = 2, which = "largest_algebraic", B = declared)
  expect_equal(fit$method, "arnoldi")
  cand <- fit$certificate$dispatch$candidates
  expect_match(cand$reason[cand$backend == "lanczos"],
               "the B inner product auto acts on positive_definite only")
  ## linop_eye() establishes definiteness by construction
  eye <- eigs(K, k = 2, which = "largest_algebraic", B = linop_eye(n))
  expect_equal(eye$method, "lanczos")
})
