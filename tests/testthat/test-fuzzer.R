## Gate 1: generated expression trees checked against a base R recomputation
## from materialised leaves.
##
## The gate calls for 1e4 trees. That is what LINOP_FUZZ_N controls, and the full
## count is run and recorded in dev_notes/GATE1.md. The default here is smaller
## so routine R CMD check stays inside CRAN's time budget; nothing about the
## check changes with the count, only how much of the space it covers.

fuzz_n <- as.integer(Sys.getenv("LINOP_FUZZ_N", "400"))

## Build a random expression of a requested shape, returning the linop and the
## dense reference side by side.
gen_expr <- function(m, n, depth, cplx) {
  leaf <- function() {
    if (m == n && runif(1) < 0.25) {
      if (runif(1) < 0.5) {
        A <- linop_eye(m); return(list(A = A, M = diag(1, m)))
      }
      d <- if (cplx) complex(real = rnorm(m), imaginary = rnorm(m)) else rnorm(m)
      return(list(A = linop_scaling(d), M = diag(d, nrow = m)))
    }
    M <- if (cplx) zmat(m, n) else rmat(m, n)
    list(A = linop(M), M = M)
  }
  if (depth <= 0) return(leaf())

  op <- sample(c("leaf", "sum", "product", "transpose", "adjoint", "conjugate", "scale"), 1)
  switch(op,
    leaf = leaf(),
    sum = {
      a <- gen_expr(m, n, depth - 1L, cplx); b <- gen_expr(m, n, depth - 1L, cplx)
      list(A = a$A + b$A, M = a$M + b$M)
    },
    product = {
      k <- sample(1:4, 1)
      a <- gen_expr(m, k, depth - 1L, cplx); b <- gen_expr(k, n, depth - 1L, cplx)
      list(A = a$A %*% b$A, M = a$M %*% b$M)
    },
    transpose = {
      a <- gen_expr(n, m, depth - 1L, cplx)
      list(A = t(a$A), M = t(a$M))
    },
    adjoint = {
      a <- gen_expr(n, m, depth - 1L, cplx)
      list(A = adjoint(a$A), M = Conj(t(a$M)))
    },
    conjugate = {
      a <- gen_expr(m, n, depth - 1L, cplx)
      list(A = Conj(a$A), M = Conj(a$M))
    },
    scale = {
      s <- if (cplx && runif(1) < 0.5) complex(real = rnorm(1), imaginary = rnorm(1)) else rnorm(1)
      a <- gen_expr(m, n, depth - 1L, cplx)
      list(A = s * a$A, M = s * a$M)
    })
}

test_that("generated expression trees agree with base R", {
  set.seed(20260727)
  worst <- 0; worst_at <- ""
  for (i in seq_len(fuzz_n)) {
    cplx <- runif(1) < 0.5
    m <- sample(1:4, 1); n <- sample(1:4, 1)
    e <- gen_expr(m, n, depth = sample(1:3, 1), cplx = cplx)

    got <- as.matrix(e$A)
    ref <- e$M
    scale <- max(max(Mod(ref)), 1e-12)
    rel <- max(Mod(got - ref)) / scale
    if (rel > worst) { worst <- rel; worst_at <- sprintf("tree %d, %dx%d, complex=%s", i, m, n, cplx) }

    ## and the action must agree, not only the materialisation
    X <- if (cplx) zmat(n, 2) else rmat(n, 2)
    rel2 <- max(Mod((e$A %*% X) - (ref %*% X))) / scale
    worst <- max(worst, rel2)
  }
  expect_lt(worst, 1e-12)
  if (worst >= 1e-12) cat("worst at:", worst_at, "\n")
})

test_that("the fuzzer covers every composite node type", {
  set.seed(11)
  seen <- character()
  collect <- function(op) {
    seen <<- c(seen, op$node)
    for (k in linop:::node_children(op)) collect(k)
  }
  for (i in seq_len(200)) {
    e <- gen_expr(sample(1:3, 1), sample(1:3, 1), depth = 3L, cplx = runif(1) < 0.5)
    collect(e$A)
  }
  for (nd in c("sum", "product", "transpose", "adjoint", "conjugate", "scale", "dense")) {
    expect_true(nd %in% seen, info = paste("fuzzer never produced node", nd))
  }
})
