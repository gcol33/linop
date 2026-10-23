# linop

*matrix-free linear operators*

linop represents a linear operator by the function that applies it, such
as a finite-difference stencil, a convolution or a projection, instead
of by the matrix of its entries.

You give [`linop()`](https://gillescolling.com/linop/reference/linop.md)
a function that applies the operator to a block of columns, and the
function that applies its conjugate transpose. The object that comes
back behaves like a matrix under `%*%`,
[`t()`](https://rdrr.io/r/base/t.html),
[`crossprod()`](https://rdrr.io/r/base/crossprod.html), `+` and `*`.
Each of these operations adds a node to an expression graph, so nothing
is computed until you apply the result to a vector or a block. An
operator also carries capabilities, such as being Hermitian or positive
definite, each with a record of the argument that established it, so you
can see what a claim rests on and decide what a solver may rely on.
linop supports:

- **Operators too large to store.** A 262,144 x 262,144 image blur, 512
  GiB as a dense matrix, is applied in 0.0095 seconds with a peak memory
  of 51 MB ([details](#operators-too-large-to-store)).
- **Solving systems and computing spectra.**
  [`solve()`](https://gillescolling.com/linop/reference/solve.linop.html)
  chooses among seven Krylov methods on the operator’s declared
  capabilities and shape, and
  [`eigs()`](https://gillescolling.com/linop/reference/eigs.html) and
  [`svds()`](https://gillescolling.com/linop/reference/svds.html)
  compute eigenpairs and singular triplets. Each result carries a
  certificate.
- **Capabilities with their evidence.** Properties propagate through the
  algebra together with the argument behind them, and a solver states
  the evidence it accepts as a
  [`requirement()`](https://gillescolling.com/linop/reference/requirement.html).
- **Checking a hand-written operator.**
  [`verify()`](https://gillescolling.com/linop/reference/verify.html)
  runs the contract every solver relies on, including the adjoint
  consistency test that catches a wrong conjugation.
- **Operators on infinite index sets.** A dimension may be `Inf`:
  [`linop_jacobi()`](https://gillescolling.com/linop/reference/linop_jacobi.html)
  builds a self-adjoint operator on the sequence space `l^2(Z)`, and
  [`eigs()`](https://gillescolling.com/linop/reference/eigs.md) runs on
  the operator itself.

Every solver is checked against closed-form truth. The propagation rules
for capabilities are checked by 9,892 brute-force soundness checks
against materialised matrices.

## Installation

``` r

# install.packages("pak")
pak::pak("gcol33/linop")
```

## Usage

1.  Build an operator with
    [`linop()`](https://gillescolling.com/linop/reference/linop.md) from
    a function and its adjoint, or from a matrix you already hold.
2.  Combine operators with ordinary matrix arithmetic. The result is a
    lazy expression.
3.  Apply it with `%*%`, or pass it to
    [`solve()`](https://rdrr.io/r/base/solve.html),
    [`eigs()`](https://gillescolling.com/linop/reference/eigs.md) or
    [`svds()`](https://gillescolling.com/linop/reference/svds.md).
4.  Check a hand-written operator with
    [`verify()`](https://gillescolling.com/linop/reference/verify.md).

The smallest operator is two functions and a shape. `S` shifts a vector
up by one position and pads with a zero; its conjugate transpose shifts
down:

``` r

library(linop)

n <- 8
up <- function(X) rbind(X[-1, , drop = FALSE], 0)
down <- function(X) rbind(0, X[-n, , drop = FALSE])
S <- linop(up, adjoint = down, dim = c(n, n))

S %*% 1:8
#> [1] 2 3 4 5 6 7 8 0
```

The function receives a block of columns, so a matrix goes through one
column at a time and a vector comes back as a vector, as with a matrix.

Arithmetic on operators builds an expression. The tree printed below is
the expression as written; nothing has been applied or allocated:

``` r

B <- crossprod(S) + 0.1 * linop_eye(n)
B
#> <linop> 8 x 8, double
#>   hermitian, symmetric, positive_semidefinite, real
#> 
#> sum
#> |- product
#> |  |- transpose
#> |  |  `- fun[8 x 8]
#> |  `- fun[8 x 8]
#> `- scale by 0.1
#>    `- I[8]
```

Applying `B` walks the tree right to left, one block at a time. A linear
system in `B` is solved with R’s own
[`solve()`](https://rdrr.io/r/base/solve.html), which returns the
solution with a certificate attached:

``` r

rhs <- as.numeric(B %*% rep(1, n))
x <- solve(B, rhs)
round(as.numeric(x), 10)
#> [1] 1 1 1 1 1 1 1 1
```

[Getting started with
linop](https://gillescolling.com/linop/articles/linop.html) builds
operators from matrices and sparse matrices, and covers the transpose,
adjoint and conjugate views and complex operators.

## Operators too large to store

A 512 x 512 image blur is written as the function that blurs one column
and reused as its own adjoint:

``` r

side <- 512L
n_img <- side * side

shift_up <- function(M) rbind(M[-1, , drop = FALSE], M[side, , drop = FALSE])
shift_down <- function(M) rbind(M[1, , drop = FALSE], M[-side, , drop = FALSE])

smooth <- function(v) {
  M <- matrix(v, side, side)
  M <- (M + shift_up(M) + shift_down(M)) / 3
  M <- t(M)
  M <- (M + shift_up(M) + shift_down(M)) / 3
  as.vector(t(M))
}

blur <- function(X) {
  Y <- X
  for (j in seq_len(ncol(X))) Y[, j] <- smooth(X[, j])
  Y
}

sym <- c(symmetric = TRUE, hermitian = TRUE)
Blur <- linop(blur, adjoint = blur, dim = c(n_img, n_img),
               properties = sym)
```

`Blur` is 262,144 x 262,144, which is 512 GiB written out as doubles.
The Tikhonov normal-equations operator for deblurring the image is one
line:

``` r

N <- crossprod(Blur) + 1e-3 * linop_eye(n_img)
img <- runif(n_img)
y <- N %*% img
```

|                                           |                   |
|-------------------------------------------|-------------------|
| Operator                                  | 262,144 x 262,144 |
| Dense storage for the same operator       | 512 GiB           |
| One apply of `Blur`                       | 0.0095 s          |
| One apply of `crossprod(Blur) + 1e-3 * I` | 0.019 s           |
| Peak R memory across the run              | 51 MB             |

The timings were taken with R 4.6.0 on Windows 11, an i9-14900K and the
reference BLAS. The expression is applied one block at a time, so the
memory figure is set by the block you hand it and by the operator’s own
working set. Wrapping a callback costs 200 ns per apply, which is 0.13%
or less of one application of an operator with 1e5 elements or more.

[The expression
graph](https://gillescolling.com/linop/articles/expressions.html) covers
the simplification rules, the cost model,
[`collapse()`](https://gillescolling.com/linop/reference/collapse.md)
and materialisation.

## Capabilities and their evidence

An operator records what is known about it together with the argument
that established it, and both travel through the algebra.
`adjoint(A) %*% A` is Hermitian whatever is claimed about `A`, so the
product node records it from its construction:

``` r

A <- linop(matrix(c(2, 1, 0, 1, 3, 1, 0, 1, 2), 3, 3))
G <- adjoint(A) %*% A
cap(G, "hermitian")
#> TRUE  (construction / identity / conf 1)
```

When two operators are declared positive definite by their author, their
sum is positive definite by a sound rule, and the declarations it rests
on stay visible underneath:

``` r

id <- function(X) X
twice <- function(X) 2 * X
pd <- c(hermitian = TRUE, positive_definite = TRUE)
P <- linop(id, adjoint = id, dim = c(4, 4), properties = pd)
Q <- linop(twice, adjoint = twice, dim = c(4, 4), properties = pd)

cap(P + Q, "positive_definite")
#> TRUE  (construction / identity / conf 1 <- [user_declaration / identity / conf 1; user_declaration / identity / conf 1])
```

A solver states what it accepts as a predicate over the whole chain.
This one accepts construction, computation and theorems, and refuses a
chain that ends in a user’s declaration:

``` r

req <- requirement(sources = c("construction", "computation", "theorem"))
evidence_satisfies(cap(P + Q, "positive_definite")$evidence, req)
#> [1] FALSE
```

Values are three-valued. `NA` means unknown and is never read as
`FALSE`, and a propagation rule that cannot prove `TRUE` returns `NA`.
[Capabilities and
evidence](https://gillescolling.com/linop/articles/capabilities.html)
lists the capabilities and the rules that propagate them.

## Checking an operator

[`verify()`](https://gillescolling.com/linop/reference/verify.md) runs
the contract that every solver in the package relies on, and reports
each check with its own status:

``` r

verify(S)
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> adjoint consistency            pass         computation      identity              1
#> view identities                pass         computation      identity              1
#> real transpose equals adjoint  pass         computation      identity              1
#> linearity                      pass         computation      identity              1
#> block consistency              pass         computation      identity              1
#> shape and dtype                pass         computation      identity              1
#> complex linearity              pass         computation      identity              1
#> declared capabilities          not_checked  probe            heuristic             -
#> purity                         pass         computation      identity              1
#> gemm agreement                 pass         computation      identity              1
#> materialisation agreement      pass         computation      identity              1
#> ------------------------------------------------------------------------------------------
#> overall                        qualified   no deterministic bound on: declared capabilities
```

The certificate keeps the checks that carry no deterministic bound in a
list of their own. A wrong conjugation in an adjoint is the most common
bug in a hand-written operator, and the adjoint consistency check fails
on it. Here the adjoint is a copy of the forward action:

``` r

Bad <- linop(up, adjoint = up, dim = c(n, n))
verify(Bad)$overall
#> [1] "fail"
```

An adapter for another storage format is one
[`linop()`](https://gillescolling.com/linop/reference/linop.md) method
for its class, followed by
[`verify()`](https://gillescolling.com/linop/reference/verify.md) on an
instance. [Writing an
adapter](https://gillescolling.com/linop/articles/adapters.html) shows
it, with node registration and the preconditioner contract.

## Solvers and eigensolvers

[`solve()`](https://rdrr.io/r/base/solve.html) reaches seven Krylov
methods for linear systems, and enforces the preconditioner contract in
each method that accepts a preconditioner.
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) is a
Hermitian eigensolver and
[`svds()`](https://gillescolling.com/linop/reference/svds.md) a singular
value solver. Each result carries a certificate, read with
[`certificate()`](https://gillescolling.com/linop/reference/certificate.md),
with the operator-norm estimate its arithmetic floor rests on.

[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) and
[`svds()`](https://gillescolling.com/linop/reference/svds.md) have two
engines. The reference engine in this package is a restarted Lanczos
iteration with full reorthogonalisation, which holds its whole basis in
each round. The [PRIMME](https://github.com/primme/primme) library
restarts implicitly and keeps a basis of fixed small size; it is
available when the companion package
[linop.primme](https://github.com/gcol33/linop.primme) is installed, and
the default `method = "auto"` then takes it. Whichever engine runs, the
reported values are measured on the operator and the certificate is the
same object, with a record of which engine ran and why.

[Solving a
system](https://gillescolling.com/linop/articles/solvers.html),
[Eigenvalues and singular
values](https://gillescolling.com/linop/articles/spectral.html) and
[Reading a
certificate](https://gillescolling.com/linop/articles/certificates.html)
cover each verb and what comes back.

## Operators on a sequence space

A dimension may be `Inf`.
[`linop_jacobi()`](https://gillescolling.com/linop/reference/linop_jacobi.md)
builds a self-adjoint operator on `l^2(Z)` whose essential spectrum is
exact by Weyl’s theorem and whose eigenvectors decay at a closed-form
rate, and
[`finite_section()`](https://gillescolling.com/linop/reference/finite_section.md)
truncates it.
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) runs on
the operator itself and chooses the truncation width: outside the window
the tail is exactly geometric, so the width that meets a tolerance is a
division. The certificate that comes back is a statement about the
operator on `l^2(Z)` itself. [Operators on a sequence
space](https://gillescolling.com/linop/articles/sequence-spaces.html)
works through an example.

## Related packages

- [eigencore](https://cran.r-project.org/package=eigencore) provides
  matrix-free block eigensolvers with an `as_operator()` adapter
  generic, complex dense arithmetic, block Krylov methods, shift-invert
  and certificates, through a constructor-object interface.
- [Matrix](https://cran.r-project.org/package=Matrix) provides stored
  sparse and structured matrix classes with direct factorisations for
  them.
- [sanic](https://cran.r-project.org/package=sanic) provides direct and
  Krylov solvers for sparse linear systems.
- [amatrix](https://github.com/bbuchsbaum/amatrix) provides
  backend-agnostic stored matrix classes, including GPU backends.

## Documentation

Start with [Getting started with
linop](https://gillescolling.com/linop/articles/linop.html): one path
from a pair of functions to an operator you apply, solve against and
check.

- [The expression
  graph](https://gillescolling.com/linop/articles/expressions.html):
  simplification rules, the cost model,
  [`collapse()`](https://gillescolling.com/linop/reference/collapse.md)
  and materialisation.
- [Capabilities and
  evidence](https://gillescolling.com/linop/articles/capabilities.html):
  what an operator knows about itself and how a solver states what it
  accepts.
- [Solving a
  system](https://gillescolling.com/linop/articles/solvers.html):
  [`solve()`](https://rdrr.io/r/base/solve.html), the seven Krylov
  methods behind it, and preconditioners.
- [Eigenvalues and singular
  values](https://gillescolling.com/linop/articles/spectral.html):
  [`eigs()`](https://gillescolling.com/linop/reference/eigs.md) and
  [`svds()`](https://gillescolling.com/linop/reference/svds.md).
- [Operators on a sequence
  space](https://gillescolling.com/linop/articles/sequence-spaces.html):
  operators with an infinite dimension and their truncation.
- [Reading a
  certificate](https://gillescolling.com/linop/articles/certificates.html):
  what comes back attached to every answer.
- [Writing an
  adapter](https://gillescolling.com/linop/articles/adapters.html):
  [`linop()`](https://gillescolling.com/linop/reference/linop.md)
  methods for your own storage format, node types and the preconditioner
  contract.
- [Full reference](https://gillescolling.com/linop/reference/): every
  function, grouped by job.

## Getting help

Questions and bug reports go to the [GitHub
issues](https://github.com/gcol33/linop/issues). A bug report is easiest
to act on with a small reproducible example, for instance made with
[reprex](https://reprex.tidyverse.org/).

## Contributing

Pull requests are welcome; see the [contributing
guide](https://github.com/gcol33/linop/blob/main/CONTRIBUTING.md).

## Support

> “Software is like sex: it’s better when it’s free.” — Linus Torvalds

I’m a PhD student who builds R packages in my free time because I
believe good tools should be free and open. I started these projects for
my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to
say thanks. It helps with my coffee addiction.

[![Buy Me A
Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

linop is licensed under MIT.

## Citation

``` text
Colling G (2026). _linop: Matrix-Free Linear Operators with
Evidence-Bearing Capabilities_. R package version 0.0.0.9000,
<https://gillescolling.com/linop/>.
```

BibTeX:

``` bibtex
@Manual{colling2026linop,
  title = {linop: Matrix-Free Linear Operators with Evidence-Bearing Capabilities},
  author = {Gilles Colling},
  year = {2026},
  note = {R package version 0.0.0.9000},
  url = {https://gillescolling.com/linop/},
}
```

Download the citation as
[BibTeX](https://gillescolling.com/linop/linop.bib) or
[RIS](https://gillescolling.com/linop/linop.ris).
