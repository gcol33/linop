# linop

*matrix-free linear operators*

[![R-CMD-check](https://github.com/gcol33/linop/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/linop/actions/workflows/R-CMD-check.yaml)
[![Dependencies](https://img.shields.io/badge/dependencies-0-brightgreen)](https://github.com/gcol33/linop/blob/main/DESCRIPTION)
[![R >= 4.4.0](https://img.shields.io/badge/R-%3E%3D%204.4.0-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Linear operators defined by their action rather than their entries.**

Give `linop()` a function that applies your operator to a block of columns, plus
the function that applies its conjugate transpose, and you get an object that
behaves like a matrix. `%*%`, `t()`, `crossprod()`, `+` and `*` all work, and
every one of them builds an expression, so nothing is formed until you ask for a
result.

```r
library(linop)

n <- 2000
S <- linop(
  function(X) rbind(X[-1, , drop = FALSE], 0),           # shift up
  adjoint = function(X) rbind(0, X[-n, , drop = FALSE]), # shift down
  dim = c(n, n)
)

S %*% rnorm(n)
B <- crossprod(S) + 0.1 * linop_eye(n)   # 2000 x 2000, nothing allocated
```

## The size that stops being a problem

A 512 x 512 image blur, written as the function that blurs one column and reused
as its own adjoint:

```r
side <- 512L
n <- side * side

up   <- function(M) rbind(M[-1, , drop = FALSE], M[side, , drop = FALSE])
down <- function(M) rbind(M[1, , drop = FALSE], M[-side, , drop = FALSE])

smooth <- function(v) {
  M <- matrix(v, side, side)
  M <- (M + up(M) + down(M)) / 3
  M <- t(M)
  M <- (M + up(M) + down(M)) / 3
  as.vector(t(M))
}

blur <- function(X) {
  Y <- X
  for (j in seq_len(ncol(X))) Y[, j] <- smooth(X[, j])
  Y
}

B <- linop(blur, adjoint = blur, dim = c(n, n),
           properties = c(symmetric = TRUE, hermitian = TRUE))
```

`B` is 262,144 x 262,144. Written out as doubles, that matrix is 512 GiB. The
Tikhonov normal-equations operator for deblurring the image is one line, and
applying it takes 19 milliseconds:

```r
N <- crossprod(B) + 1e-3 * linop_eye(n)
N %*% img
```

| | |
|---|---|
| Operator | 262,144 x 262,144 |
| Dense storage for the same operator | 512 GiB |
| One apply of `B` | 0.0095 s |
| One apply of `crossprod(B) + 1e-3 * I` | 0.019 s |
| Peak R memory across the run | 51 MB |

R 4.6.0, Windows 11, i9-14900K, reference BLAS. The expression is applied right
to left, one block at a time, so the memory figure is set by the block you hand
it and by the operator's own working set.

Wrapping a callback costs 200 ns per apply, measured during design against a
bare closure call. Against an operator of 1e5 elements or more that is 0.13% of
one application or less, which is why the package carries no compiled code.

## Capabilities carry their evidence

An operator records what is known about it together with the argument that
established it. Both halves travel through the algebra.

```r
A <- linop(matrix(c(2, 1, 0, 1, 3, 1, 0, 1, 2), 3, 3))
G <- adjoint(A) %*% A

cap(G, "hermitian")
#> TRUE  (construction / identity / conf 1)
```

The evidence has nothing underneath it because `adjoint(X) %*% X` is Hermitian
whatever anyone claims about `X`. The product node recognises the shape and
stamps it unconditionally.

Now the case the design exists for. Two operators are declared positive definite
by their author, and their sum is positive definite by a rule that is
unconditionally sound:

```r
P <- linop(function(X) X, adjoint = function(X) X, dim = c(4, 4),
           properties = c(hermitian = TRUE, positive_definite = TRUE))
Q <- linop(function(X) 2 * X, adjoint = function(X) 2 * X, dim = c(4, 4),
           properties = c(hermitian = TRUE, positive_definite = TRUE))

cap(P + Q, "positive_definite")
#> TRUE  (construction / identity / conf 1 <- [user_declaration / identity / conf 1;
#>                                             user_declaration / identity / conf 1])
```

The rule is sound, so the value is `TRUE`. What the rule rests on stays visible
underneath it. A solver states what it will accept as a predicate over the whole
chain:

```r
req <- requirement(sources = c("construction", "computation", "theorem"))
evidence_satisfies(cap(P + Q, "positive_definite")$evidence, req)
#> FALSE
```

Values are three-valued. `NA` means unknown and is never read as `FALSE`. A
propagation rule that cannot prove `TRUE` returns `NA`.

## Checking an operator

`verify()` answers one question: is this object's claim valid. It runs the
contract that every solver in the package relies on.

```r
verify(B)
```

```
check                          status       source           guarantee             conf
------------------------------------------------------------------------------------------
adjoint consistency            pass         computation      identity              1
view identities                pass         computation      identity              1
real transpose equals adjoint  pass         computation      identity              1
linearity                      pass         computation      identity              1
block consistency              pass         computation      identity              1
shape and dtype                pass         computation      identity              1
complex linearity              pass         computation      identity              1
declared capabilities          pass         probe            heuristic             -
purity                         pass         computation      identity              1
gemm agreement                 pass         computation      identity              1
materialisation agreement      pass         computation      identity              1
------------------------------------------------------------------------------------------
overall                        pass   no deterministic bound on: declared capabilities
```

Every line reports its own status, and the certificate keeps the checks that
carry no deterministic bound in a list of their own. A wrong conjugation in an
adjoint is the most common bug in a hand-written operator, and it lands on the
first line:

```
adjoint consistency            fail         computation      identity              1
...
failures:
  adjoint consistency: max relative gap 9.546e-01
```

## Adapting your own storage format

Implement one method. Everything else follows from it.

```r
linop.my_format <- function(x, ...) {
  linop(
    function(X) my_multiply(x, X),
    adjoint = function(X) my_multiply_adjoint(x, X),
    dim = c(my_nrow(x), my_ncol(x))
  )
}
```

Then run `verify()` on an instance and read the certificate. Storage that does
not fit the callback shape registers a node type instead, through the same door
the built-in nodes use:

```r
linop_register_node("scaled_identity",
  apply       = function(op, X, mode) op$args$a * X,
  materialize = function(op) diag(op$args$a, op$dim[1L]),
  cost        = function(op) op$dim[1L])
```

`linop_nodes()` lists what is registered. An operator that came from a continuous
problem can carry an opaque envelope back to its provider with
`set_provenance()`, which the core never looks inside.

## Installation

```r
install.packages("pak")
pak::pak("gcol33/linop")
```

Installing `linop` installs one package. There are no `Imports` and no compiled
code, so there is no toolchain to have working first. `Matrix` sits in `Suggests`
and is used for sparse leaves and `as_sparse()`.

R 4.4.0 or newer. The floor is set by `crossprod()` and `tcrossprod()`, which
became S3 generic in R 4.4.0, and are what let one `matrixOps` method cover all
three matrix products.

## Documentation

- [Getting started](https://gillescolling.com/linop/articles/linop.html)
- [The expression graph](https://gillescolling.com/linop/articles/expressions.html)
- [Capabilities and evidence](https://gillescolling.com/linop/articles/capabilities.html)
- [Writing an adapter](https://gillescolling.com/linop/articles/adapters.html)

## Status

The finite object model, the expression graph, the capability and evidence
lattices, the adapters and the conformance suite are in place, with 10,403
assertions across 120 tests. The propagation suite alone is 9,892 brute-force
soundness checks against materialised matrices.

The solver layer is under way. Conjugate gradients, the certificate that belongs
to a result rather than to an operator, and the operator-norm estimate its
arithmetic floor rests on are written and checked against closed-form truth. The
remaining six Krylov methods and the `solve()`, `eigs()` and `svds()` front doors
follow, and the front doors open once `method = "auto"` has a choice to make.

## Related packages

- [eigencore](https://cran.r-project.org/package=eigencore) provides matrix-free
  block eigensolvers with an `as_operator()` adapter generic, complex dense
  arithmetic, block Krylov methods, shift-invert and certificates, through a
  constructor-object interface.
- [Matrix](https://cran.r-project.org/package=Matrix) provides stored sparse and
  structured matrix classes with direct factorisations for them.
- [sanic](https://cran.r-project.org/package=sanic) provides direct and Krylov
  solvers for sparse linear systems.
- [amatrix](https://github.com/bbuchsbaum/amatrix) provides backend-agnostic
  stored matrix classes, including GPU backends.

## Support

> "Software is like sex: it's better when it's free." - Linus Torvalds

I'm a PhD student who builds R packages in my free time because I believe good
tools should be free and open. I started these projects for my own work and
figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to say
thanks. It helps with my coffee addiction.

[![Buy Me A Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

MIT (see the LICENSE file)

## Citation

```bibtex
@software{linop,
  author = {Colling, Gilles},
  title  = {linop: Matrix-Free Linear Operators with Evidence-Bearing Capabilities},
  year   = {2026},
  url    = {https://github.com/gcol33/linop}
}
```
