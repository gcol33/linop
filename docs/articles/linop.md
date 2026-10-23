# Getting started with linop

## Installation

``` r

install.packages("pak")
pak::pak("gcol33/linop")
```

Installing `linop` installs one package. There are no `Imports` and no
compiled code, so there is no toolchain to have working first. `Matrix`
sits in `Suggests` and is used for sparse leaves and
[`as_sparse()`](https://gillescolling.com/linop/reference/as_sparse.md).

R 4.4.0 or newer. The floor is set by
[`crossprod()`](https://rdrr.io/r/base/crossprod.html) and
[`tcrossprod()`](https://rdrr.io/r/base/crossprod.html), which became S3
generic in R 4.4.0, and are what let one `matrixOps` method cover all
three matrix products.

## What a linop is

A matrix in R is a rectangle of numbers. Many of the linear maps that
turn up in practice are described more directly by what they do than by
the numbers they would contain: a finite-difference stencil, a
convolution, a projection onto a subspace, the Jacobian of a simulation,
a Gram matrix of a design you already hold in another form. Writing out
the entries is possible for all of them and affordable for few.

A `linop` is that map as an object. You supply the action, `linop`
supplies everything a matrix would have given you.

``` r

n <- 8
S <- linop(
  function(X) rbind(X[-1, , drop = FALSE], 0),
  adjoint = function(X) rbind(0, X[-n, , drop = FALSE]),
  dim = c(n, n)
)
S
#> <linop> 8 x 8, double
#>   real
#> 
#> fun[8 x 8]
```

`S` shifts a vector up by one position and pads with a zero. Its
conjugate transpose shifts down. Two closures and a shape, and the
object is complete.

``` r

v <- 1:8
S %*% v
#> [1] 2 3 4 5 6 7 8 0
```

A vector goes in and a vector comes out, as it would with a matrix.
Handing it a matrix gives a matrix, one column at a time:

``` r

S %*% cbind(1:8, rep(1, 8))
#>      [,1] [,2]
#> [1,]    2    1
#> [2,]    3    1
#> [3,]    4    1
#> [4,]    5    1
#> [5,]    6    1
#> [6,]    7    1
#> [7,]    8    1
#> [8,]    0    0
```

## Building one

[`linop()`](https://gillescolling.com/linop/reference/linop.md) is a
single generic dispatching on its first argument, with methods for the
inputs below and a clear error for anything else.

A base matrix wraps the storage you already have:

``` r

M <- matrix(c(2, 1, 0,
              1, 3, 1,
              0, 1, 2), 3, 3)
A <- linop(M)
A
#> <linop> 3 x 3, double
#>   hermitian, symmetric, not diagonal, not triangular_upper, not triangular_lower, real
#> 
#> dense[3 x 3]
```

Notice what it says about itself. The dense leaf checks symmetry,
diagonality and triangularity exactly, on data it already holds, and
records the outcome. That check is
[`all.equal()`](https://rdrr.io/r/base/all.equal.html) against the
transpose, not a guess.

A `Matrix` from the **Matrix** package wraps sparse storage:

``` r

Sp <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 2, 3), x = c(2, 3, 4))
linop(Sp) %*% c(1, 1, 1)
#> [1] 2 3 4
```

A function is the callback form used above. `dim` is required, `adjoint`
is optional, and omitting the adjoint means the transposed modes raise
an error naming what is missing rather than returning something wrong.

``` r

Twice <- linop(function(X) 2 * X, dim = c(4, 4))
t(Twice) %*% rep(1, 4)
#> Error:
#> ! this operator has no adjoint; supply adjoint = to linop() to use mode 'T'
```

An object with no
[`linop()`](https://gillescolling.com/linop/reference/linop.md) method
gets an error that names the method to write:

``` r

linop(list(1, 2))
#> Error:
#> ! no linop() method for class list
#>   Adapter authors: implement linop.list() returning a linop, then run verify().
```

Two named constructors cover the operators that appear constantly inside
expressions:

``` r

linop_eye(4)
#> <linop> 4 x 4, double
#>   hermitian, symmetric, positive_definite, positive_semidefinite, unitary, diagonal, triangular_upper, triangular_lower, full_rank, real
#> 
#> I[4]
linop_scaling(c(2, 4, 8))
#> <linop> 3 x 3, double
#>   hermitian, symmetric, positive_definite, positive_semidefinite, diagonal, triangular_upper, triangular_lower, full_rank, real
#> 
#> diag[3]
```

## The algebra is lazy

Every arithmetic operation returns a new operator. Nothing is applied,
nothing is allocated.

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

The tree printed underneath is the expression as written. `crossprod(S)`
became a product of a transpose view and the original leaf; the scaled
identity became a `scale` over `I`; the sum holds both. Applying `B`
walks that tree right to left, one block at a time.

``` r

B %*% rep(1, n)
#> [1] 0.1 1.1 1.1 1.1 1.1 1.1 1.1 1.1
```

Everything R offers for matrices is routed through one `matrixOps` group
method, so all three products behave the same way:

``` r

X <- matrix(rnorm(n * 2), n, 2)
all.equal(crossprod(S) %*% X, t(S) %*% (S %*% X))
#> [1] TRUE
all.equal(tcrossprod(S) %*% X, S %*% (t(S) %*% X))
#> [1] TRUE
```

Two operations are deliberately absent, because for matrices they mean
something this package cannot mean. Adding a scalar to an operator would
have to guess between `A + c` elementwise and `A + cI`, and `*` between
two operators is the elementwise product:

``` r

A + 1
#> Error:
#> ! adding a scalar to an operator is ambiguous.
#>   For A + c*I write  A + c * linop_eye(nrow(A))
```

``` r

A * A
#> Error:
#> ! '*' between two operators is elementwise for matrices and is not defined here.
#>   For the matrix product write  A %*% B
```

## Transpose, adjoint and conjugate are three things

For a real operator, the transpose and the conjugate transpose have the
same action. They stay separate objects anyway, so that the expression
you printed is the expression you wrote, and so that a real subgraph
inside a complex expression keeps the distinction at the point where it
starts to matter.

``` r

Rect <- linop(matrix(1:6, 2, 3))
t(Rect)$node
#> [1] "transpose"
adjoint(Rect)$node
#> [1] "adjoint"
Conj(Rect)$node
#> [1] "conjugate"
```

Composing them collapses through the identities you would expect, at
construction time:

``` r

t(t(Rect))$node              # involution
#> [1] "dense"
Conj(adjoint(Rect))$node     # adjoint then conjugate is the transpose
#> [1] "transpose"
adjoint(t(Rect))$node        # transpose then adjoint is the conjugate
#> [1] "conjugate"
```

Internally these are four apply modes in the BLAS sense, composed
through a Klein four-group: `N` for `A X`, `T` for `A^T X`, `C` for
`A^H X`, `R` for `conj(A) X`. A callback leaf supplies two of them and
the other two are derived by conjugation, which is why nobody writes
four functions.

## Complex operators, and one thing to know about crossprod

An operator is `"double"` or `"complex"`. No path silently drops an
imaginary part: a real operator applied to a complex block promotes the
result.

``` r

dtype(A)
#> [1] "double"
A %*% c(1 + 2i, 0, 0)
#> [1] 2+4i 1+2i 0+0i
```

Base R’s `crossprod(A)` is `A^T A`, the transpose and not the conjugate
transpose. For a complex operator that product is symmetric and is
generally not Hermitian. `linop` follows base R, because the promise of
the package is that a `linop` behaves like a matrix:

``` r

Z <- linop(matrix(c(1 + 1i, 2, 0, 3 - 2i), 2, 2))
cap(crossprod(Z), "symmetric")$value
#> [1] TRUE
cap(crossprod(Z), "hermitian")$value
#> [1] NA
```

The Hermitian Gram matrix is spelled out:

``` r

cap(adjoint(Z) %*% Z, "hermitian")$value
#> [1] TRUE
```

## Looking at an expression

[`print()`](https://rdrr.io/r/base/print.html) shows the tree.
[`explain()`](https://gillescolling.com/linop/reference/explain.md) adds
the shape, dtype and cost of every node:

``` r

explain(B)
#> operator 8 x 8, dtype double, estimated 264 flops per column
#> 7 nodes, depth 3
#> 
#>  depth      node nrow ncol  dtype cost
#>      0       sum    8    8 double  264
#>      1   product    8    8 double  256
#>      2 transpose    8    8 double  128
#>      3       fun    8    8 double  128
#>      2       fun    8    8 double  128
#>      1     scale    8    8 double    8
#>      2  identity    8    8 double    0
```

The cost column is flops per column, and it is an ordering heuristic.
Measurement during design found operators where a wider block is
*slower* per column, which a flop count cannot express, so this number
orders candidate plans and never predicts wall time.

[`summary()`](https://rdrr.io/r/base/summary.html) reports the
capability table with the evidence behind each entry:

``` r

summary(B)$caps[, c("capability", "value")]
#>               capability value
#> 1              hermitian  TRUE
#> 2              symmetric  TRUE
#> 3      positive_definite    NA
#> 4  positive_semidefinite  TRUE
#> 5                unitary    NA
#> 6               diagonal    NA
#> 7       triangular_upper    NA
#> 8       triangular_lower    NA
#> 9              full_rank    NA
#> 10                  real  TRUE
```

## Forming a matrix, when you mean to

Materialisation is an explicit request with a guard, because the entire
point of the object is that it normally does not happen.

``` r

as.matrix(A)
#>      [,1] [,2] [,3]
#> [1,]    2    1    0
#> [2,]    1    3    1
#> [3,]    0    1    2
round(as.matrix(B), 3)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]  0.1  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#> [2,]  0.0  1.1  0.0  0.0  0.0  0.0  0.0  0.0
#> [3,]  0.0  0.0  1.1  0.0  0.0  0.0  0.0  0.0
#> [4,]  0.0  0.0  0.0  1.1  0.0  0.0  0.0  0.0
#> [5,]  0.0  0.0  0.0  0.0  1.1  0.0  0.0  0.0
#> [6,]  0.0  0.0  0.0  0.0  0.0  1.1  0.0  0.0
#> [7,]  0.0  0.0  0.0  0.0  0.0  0.0  1.1  0.0
#> [8,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.1
```

The guard is on the entry count and refuses politely:

``` r

big <- linop(function(X) X, adjoint = function(X) X, dim = c(1e5, 1e5))
as.matrix(big)
#> Error:
#> ! refusing to materialise a 100000 x 100000 operator (1e+10 entries, about 74.5 GB).
#>   Raise the guard with as.matrix(A, max_entries = ...) if you mean it.
```

[`collapse()`](https://gillescolling.com/linop/reference/collapse.md) is
the middle path. It walks the expression and replaces subexpressions
that are cheaper to hold as a dense leaf, and it reports what it did, so
a wrong cost model costs a suboptimal plan and never a surprise
allocation.

``` r

K <- linop(matrix(rnorm(200 * 200), 200, 200))
E <- crossprod(K) + 0.5 * linop_eye(200)
linop_cost(E)
#> [1] 160200

Ec <- collapse(E, expected_applies = 500)
linop_cost(Ec)
#> [1] 80000
```

``` r

report <- attr(Ec, "report")
class(report) <- "linop_collapse_report"
report
#> collapsed:     sum of 200 x 200
#> retained:      nothing
#> memory added:  312.5 KB
#> break-even:    1 applications
```

## Checking that it is what it claims

[`verify()`](https://gillescolling.com/linop/reference/verify.md) runs
the contract. It is the function to reach for after writing a
hand-rolled operator, and the one an adapter author runs before
shipping.

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

The interesting line is the first. `adjoint consistency` tests
`<y, A x> = <A^H y, x>` over random probes, which is what catches a
missing or misplaced conjugation. Here is the same operator with its
adjoint replaced by a copy of the forward action:

``` r

Bad <- linop(function(X) rbind(X[-1, , drop = FALSE], 0),
             adjoint = function(X) rbind(X[-1, , drop = FALSE], 0),
             dim = c(6, 6))
verify(Bad)$overall
#> [1] "fail"
```

The full certificate, and what each check is for, is covered in
[`vignette("adapters")`](https://gillescolling.com/linop/articles/adapters.md).

## Solving

An operator you can apply is an operator you can solve against.
[`solve()`](https://rdrr.io/r/base/solve.html) is R’s own generic, so it
needs no new verb:

``` r

Lap <- linop(function(X) 2 * X - rbind(X[-1, , drop = FALSE], 0) -
                         rbind(0, X[-8, , drop = FALSE]), dim = c(8, 8))
rhs <- as.numeric(Lap %*% rep(1, 8))
round(solve(Lap, rhs), 12)
#> [1] 1 1 1 1 1 1 1 1
#> attr(,"certificate")
#> solved by gmres (no capability is established, and gmres requires none)
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         computation      estimate              -
#> residual                       pass         computation      identity              1
#> backward error                 pass         computation      estimate              -
#> convergence                    pass         computation      identity              1
#> forward error                  not_checked  computation      identity              -
#> ------------------------------------------------------------------------------------------
#> overall                        qualified   no deterministic bound on: arithmetic floor, backward error, forward error
```

[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) and
[`svds()`](https://gillescolling.com/linop/reference/svds.md) are the
other two, for eigenpairs and singular triplets. All three return a
certificate alongside the answer.

## Where to next

- [`vignette("expressions")`](https://gillescolling.com/linop/articles/expressions.md)
  covers the expression graph: the simplification rules, the cost model,
  [`collapse()`](https://gillescolling.com/linop/reference/collapse.md)
  and materialisation.

- [`vignette("capabilities")`](https://gillescolling.com/linop/articles/capabilities.md)
  covers what an operator knows about itself and how a solver states
  what it will accept.

- [`vignette("solvers")`](https://gillescolling.com/linop/articles/solvers.md)
  covers [`solve()`](https://rdrr.io/r/base/solve.html) and the seven
  Krylov methods behind it, with preconditioners.

- [`vignette("spectral")`](https://gillescolling.com/linop/articles/spectral.md)
  covers [`eigs()`](https://gillescolling.com/linop/reference/eigs.md)
  and [`svds()`](https://gillescolling.com/linop/reference/svds.md).

- [`vignette("certificates")`](https://gillescolling.com/linop/articles/certificates.md)
  covers what comes back attached to every answer.

- [`vignette("adapters")`](https://gillescolling.com/linop/articles/adapters.md)
  covers writing `linop.myclass()`, registering a node type, and the
  preconditioner contract.
