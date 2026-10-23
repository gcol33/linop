# Capabilities and evidence

A solver needs to know things about an operator that it cannot afford to
find out. Conjugate gradients converges on a Hermitian positive definite
operator and can silently produce nonsense on one that is neither, and
checking positive definiteness costs more than the solve. So the
operator has to carry the claim.

The moment it carries a claim, a second question arrives: on whose word.
A property established by construction, one measured exactly on stored
data, one asserted by whoever wrote the adapter, and one guessed from
ten random probes are all “positive definite is TRUE”, and they are not
interchangeable. This vignette is about keeping the two apart and
keeping both correct through the algebra.

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

## Ten capabilities, three values

``` r

caps <- c("hermitian", "symmetric", "positive_definite", "positive_semidefinite",
          "unitary", "diagonal", "triangular_upper", "triangular_lower",
          "full_rank", "real")
```

`hermitian` and `symmetric` are separate entries. For a complex operator
`A = A^H` and `A = A^T` are different statements, and a solver that
takes one for the other takes a wrong branch. For a real operator they
coincide, and the package records that coincidence as a rule instead of
merging the fields.

Each capability holds `TRUE`, `FALSE` or `NA`. `NA` means unknown, and
it is never read as `FALSE`.

``` r

S <- linop(function(X) X[c(2, 1), , drop = FALSE],
           adjoint = function(X) X[c(2, 1), , drop = FALSE],
           dim = c(2, 2))
cap(S, "positive_definite")
#> NA  (<none>)
```

That operator is a swap. It is Hermitian and it is not positive
definite, and `linop` says neither, because nothing has established
either one. Propagation is conservative in the same direction: a rule
that cannot prove `TRUE` returns `NA` and never `FALSE`.

[`cap()`](https://gillescolling.com/linop/reference/cap.md) reads one
entry. [`summary()`](https://rdrr.io/r/base/summary.html) reads them
all:

``` r

A <- linop(matrix(c(2, 1, 0,
                    1, 3, 1,
                    0, 1, 2), 3, 3))
summary(A)$caps[, c("capability", "value")]
#>               capability value
#> 1              hermitian  TRUE
#> 2              symmetric  TRUE
#> 3      positive_definite    NA
#> 4  positive_semidefinite    NA
#> 5                unitary    NA
#> 6               diagonal FALSE
#> 7       triangular_upper FALSE
#> 8       triangular_lower FALSE
#> 9              full_rank    NA
#> 10                  real  TRUE
```

The dense leaf filled those in by checking the stored matrix exactly.
Symmetry came from `all.equal(M, t(M))`, diagonality from testing the
off-diagonal entries. Positive definiteness is absent because
establishing it would have cost a factorisation the constructor did not
run.

## Evidence is three independent fields

``` r

evidence(source = "computation", guarantee = "identity", confidence = 1)
#> computation / identity / conf 1
```

**source** is where the claim came from: `construction`,
`adapter_contract`, `user_declaration`, `probe`, `computation`,
`theorem`.

**guarantee** is the kind of statement being made: `identity`,
`deterministic_bound`, `probabilistic_bound`, `estimate`, `heuristic`.

**confidence** is `1`, a stated probability, or `NA`.

The three are independent on purpose, and none of them is a rank. A
randomised trace estimator carries a `probabilistic_bound` with a
rigorous confidence; a deterministic heuristic carries no bound at all;
an adapter contract can be worth more than a hundred probes. Collapsing
that into one “strength” number loses exactly the distinctions a solver
needs.

The fourth field is what makes it work: `depends_on`, a list of the
evidence a rule relied on. It is empty when the rule is unconditional.

## Declaring what you know

`properties` takes a named logical vector and records each entry as a
user declaration.

``` r

P <- linop(function(X) X, adjoint = function(X) X, dim = c(4, 4),
           properties = c(hermitian = TRUE, positive_definite = TRUE))
cap(P, "positive_definite")
#> TRUE  (user_declaration / identity / conf 1)
```

A bare character vector is shorthand for all-`TRUE`, and an unknown name
is an error naming the ten that exist.

``` r

linop(matrix(1:4, 2, 2), properties = c(hermitiann = TRUE))
#> Error:
#> ! unknown properties: hermitiann
#>   known: hermitian, symmetric, positive_definite, positive_semidefinite, unitary, diagonal, triangular_upper, triangular_lower, full_rank, real
```

Some implications hold by definition and are closed once at
construction, so declaring the strong form does not oblige you to also
declare the weak one:

``` r

cap(P, "positive_semidefinite")
#> TRUE  (construction / identity / conf 1 <- [user_declaration / identity / conf 1])
```

Positive definite implies positive semidefinite; diagonal implies both
triangularities; a real Hermitian operator is symmetric and the
converse. Each closure records what it depended on. A `"double"`
operator is `real` by construction, unconditionally.

``` r

cap(linop(matrix(1:4, 2, 2)), "real")
#> TRUE  (construction / identity / conf 1)
```

Declaring `hermitian` or `symmetric` on a non-square operator is an
error rather than a stored contradiction:

``` r

linop(matrix(1:6, 2, 3), properties = c(symmetric = TRUE))
#> Error:
#> ! a hermitian or symmetric operator must be square
```

## Propagation

Every composite recomputes its capabilities from its operands. The rules
are sound in the direction that matters: they claim `TRUE` only where
the conclusion follows, and answer `NA` everywhere else.

**Sums.** `hermitian`, `symmetric`, `positive_definite`,
`positive_semidefinite`, `real`, `diagonal` and both triangularities
survive addition when every term has them.

``` r

D1 <- linop_scaling(c(1, 2, 3))
D2 <- linop_scaling(c(4, 5, 6))
cap(D1 + D2, "diagonal")$value
#> [1] TRUE
cap(D1 + D2, "positive_definite")$value
#> [1] TRUE
```

**Scaling.** A real scalar preserves `hermitian`; any scalar preserves
`symmetric`; a positive real scalar preserves definiteness. A nonzero
scalar preserves `full_rank`, and structure such as `diagonal` survives
any scalar.

``` r

Pd <- linop(function(X) X, adjoint = function(X) X, dim = c(3, 3),
            properties = c(hermitian = TRUE, positive_definite = TRUE))

cap(2 * Pd, "positive_definite")$value
#> [1] TRUE
cap(-2 * Pd, "positive_definite")$value    # the rule cannot conclude
#> [1] NA
cap(-2 * Pd, "hermitian")$value            # a real scalar keeps this one
#> [1] TRUE
```

A diagonal leaf answers `FALSE` to the same question, and the difference
is worth understanding. Scaling a `diag` is a simplification rule: the
result is a fresh diagonal leaf, whose constructor computes definiteness
exactly from the stored vector rather than inheriting a propagated
claim.

``` r

cap(-2 * D1, "positive_definite")$value
#> [1] FALSE
(-2 * D1)$node
#> [1] "diag"
```

**Transpose, adjoint, conjugate.**
[`adjoint()`](https://gillescolling.com/linop/reference/adjoint.md)
preserves `hermitian`, definiteness, `unitary`, `diagonal`, `full_rank`
and `real`, and preserves `symmetric` when the operand is also real.
[`t()`](https://rdrr.io/r/base/t.html) preserves `symmetric` and
preserves `hermitian` when the operand is real. Both flip the
triangularities. [`Conj()`](https://rdrr.io/r/base/complex.html)
preserves everything.

``` r

U <- linop(matrix(c(1, 0, 2, 1), 2, 2))      # upper triangular
cap(U, "triangular_upper")$value
#> [1] TRUE
cap(t(U), "triangular_upper")$value
#> [1] NA
cap(t(U), "triangular_lower")$value
#> [1] TRUE
```

**Products.** A product of two arbitrary operators tells you nothing,
and the rules say so. Two exceptions are structural, and they hold
whatever is known about the operand:

``` r

K <- linop(matrix(rnorm(12), 4, 3))
G <- adjoint(K) %*% K
cap(G, "hermitian")
#> TRUE  (construction / identity / conf 1)
cap(G, "positive_semidefinite")$value
#> [1] TRUE
cap(G, "positive_definite")$value            # NA until K is known full rank
#> [1] NA
```

`adjoint(X) %*% X` is Hermitian and positive semidefinite for every `X`,
so the evidence has nothing underneath it. Positive definiteness needs
`X` to have full rank, so it appears only when that is known, and then
it records the dependency.

The product node finds this by structural recognition, comparing the two
factors rather than testing the result. It is conservative: it may fail
to notice that two separately built operators are equal, and it never
claims that two different operators are the same.

``` r

Kr <- linop(matrix(rnorm(12), 4, 3))
cap(crossprod(Kr), "symmetric")$value        # t(A) A is symmetric for any A
#> [1] TRUE
cap(crossprod(Kr), "hermitian")$value        # and Hermitian too, since A is real
#> [1] TRUE
```

For a complex operator the two come apart, which is the point of keeping
them separate:

``` r

Kz <- linop(matrix(c(1 + 1i, 2, 0, 3 - 2i), 2, 2))
cap(crossprod(Kz), "symmetric")$value        # A^T A is symmetric
#> [1] TRUE
cap(crossprod(Kz), "hermitian")$value        # and is not Hermitian
#> [1] NA
cap(adjoint(Kz) %*% Kz, "hermitian")$value   # A^H A is
#> [1] TRUE
```

Structure that survives any product is carried through: a product of
diagonals is diagonal, of upper triangulars is upper triangular, of
unitaries is unitary, of reals is real.

``` r

cap(D1 %*% D2, "diagonal")$value
#> [1] TRUE
```

## The case this design exists for

Take two operators whose author declared them positive definite. The sum
rule is unconditionally sound. So the sum is positive definite, by
construction:

``` r

P <- linop(function(X) X, adjoint = function(X) X, dim = c(4, 4),
           properties = c(hermitian = TRUE, positive_definite = TRUE))
Q <- linop(function(X) 2 * X, adjoint = function(X) 2 * X, dim = c(4, 4),
           properties = c(hermitian = TRUE, positive_definite = TRUE))

cap(P + Q, "positive_definite")
#> TRUE  (construction / identity / conf 1 <- [user_declaration / identity / conf 1; user_declaration / identity / conf 1])
```

Read the evidence. The source is `construction`, which is true and is
not the whole story: underneath it sit the two user declarations the
rule consumed. Had the composite reported a bare `construction`, two
unverified assertions would have been laundered into a fact by one
addition, and a solver reading only the top-level source would accept
them.

## Requirements

A solver states what it will accept. The statement is a predicate over
evidence, not a threshold on a strength number.

``` r

strict <- requirement(sources = c("construction", "computation", "theorem"),
                      guarantees = "identity")
```

[`evidence_satisfies()`](https://gillescolling.com/linop/reference/evidence_satisfies.md)
checks the evidence against it, and recurses into `depends_on`, so a
requirement the inputs would have failed is failed by anything built
from them.

``` r

evidence_satisfies(cap(P, "positive_definite")$evidence, strict)
#> [1] FALSE
evidence_satisfies(cap(P + Q, "positive_definite")$evidence, strict)
#> [1] FALSE
```

The Gram matrix passes the same requirement, because there is nothing
underneath its evidence to fail it:

``` r

evidence_satisfies(cap(G, "positive_semidefinite")$evidence, strict)
#> [1] TRUE
```

A solver willing to take the author’s word says so explicitly, and then
the same chain passes:

``` r

lenient <- requirement(sources = c("construction", "computation", "theorem",
                                   "user_declaration"))
evidence_satisfies(cap(P + Q, "positive_definite")$evidence, lenient)
#> [1] TRUE
```

`min_confidence` adds a numeric floor, applied at every level of the
chain:

``` r

probabilistic <- evidence("probe", "probabilistic_bound", confidence = 0.99)
evidence_satisfies(probabilistic, requirement(min_confidence = 0.999))
#> [1] FALSE
evidence_satisfies(probabilistic, requirement(min_confidence = 0.95))
#> [1] TRUE
```

This is the shape of the contract Phase 2 solvers will use: CG will
require Hermitian positive definiteness under a requirement it
publishes, and an operator whose claim rests on an unverified
declaration will be turned away by name rather than converge to
something wrong.

## Probing a declaration

[`verify()`](https://gillescolling.com/linop/reference/verify.md) tests
declared capabilities against the operator’s actual behaviour. Agreement
is recorded beside the declaration as a probe and never promoted to
proof; a contradiction fails the certificate.

``` r

Good <- linop(function(X) X, adjoint = function(X) X, dim = c(4, 4),
              properties = c(hermitian = TRUE))
verify(Good)$probes
#>   capability claimed gap    verdict
#> 1  hermitian    TRUE   0 consistent
#> 2  symmetric    TRUE   0 consistent
```

``` r

Liar <- linop(function(X) rbind(X[-1, , drop = FALSE], 0),
              adjoint = function(X) rbind(0, X[-4, , drop = FALSE]),
              dim = c(4, 4), properties = c(hermitian = TRUE))
verify(Liar)$probes
#>   capability claimed      gap      verdict
#> 1  hermitian    TRUE 1.664972 contradicted
#> 2  symmetric    TRUE 3.203087 contradicted
verify(Liar)$overall
#> [1] "fail"
```

Positive definiteness is the one case where a probe can only ever
disagree. A positive Rayleigh quotient on ten random vectors is not
evidence of definiteness, so the check looks for a negative one and
reports `consistent` otherwise. The certificate records the guarantee as
`heuristic` for exactly this reason.

## Where to next

- [`vignette("adapters")`](https://gillescolling.com/linop/articles/adapters.md)
  covers the rest of the certificate, and the contracts that
  preconditioners and solve objects carry.

- [`vignette("expressions")`](https://gillescolling.com/linop/articles/expressions.md)
  covers the tree these capabilities travel up.

- [`vignette("solvers")`](https://gillescolling.com/linop/articles/solvers.md)
  covers what `method = "auto"` does with the evidence above, which is
  where these three fields decide something visible.

- [`vignette("certificates")`](https://gillescolling.com/linop/articles/certificates.md)
  covers the same three fields underneath a computed result, including a
  bound that inherits the evidence of the capability it rests on.
