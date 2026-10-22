# Story: linop

## Claim
Matrix-free linear operators for R: an operator defined by its action behaves like a matrix, composes lazily, and carries capabilities with the evidence that established them.

## Audience
R users who solve systems or compute spectra of operators too large or too structured to store (stencils, convolutions, Jacobians), and authors of operator classes who want an adapter. They know base R matrix arithmetic and basic numerical linear algebra.

## Relation
- to: none
- phrase: none
- home: README, Related packages (eigencore, Matrix, sanic, amatrix)

## Arcs

### A1 Behaves like a matrix, without forming one
- question: can I use my function as if it were the matrix?
- answer: `linop()` from an action and its adjoint; `%*%`, `t()`, `crossprod()`, `+`, `*` build a lazy expression
- payoff: operators far larger than memory, applied block by block
- evidence: README blur example; vignettes/expressions.Rmd
- limits: scalar addition and elementwise products are not defined
- home: vignettes/linop.Rmd (Get started), vignettes/expressions.Rmd
- also: README, DESCRIPTION

### A2 Capabilities carry their evidence
- question: can a solver trust that this operator is Hermitian or positive definite?
- answer: three-valued capabilities with an evidence chain; solvers state a `requirement()`
- payoff: the reader decides which claims a solver may rely on
- evidence: propagation test suite (brute-force soundness checks)
- home: vignettes/capabilities.Rmd
- also: README

### A3 Solving and spectra with certificates
- question: can I solve and compute eigenvalues on it?
- answer: `solve()` over seven Krylov methods, `eigs()`, `svds()`, each with a certificate
- limits: `eigs()` and `svds()` are reference implementations (no implicit restart)
- home: vignettes/solvers.Rmd, vignettes/spectral.Rmd, vignettes/certificates.Rmd
- also: README

### A4 Checking an operator
- question: is my hand-written adjoint right?
- answer: `verify()` runs the conformance contract
- home: vignettes/adapters.Rmd
- also: README, Get started

### A5 Infinite index sets
- question: can the operator live on l^2(Z)?
- answer: `linop_jacobi()`, `finite_section()`, `eigs()` choosing the width
- home: vignettes/sequence-spaces.Rmd
- also: README

## Non-goals
Stored sparse or dense matrix classes and direct factorisations (Matrix covers them).

## Documents
- DESCRIPTION: claim, A1, A2, conformance suite
- README / landing: claim, A1 to A5 one line each, first example, routing
- Get started (linop.Rmd): one path from two closures to apply, solve and verify
- articles: expressions, capabilities, solvers, spectral, sequence-spaces, certificates (explanation/how-to); adapters (how-to)

## Open
- Break line kept from the previous README ("matrix-free linear operators").
