# Package index

## Build an operator

One generic dispatching on its first argument, plus named constructors
for the operators that appear constantly inside expressions.

- [`linop()`](https://gillescolling.com/linop/reference/linop.md) :
  Build a linear operator
- [`linop_eye()`](https://gillescolling.com/linop/reference/linop_eye.md)
  : Identity operator
- [`linop_scaling()`](https://gillescolling.com/linop/reference/linop_scaling.md)
  : Diagonal operator
- [`is_linop()`](https://gillescolling.com/linop/reference/is_linop.md)
  : Is this a linop?
- [`dtype()`](https://gillescolling.com/linop/reference/dtype.md) : The
  scalar type of an operator

## Operators on a sequence space

An operator whose index set is infinite, the truncation that makes it
computable, and the rate its eigenvectors decay at. What comes back from
the truncation is an ordinary operator, and what comes back from a
spectral computation on it is a statement about the operator underneath.

- [`linop_jacobi()`](https://gillescolling.com/linop/reference/linop_jacobi.md)
  : A Jacobi operator on the sequence space
- [`finite_section()`](https://gillescolling.com/linop/reference/finite_section.md)
  : Truncate an operator on a sequence space
- [`decay_rate()`](https://gillescolling.com/linop/reference/decay_rate.md)
  : How fast an eigenvector decays

## Views

The conjugate transpose, distinct from the transpose.

- [`adjoint()`](https://gillescolling.com/linop/reference/adjoint.md) :
  The conjugate transpose

## Inspect and cost an expression

Read the structure of a lazy expression, estimate its cost, and trade
cost for memory where that pays.

- [`explain()`](https://gillescolling.com/linop/reference/explain.md) :
  Report the structure and cost of an expression
- [`linop_cost()`](https://gillescolling.com/linop/reference/linop_cost.md)
  : Estimated flops per column
- [`collapse()`](https://gillescolling.com/linop/reference/collapse.md)
  : Collapse cheap subexpressions into dense leaves

## Form a matrix, deliberately

Explicit materialisation, guarded on size.

- [`as.matrix(`*`<linop>`*`)`](https://gillescolling.com/linop/reference/as.matrix.linop.md)
  : Materialise an operator as a dense matrix
- [`as_sparse()`](https://gillescolling.com/linop/reference/as_sparse.md)
  : Materialise an operator as a sparse matrix

## Capabilities and evidence

What an operator knows about itself, the argument that established it,
and what a solver will accept.

- [`cap()`](https://gillescolling.com/linop/reference/cap.md) : Read a
  capability off an operator
- [`capability()`](https://gillescolling.com/linop/reference/capability.md)
  : A capability value together with its evidence
- [`evidence()`](https://gillescolling.com/linop/reference/evidence.md)
  : Describe the argument behind a capability value
- [`evidence_satisfies()`](https://gillescolling.com/linop/reference/evidence_satisfies.md)
  : Does this evidence satisfy a requirement?
- [`requirement()`](https://gillescolling.com/linop/reference/requirement.md)
  : State what evidence a solver will accept

## Check a claim

The conformance suite an adapter author runs before shipping.

- [`verify()`](https://gillescolling.com/linop/reference/verify.md) :
  Check that a claim is valid

## Read a certificate

Each verb puts its certificate where its own return value has room. This
reads it off any of them.

- [`certificate()`](https://gillescolling.com/linop/reference/certificate.md)
  : Read the certificate off a result

## Solve a system

One verb over seven Krylov methods, choosing on declared capabilities
and shape, returning what a matrix solve returns with a certificate
attached.

- [`solve(`*`<linop>`*`)`](https://gillescolling.com/linop/reference/solve.linop.md)
  : Solve a linear system

## Eigenvalues and singular values

Restarted Lanczos with full reorthogonalisation, and the same loop on
the augmented operator. Each pair carries a backward error and a bound
on how far the value can be from an exact one.

- [`eigs()`](https://gillescolling.com/linop/reference/eigs.md) :
  Eigenvalues and eigenvectors of a hermitian operator
- [`svds()`](https://gillescolling.com/linop/reference/svds.md) :
  Singular values and vectors of an operator

## Preconditioners

One direction internally, and a contract each method enforces at the
point the preconditioner is accepted.

- [`preconditioner()`](https://gillescolling.com/linop/reference/preconditioner.md)
  : Build a preconditioner
- [`as_preconditioner()`](https://gillescolling.com/linop/reference/as_preconditioner.md)
  : Treat an approximate inverse as a preconditioner
- [`as_preconditioner_inverse()`](https://gillescolling.com/linop/reference/as_preconditioner_inverse.md)
  : Treat an operator P ~ A^-1 as a preconditioner
