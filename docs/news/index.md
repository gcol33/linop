# Changelog

## linop 0.0.0.9000

First development release. The finite object model and the expression
graph.

### Operators

- [`linop()`](https://gillescolling.com/linop/reference/linop.md) builds
  an operator from a base matrix, a `Matrix`, a callback, or any class
  with a [`linop()`](https://gillescolling.com/linop/reference/linop.md)
  method.
  [`linop_eye()`](https://gillescolling.com/linop/reference/linop_eye.md)
  and
  [`linop_scaling()`](https://gillescolling.com/linop/reference/linop_scaling.md)
  cover the identity and diagonal operators.
- `%*%`, [`crossprod()`](https://rdrr.io/r/base/crossprod.html),
  [`tcrossprod()`](https://rdrr.io/r/base/crossprod.html), `+`, `-`,
  `*`, `/`, [`t()`](https://rdrr.io/r/base/t.html),
  [`Conj()`](https://rdrr.io/r/base/complex.html) and
  [`adjoint()`](https://gillescolling.com/linop/reference/adjoint.md)
  build a lazy expression graph. Nothing is applied until a block
  arrives.
- Thirteen node types: the `dense`, `sparse`, `fun`, `identity`, `diag`
  and `jacobi` leaves, and the `transpose`, `adjoint`, `conjugate`,
  `scale`, `sum`, `product` and `section` composites.
- Four apply modes composed through a Klein four-group, so transpose,
  adjoint and conjugate stay three distinct things.
- Operators are `"double"` or `"complex"` and no path downcasts. A real
  operator applied to a complex block promotes its result.

### Capabilities

- Ten capabilities, each three-valued. `NA` means unknown and is never
  read as `FALSE`.
- Every value carries
  [`evidence()`](https://gillescolling.com/linop/reference/evidence.md):
  a source, a guarantee, a confidence, and the evidence it depends on.
  [`requirement()`](https://gillescolling.com/linop/reference/requirement.md)
  and
  [`evidence_satisfies()`](https://gillescolling.com/linop/reference/evidence_satisfies.md)
  let a consumer state what it accepts as a predicate over the whole
  chain.
- Propagation rules for sums, scales, the three views and products,
  sound in the direction that matters: `TRUE` only where the conclusion
  follows, `NA` elsewhere.

### Solvers

- [`solve()`](https://rdrr.io/r/base/solve.html) solves `a x = b` by one
  of seven Krylov methods: CG, MINRES, GMRES, FGMRES, BiCGSTAB, LSQR and
  LSMR. `method = "auto"` chooses from the declared capabilities and
  records why; naming a method is the caller asserting their own
  declaration and applies no evidence minimum.
- Several right-hand sides run in lockstep rather than one after
  another, at one block apply per step.
- [`eigs()`](https://gillescolling.com/linop/reference/eigs.md) returns
  eigenpairs of a hermitian operator by Lanczos with full
  reorthogonalisation, locking and thick restarts, and takes `sigma` for
  the eigenvalues nearest a shift. It needs no adjoint.
- [`svds()`](https://gillescolling.com/linop/reference/svds.md) returns
  singular triplets by Golub-Kahan-Lanczos with full
  reorthogonalisation, locking and thick restarts.
- `eigs(method = "primme")` and `svds(method = "primme")` run the PRIMME
  library when the `linop.primme` engine package is installed, and
  `method = "auto"` takes it then, falling back to the reference
  otherwise. PRIMME restarts implicitly, so its basis stays at a size of
  its own choosing whatever the problem. The reported values are
  Rayleigh quotients measured on `A` and the certificate is the same
  whichever engine ran; its `dispatch` lists every engine considered and
  why each was taken or passed over.
  `options(linop.backend_priority = c(primme = -1))` prefers the
  reference.
- Every result carries a certificate: what residual was reached, what
  backward error it implies, and what argument establishes it.
  `target identity` is `not_checked`, because a residual does not say
  which eigenvalue was found.
- [`certificate()`](https://gillescolling.com/linop/reference/certificate.md)
  reads it off whichever result carries it: an attribute on a solution,
  a field of a spectral result, or the certificate
  [`verify()`](https://gillescolling.com/linop/reference/verify.md)
  returns directly. An operator has none, and the refusal names
  [`verify()`](https://gillescolling.com/linop/reference/verify.md).
- The eigenpair certificate’s `forward error` is a real bound,
  `min_j |theta - lambda_j| <= ||r|| / ||x||` by Weyl, and it records
  the hermitian declaration it rests on.
- Non-convergence returns a `fail` certificate naming what stopped the
  run. A contradicted declaration is what throws.
- The reference
  [`eigs()`](https://gillescolling.com/linop/reference/eigs.md) and
  [`svds()`](https://gillescolling.com/linop/reference/svds.md) end a
  run as stalled only when a round neither moves a wanted Ritz value
  beyond the arithmetic floor nor beats the best residual a wanted pair
  has held. A residual rises and falls across thick restarts while the
  Ritz value converges, and the smallest residual alone ended runs
  early: [`eigs()`](https://gillescolling.com/linop/reference/eigs.md)
  on the 1-D Laplacian of size 150, three largest at the defaults,
  stopped at 83 applies with an error of 6.6e-3 and now converges in
  237.
- The default `maxit` of
  [`eigs()`](https://gillescolling.com/linop/reference/eigs.md) and
  [`svds()`](https://gillescolling.com/linop/reference/svds.md) is
  `20 n`, with `n` the dimension the recurrence runs in (`min(m, n)` for
  [`svds()`](https://gillescolling.com/linop/reference/svds.md)), in
  place of `max(300, 20 k)`. At the default `ncv` the applies a crowded
  end of the spectrum needs grow with `n`: the three largest eigenvalues
  of the 1-D Laplacian of size 200 take 322, the largest at size 2000
  takes 14255. The rule covers the worst converging run of a sweep over
  Laplacians, random and clustered spectra and difference operators by
  2.5x for the reference engine and 1.7x for PRIMME. The small end of a
  spectrum spread over four decades needs up to 100 `n` and is left to
  `sigma` (`dev_notes/spikes/eigs-budget-sweep.R`).
- A reference run whose budget runs out in a round that locked a pair
  reports each pair once. The locked pair was also taken again from that
  round’s unconverged pairs, returning a duplicate in place of the next
  one.

### Operators on a sequence space

- A dimension is two non-negative whole numbers or `Inf`. `Inf` is the
  extent of the index set rather than a flag, so the conformability test
  in a product, the squareness test for a hermitian operator and
  [`rev()`](https://rdrr.io/r/base/rev.html) for a transpose all keep
  working with no special case.
- [`linop_jacobi()`](https://gillescolling.com/linop/reference/linop_jacobi.md)
  builds a self-adjoint Jacobi operator on `l^2(Z)` from its diagonal
  and off-diagonal sequences, both real and eventually constant. Its
  essential spectrum is exact by Weyl rather than estimated.
- [`finite_section()`](https://gillescolling.com/linop/reference/finite_section.md)
  truncates it to an ordinary operator, as a node holding the operator
  it truncates, so the printed tree shows the infinite one underneath
  the finite one.
- [`decay_rate()`](https://gillescolling.com/linop/reference/decay_rate.md)
  returns the rate an eigenvector decays at outside the window, closed
  form in the eigenvalue alone, and `NA` inside the essential spectrum
  where nothing decays.
- [`eigs()`](https://gillescolling.com/linop/reference/eigs.md) on a
  finite section certifies a statement about the operator that was
  truncated. `forward error` bounds the distance from the computed value
  to the spectrum of that operator, and it is the first row in the
  package to rest on nothing declared. `isolation` refuses a value that
  has not cleared the essential spectrum, which is what separates an
  eigenvalue from a discretisation of continuous spectrum.
- [`eigs()`](https://gillescolling.com/linop/reference/eigs.md) also
  runs on the operator itself and chooses the width. The eigenvector is
  exactly geometric outside the window, so the truncation term at one
  width gives it at every other and the width that meets a tolerance is
  a division rather than a search: one section to measure and one to
  confirm, whatever the rate. The widths tried are on the result and the
  reason the search stopped is on the certificate.
- Every other numeric verb refuses an infinite operator by name and
  points at
  [`finite_section()`](https://gillescolling.com/linop/reference/finite_section.md).

### Checking

- [`verify()`](https://gillescolling.com/linop/reference/verify.md) runs
  an eleven-check conformance suite over an operator and returns a
  certificate. Declared capabilities are probed, and a contradiction
  fails.
- The `complex linearity` check tests `A(x + iy) = A(x) + i A(y)`. A
  real operator that discards the imaginary part of its input satisfies
  every other check, because the promotion in the apply path upcasts the
  real result it returned.
- [`preconditioner()`](https://gillescolling.com/linop/reference/preconditioner.md)
  carries a contract that each method enforces at the point the
  preconditioner is accepted.

### Inspection

- [`print()`](https://rdrr.io/r/base/print.html) shows the expression
  tree,
  [`explain()`](https://gillescolling.com/linop/reference/explain.md)
  adds shape and cost per node,
  [`summary()`](https://rdrr.io/r/base/summary.html) reports the
  capability table with evidence.
- [`collapse()`](https://gillescolling.com/linop/reference/collapse.md)
  replaces subexpressions that are cheaper as dense leaves, and reports
  what it did and the number of applications at which the memory pays
  for itself.
- [`as.matrix()`](https://rdrr.io/r/base/matrix.html) and
  [`as_sparse()`](https://gillescolling.com/linop/reference/as_sparse.md)
  materialise on request, guarded on size.

### Notes

- No `Imports` and no compiled code. `Matrix` is suggested and is used
  for sparse leaves and
  [`as_sparse()`](https://gillescolling.com/linop/reference/as_sparse.md).
  The PRIMME engine is compiled and lives in the separate `linop.primme`
  package, so installing `linop` compiles nothing.
- R 4.4.0 is the floor, set by
  [`crossprod()`](https://rdrr.io/r/base/crossprod.html) and
  [`tcrossprod()`](https://rdrr.io/r/base/crossprod.html) becoming S3
  generic.
