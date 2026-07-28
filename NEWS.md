# linop 0.0.0.9000

First development release. The finite object model and the expression graph.

## Operators

* `linop()` builds an operator from a base matrix, a `Matrix`, a callback, or any
  class with a `linop()` method. `linop_eye()` and `linop_scaling()` cover the
  identity and diagonal operators.
* `%*%`, `crossprod()`, `tcrossprod()`, `+`, `-`, `*`, `/`, `t()`, `Conj()` and
  `adjoint()` build a lazy expression graph. Nothing is applied until a block
  arrives.
* Eleven node types: the `dense`, `sparse`, `fun`, `identity` and `diag` leaves,
  and the `transpose`, `adjoint`, `conjugate`, `scale`, `sum` and `product`
  composites. `linop_register_node()` adds more through the same door.
* Four apply modes composed through a Klein four-group, so transpose, adjoint
  and conjugate stay three distinct things.
* Operators are `"double"` or `"complex"` and no path downcasts. A real operator
  applied to a complex block promotes its result.

## Capabilities

* Ten capabilities, each three-valued. `NA` means unknown and is never read as
  `FALSE`.
* Every value carries `evidence()`: a source, a guarantee, a confidence, and the
  evidence it depends on. `requirement()` and `evidence_satisfies()` let a
  consumer state what it accepts as a predicate over the whole chain.
* Propagation rules for sums, scales, the three views and products, sound in the
  direction that matters: `TRUE` only where the conclusion follows, `NA`
  elsewhere.

## Solvers

* `solve()` solves `a x = b` by one of seven Krylov methods: CG, MINRES, GMRES,
  FGMRES, BiCGSTAB, LSQR and LSMR. `method = "auto"` chooses from the declared
  capabilities and records why; naming a method is the caller asserting their own
  declaration and applies no evidence minimum.
* Several right-hand sides run in lockstep rather than one after another, at one
  block apply per step.
* `eigs()` returns eigenpairs of a hermitian operator by Lanczos with full
  reorthogonalisation, locking and thick restarts, and takes `sigma` for the
  eigenvalues nearest a shift. It needs no adjoint.
* `svds()` returns singular triplets by Golub-Kahan-Lanczos with full
  reorthogonalisation, locking and thick restarts.
* Every result carries a certificate: what residual was reached, what backward
  error it implies, and what argument establishes it. `target identity` is
  `not_checked`, because a residual does not say which eigenvalue was found.
* The eigenpair certificate's `forward error` is a real bound,
  `min_j |theta - lambda_j| <= ||r|| / ||x||` by Weyl, and it records the
  hermitian declaration it rests on.
* Non-convergence returns a `fail` certificate naming what stopped the run. A
  contradicted declaration is what throws.

## Checking

* `verify()` runs an eleven-check conformance suite over an operator and returns
  a certificate. Declared capabilities are probed, and a contradiction fails.
* The `complex linearity` check tests `A(x + iy) = A(x) + i A(y)`. A real
  operator that discards the imaginary part of its input satisfies every other
  check, because the promotion in the apply path upcasts the real result it
  returned.
* `preconditioner()` carries a contract that each method enforces at the point
  the preconditioner is accepted.

## Inspection

* `print()` shows the expression tree, `explain()` adds shape and cost per node,
  `summary()` reports the capability table with evidence.
* `collapse()` replaces subexpressions that are cheaper as dense leaves, and
  reports what it did and the number of applications at which the memory pays
  for itself.
* `as.matrix()` and `as_sparse()` materialise on request, guarded on size.

## Notes

* No `Imports` and no compiled code. `Matrix` is suggested and is used for
  sparse leaves and `as_sparse()`.
* R 4.4.0 is the floor, set by `crossprod()` and `tcrossprod()` becoming S3
  generic.
