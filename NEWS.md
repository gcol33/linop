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

## Checking

* `verify()` runs a ten-check conformance suite over an operator and returns a
  certificate. Declared capabilities are probed, and a contradiction fails.
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
