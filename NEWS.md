# linop 0.0.0.9000

First development release. The finite object model and the expression graph.

## Operators

* `linop()` builds an operator from a base matrix, a `Matrix`, a callback, or any
  class with a `linop()` method. `linop_eye()` and `linop_scaling()` cover the
  identity and diagonal operators.
* `%*%`, `crossprod()`, `tcrossprod()`, `+`, `-`, `*`, `/`, `t()`, `Conj()` and
  `adjoint()` build a lazy expression graph. Nothing is applied until a block
  arrives.
* Thirteen node types: the `dense`, `sparse`, `fun`, `identity`, `diag` and
  `jacobi` leaves, and the `transpose`, `adjoint`, `conjugate`, `scale`, `sum`,
  `product` and `section` composites.
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
* `eigs(method = "arnoldi")` returns eigenpairs of any square operator by
  Arnoldi with Krylov-Schur restarting and locking, and `method = "auto"` takes
  it for an operator not established hermitian, which `eigs()` used to refuse.
  A real operator keeps real arithmetic and returns its conjugate pairs
  together. `which` gains `"largest_real"`, `"smallest_real"`,
  `"largest_imaginary"` and `"smallest_imaginary"` (RSpectra's `"LR"`, `"SR"`,
  `"LI"`, `"SI"`), and `sigma` may be complex. Its certificate reports the
  backward error and leaves the forward error and orthogonality rows
  `not_checked`, since neither holds for a non-normal operator.
* A callback declared real is handed the real and imaginary parts of a complex
  block as one real block, so it promotes without having to accept complex
  input.
* `svds()` returns singular triplets by Golub-Kahan-Lanczos with full
  reorthogonalisation, locking and thick restarts.
* `eigs(method = "primme")` and `svds(method = "primme")` run the PRIMME
  library when the `linop.primme` engine package is installed, and
  `method = "auto"` takes it then, falling back to the reference otherwise.
  PRIMME restarts implicitly, so its basis stays at a size of its own choosing
  whatever the problem. The reported values are Rayleigh quotients measured on
  `A` and the certificate is the same whichever engine ran; its `dispatch`
  lists every engine considered and why each was taken or passed over.
  `options(linop.backend_priority = c(primme = -1))` prefers the reference.
* `linop_space("fourier", derivative, multiplier)` is an operator on functions
  on the circle, held as a banded operator on their Fourier coefficients: a
  constant-coefficient differential operator plus multiplication by a
  trigonometric polynomial. `eigs()` runs on it directly, chooses the
  truncation, and certifies against the operator on the circle with an exact
  truncation term; `finite_section()` truncates it by hand.
* `linop_space("chebyshev")` takes polynomial coefficients on every derivative
  term (`derivative = list(a_0, a_1, ...)` in Chebyshev coefficients); the
  adjoint and self-adjointness follow with the coefficients' derivatives in
  the boundary form. `method = "galerkin"` covers `-(k u')' + q u` with
  Neumann and Robin conditions as well as Dirichlet.
* `eigs()` finds every copy of a repeated eigenvalue among the `k` it returns:
  after the search it deflates what it found, searches again from a fresh
  direction, and swaps in any pair that ranks better. Lanczos, and Arnoldi on a
  hermitian operator. It costs about one more converged pair.
* The domain of a `linop_space("chebyshev")` operator is its boundary
  conditions, which may be mixed (`weights = `). `adjoint()` returns the
  adjoint operator with the adjoint domain, by integration by parts, and the
  operator records whether it is self-adjoint in L2. `adjoint()` on Fourier,
  Hermite and Laguerre operators returns the adjoint operator too.
  `eigs(method = "galerkin")` discretises the weak form of a self-adjoint
  second-order Dirichlet problem and returns Courant-Fischer upper bounds,
  lifted to Chebyshev coefficients.
* `+` and scalar `*` on function-space operators of one kind merge them into one
  operator, whose order and range space are the larger of the terms'; a
  Chebyshev term may be built without boundary conditions, and the sum carries
  the problem's.
* `linop_space("hermite")` and `linop_space("laguerre")`: operators on the line
  in Hermite functions and on the half line in Laguerre functions, banded
  because `x` is tridiagonal in both, with `-d^2/dx^2 + V` (`V` a polynomial)
  and the Laguerre operator `-(x u')' + x u / 4` respectively. Truncated and
  certified like the circle's.
* `linop_jacobi()` takes periodic limits: `diagonal_limit` and
  `offdiagonal_limit` may be one period of values, the essential spectrum is
  then a union of bands, and the certificate measures isolation against the
  nearest band. A point a finite section invents in a gap comes back with a
  failing truncation row, since its eigenvector sits at the cut; a genuine gap
  eigenvalue certifies.
* `linop_space("chebyshev", derivative, multiplier, boundary)` is a
  differential operator on `[-1, 1]` with polynomial coefficients and boundary
  conditions, held in the ultraspherical method: banded from Chebyshev
  coefficients to those of `C^(p)`. `solve()` takes the right-hand side's
  Chebyshev coefficients and `eigs(which = "smallest")` the smallest eigenvalues;
  both choose the number of coefficients and certify the exact residual on the
  untruncated system.
* `linop()` takes a `DelayedMatrix` (DelayedArray, and so HDF5-backed matrices
  from HDF5Array), a two-dimensional `torch_tensor` of type float64 or
  complex128, and a `big.matrix` (bigmemory, file-backed or not), applied by
  their own block-processed or on-device products, or for a `big.matrix` a
  block of rows at a time. A
  single-precision tensor is refused, since its apply error is below nothing
  the certificates can see.
* A dense real symmetric leaf of order up to 200 proves positive definiteness
  by a verified Cholesky factorisation (Rump's approach, Higham's backward
  error bound), so `solve(method = "auto")` reaches CG on it; until now nothing
  in the package reached that branch. A proof that fails leaves the capability
  unknown, never false.
* Structured operators: `linop_kron()`, `linop_blockdiag()`, `linop_lowrank()`,
  `linop_perm()` and `linop_power()`, `cbind()` and `rbind()` on linops, and
  `solve(a)` with no right-hand side for the exact inverse where one exists (a
  stored matrix, or a structure whose inverse is known). Each is a node, applied
  matrix-free, with capabilities propagated from its parts and its norm read off
  the structure where that is exact.
* `eigs(method = "auto")` acts on a hermitian capability only at computed,
  constructed or proven evidence, as `solve()` does for MINRES and CG: a bare
  `hermitian = TRUE` declaration on a callback now runs Arnoldi under auto, and
  Lanczos or PRIMME by name. For a pencil, `B`'s hermitian and positive definite
  evidence meets CG's minimum. The certificate still reads the declaration.
* Arnoldi on a hermitian operator returns real values in every storage mode,
  and Schur vectors as eigenvectors, so a repeated eigenvalue comes back with an
  orthonormal basis of its eigenspace.
* Fixed: `adjoint(X) %*% X` was stamped positive definite for a wide `X`
  declared `full_rank`, whose Gram matrix is singular. Definiteness now needs
  full rank on the side the Gram matrix lives on.
* `eigs(B = )` solves the generalized problem `A x = lambda B x`. For `A`
  hermitian and `B` hermitian positive definite it runs Lanczos on `B^-1 A` in
  the `B` inner product and returns `B`-orthonormal vectors, with a forward
  error bound through a solve with `B`; otherwise Arnoldi on `B^-1 A`. A shift
  runs on `(A - sigma B)^-1 B`. The backward error perturbs both operators,
  `||A x - theta B x|| / ((||A|| + |theta| ||B||) ||x||)`.
* `eigs(method = "rspectra")` and `svds(method = "rspectra")` run the Spectra
  library through RSpectra, when it is installed, on a real operator: its
  implicitly restarted Arnoldi, or Lanczos for an operator established
  hermitian, and its largest singular values. A complex operator is refused by
  name, because RSpectra would coerce it and answer for `Re(A)`. A shift runs
  on `eigs()`'s own `(A - sigma I)^-1`. `maxit` counts applies for it as for
  every engine. `method = "auto"` never takes it, since it stops on its own
  rule rather than the certificate's.
* Every result carries a certificate: what residual was reached, what backward
  error it implies, and what argument establishes it. `target identity` is
  `not_checked`, because a residual does not say which eigenvalue was found.
* `certificate()` reads it off whichever result carries it: an attribute on a
  solution, a field of a spectral result, or the certificate `verify()` returns
  directly. An operator has none, and the refusal names `verify()`.
* The eigenpair certificate's `forward error` is a real bound,
  `min_j |theta - lambda_j| <= ||r|| / ||x||` by Weyl, and it records the
  hermitian declaration it rests on.
* Non-convergence returns a `fail` certificate naming what stopped the run. A
  contradicted declaration is what throws.
* The reference `eigs()` and `svds()` end a run as stalled only when a round
  neither moves a wanted Ritz value beyond the arithmetic floor nor beats the
  best residual a wanted pair has held. A residual rises and falls across thick
  restarts while the Ritz value converges, and the smallest residual alone
  ended runs early: `eigs()` on the 1-D Laplacian of size 150, three largest
  at the defaults, stopped at 83 applies with an error of 6.6e-3 and now
  converges in 237.
* The default `maxit` of `eigs()` and `svds()` is `20 n`, with `n` the
  dimension the recurrence runs in (`min(m, n)` for `svds()`), in place of
  `max(300, 20 k)`. At the default `ncv` the applies a crowded end of the
  spectrum needs grow with `n`: the three largest eigenvalues of the 1-D
  Laplacian of size 200 take 322, the largest at size 2000 takes 14255. The
  rule covers the worst converging run of a sweep over Laplacians, random and
  clustered spectra and difference operators by 1.6x for the reference engine
  and 1.7x for PRIMME. The small end of a spectrum spread over four decades
  needs up to 100 `n` and is left to `sigma`
  (`dev_notes/spikes/eigs-budget-sweep.R`).
* The reference `eigs()` and `svds()` build the projected matrix from the
  coefficients their reorthogonalisation removes, `V^H A V` and `U^H A V` as
  measured, instead of the tridiagonal and arrow form exact arithmetic gives
  after a thick restart. The retained Ritz vectors are formed with a rounding
  error at every restart, which the assumed form took as exact; at the small
  end of a crowded spectrum that error decided convergence. The smallest
  eigenvalues of the 1-D Laplacian of size 500 now converge for k = 1, 3, 5
  and 10 (0, 2, 4 and 9 pairs did), the residual the run plateaus at for
  size 1000 falls from 20.5 to 13.6 eps ||A||, and runs that converged
  before take the same number of applies.
* The reference `svds()` bidiagonalises `A^H` when `A` is wider than tall.
  Started in a domain larger than the codomain, the run picked up the start
  vector's component in `null(A)`, which is not a singular direction, and
  returned it as a zero singular value with no left vector to match:
  `svds(diff_1d(100), k = 1, which = "smallest")` stopped with `d = 2.7e-17`
  and a `fail` certificate where the smallest singular value is `3.1e-2`. A
  caller's `v0` is still in the domain of `A`.
* `svds()` through PRIMME takes PRIMME's estimate of `||A||` into the lower
  bound its convergence test reads, as `eigs()` through PRIMME already did. The
  bound from the measured triplets alone could sit well below `||A||` at the
  small end, putting the residual row's floor below what a converged triplet
  reaches: a 60 x 60 operator with smallest singular value 1.67e-8 ran to the
  budget with the right triplet in hand and now converges in 651 applies.
* The `enclosure` row of a finite section carries the arithmetic floor at the
  Temple end as well as at the variational one. Both ends are inequalities
  about the computed Rayleigh quotient, and at `n = 40` the Temple half-width
  (2.8e-17) is below its rounding.
* A reference run whose budget runs out in a round that locked a pair reports
  each pair once. The locked pair was also taken again from that round's
  unconverged pairs, returning a duplicate in place of the next one.

## Operators on a sequence space

* A dimension is two non-negative whole numbers or `Inf`. `Inf` is the extent of
  the index set rather than a flag, so the conformability test in a product, the
  squareness test for a hermitian operator and `rev()` for a transpose all keep
  working with no special case.
* `linop_jacobi()` builds a self-adjoint Jacobi operator on `l^2(Z)` from its
  diagonal and off-diagonal sequences, both real and eventually constant. Its
  essential spectrum is exact by Weyl rather than estimated.
* `finite_section()` truncates it to an ordinary operator, as a node holding the
  operator it truncates, so the printed tree shows the infinite one underneath
  the finite one.
* `decay_rate()` returns the rate an eigenvector decays at outside the window,
  closed form in the eigenvalue alone, and `NA` inside the essential spectrum
  where nothing decays.
* `eigs()` on a finite section certifies a statement about the operator that was
  truncated. `forward error` bounds the distance from the computed value to the
  spectrum of that operator, and it is the first row in the package to rest on
  nothing declared. `isolation` refuses a value that has not cleared the
  essential spectrum, which is what separates an eigenvalue from a
  discretisation of continuous spectrum.
* `eigs()` also runs on the operator itself and chooses the width. The
  eigenvector is exactly geometric outside the window, so the truncation term at
  one width gives it at every other and the width that meets a tolerance is a
  division rather than a search: one section to measure and one to confirm,
  whatever the rate. The widths tried are on the result and the reason the
  search stopped is on the certificate.
* Every other numeric verb refuses an infinite operator by name and points at
  `finite_section()`.

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
  sparse leaves and `as_sparse()`. The PRIMME engine is compiled and lives in
  the separate `linop.primme` package, so installing `linop` compiles nothing.
* R 4.4.0 is the floor, set by `crossprod()` and `tcrossprod()` becoming S3
  generic.
