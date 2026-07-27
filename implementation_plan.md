# linop: implementation plan

Status: draft
Date: 2026-07-27
Target: CRAN, plus a JOSS or MEE paper once there is demonstrated research use

---

## 1. What linop is

An interoperability contract for matrix-free numerical linear algebra in R, in three
layers.

**Layer 1, the finite contract.** A `linop` behaves like a matrix under R's normal
arithmetic, carries evidence-bearing capabilities, composes lazily into a typed expression
graph, and is consumed by iterative solvers, least-squares solvers, eigensolvers and
singular-value solvers without ever being materialised.

**Layer 2, production engines.** Vendored numerical libraries behind a backend registry
with a deterministic selection protocol.

**Layer 3, Hilbert-space operators.** Operators between infinite-dimensional spaces,
never solved directly but discretised into Layer 1 objects under an explicit
discretisation, with discretisation error certified alongside algebraic error.

The near-term defensible asset is Layer 1: a rigorously tested interoperability contract
that another package can adopt in its first year. Layer 3 has the higher intellectual
ceiling and demonstrates that the contract is the numerical target of a discretisation
system rather than a wrapper around finite matrices.

Non-goals for v1: dense direct factorisations, distributed operators, GPU backends, PDE
geometry beyond one dimension.

---

## 1.1 The public surface

Everything in sections 4 through 8 is machinery. This is what a user types.

> **A `linop` behaves like a matrix. Use ordinary matrix operations and ordinary numerical
> verbs.**

```r
A <- linop(M)

A %*% x
t(A); adjoint(A); Conj(A)
A + B;  A %*% B;  2 * A;  crossprod(A)

x   <- solve(A, b)
fit <- eigs(A, k = 10)
sv  <- svds(A, k = 10)
verify(A)
```

Six visible names, plus R's own generics. Everything else is progressive disclosure.

### `linop()` is one generic, dispatching on its first argument

```r
linop(M)                                        # matrix, Matrix, sparse
linop(function(X) ..., adjoint = , dim = )      # callback
linop(x)                                        # third party writes linop.myclass()
```

No separate `as_linop()`. Dispatching on the first argument covers the adapter case and
the callback case in one verb, since a callback constructor whose first argument *is* the
apply function reads the same way an adapter does. Third-party packages implement
`linop.myclass()` and nothing else.

Declared properties arrive as a named logical, not a character vector:

```r
linop(f, adjoint = g, dim = c(n, n),
      properties = c(hermitian = TRUE, positive_definite = TRUE))
```

Named logical rather than `c("hermitian", "positive_definite")` because the capability
system is three-valued everywhere else, and a character vector cannot say `FALSE`.
Declaring an operator explicitly non-hermitian is useful information that saves a probe. A
bare character vector is accepted as shorthand for all-`TRUE`. The constructor stamps
`source = "user_declaration"` internally; users never build evidence objects.

### `%*%` returns what a matrix would return

```r
A %*% x    # x numeric  -> numeric.   Computed now.
A %*% B    # B a linop  -> linop.     Lazy.
```

Dispatch is on the second argument. Returning a lazy node where a matrix would have
returned a vector is the single fastest way to break the mental model, so it does not
happen.

### `solve()`, not `solve_iterative()`

```r
x <- solve(A, b)
x <- solve(A, b, method = "auto", preconditioner = "jacobi", tol = 1e-8)
```

If the operator declares a direct solve capability, `solve()` uses it. Otherwise it solves
iteratively. Either way the result records which happened and what residual was achieved,
reachable through `summary(x)`.

This reverses an earlier decision in this plan, which had `solve()` error without a direct
capability and point at `solve_iterative()`. Refusing was the wrong way to honour "nothing
silent": reporting achieves it without spending a verb, and for a `linop` there is usually
no factorization available, so iterative is the only meaning `solve()` can have.

`solve()` returns a plain numeric vector or matrix with the certificate in an attribute;
`solve(..., details = TRUE)` returns the full object. The tradeoff is real and documented:
arithmetic on the result drops the attribute. A classed object would survive arithmetic
and would break the "behaves like a matrix" promise in exactly the ways R users dislike,
so the attribute is the better side to err on.

**Rectangular input.** Base R's `solve()` requires square, and least squares is a
different mathematical request. `solve()` on a rectangular `linop` errors and names the
least-squares methods, which are then reachable through the same verb:

```r
solve(A, b, method = "lsmr")
```

That keeps base semantics, keeps the verb budget, and makes the user state their intent
exactly where the meaning changes.

### `eigs()` and `svds()` take strings

```r
eigs(A, k = 5, which = "largest")
eigs(A, k = 5, which = "smallest")
eigs(A, k = 5, sigma = 0)          # shift-invert
eigs(A, k = 5, B = M)              # generalized
svds(A, k = 10, which = "smallest")
```

Not `target = largest()` or `method = shift_invert(...)`. Constructor objects are
architecturally tidy and impose a miniature language on the user.

`which` also accepts RSpectra's short forms (`"LM"`, `"SM"`, and so on), because the
incumbent's vocabulary is what users already have in their fingers and lowering the
switching cost matters more here than vocabulary purity.

### `verify()` is one generic

```r
verify(A)         # does this operator satisfy the contract?
verify(fit, A)    # is this result valid for this operator?
```

Internally these run different checks. Publicly they answer one question: is this object's
claim valid. No separate `linop_verify()` and `certify()`.

### No public planner in v0.1

The solver chooses and reports:

```r
fit <- eigs(A, k = 5)
summary(fit)
```

```
method:  block Lanczos
backend: primme
reason:  hermitian operator, smallest eigenvalues requested
```

An expert overrides directly, `eigs(A, k = 5, method = "lanczos", backend = "primme")`. A
separate plan object adds a workflow step without solving a user problem. `plan_eigs()`
remains as a developer and debugging tool documented with the backend protocol.

### The three classes stay internal

`linop` / `linsolve` / `preconditioner` is mathematically essential and is not three
concepts every user learns.

```r
solve(A, b, preconditioner = "jacobi")     # common
P <- preconditioner(apply_inverse = ..., fixed = TRUE, ...)   # advanced, public
solve(A, b, preconditioner = P)
```

`solver()` is **not public in v0.1**. The case that would justify it is a sequence of
right-hand sides against one operator with an expensive preconditioner setup, and that is
already covered: the reusable object is the preconditioner, not the solver. Multiple
simultaneous right-hand sides go through the block interface as a matrix `b`. `solver()`
becomes public when inexact shift-invert or a PRIMME warm-start workflow creates a real
need. Designing the class now does not oblige exposing it now.

This defer is only safe because `preconditioner()` is public and reusable. If that
changes, `solver()` comes forward with it.

### Hilbert reuses the same verbs

```r
H   <- hilbert_operator(...)
fit <- eigs(H, k = 3, scheme = finite_section())
```

The `eigs()` method for a Hilbert operator checks whether the request is justified,
discretises, calls finite `eigs()`, lifts, and certifies against the original. Unjustified
requests error and name the weaker verb:

```r
spectral_approximation(H, scheme, n)
```

That verb exists because it returns a fundamentally weaker claim, which is the one reason
in this plan that justifies a parallel name. Everything else shares the finite API.

### API budget

```
common      linop  adjoint  solve  eigs  svds  verify
generics    %*%  +  -  *  t  Conj  crossprod  tcrossprod  dim  nrow  ncol
            solve  print  summary  as.matrix
advanced    preconditioner  spectral_approximation  collapse  explain
developer   linop_register_node  linop_register_backend  plan_eigs
            provenance_*  evidence_satisfies  capability  evidence
```

Developer names live in separate documentation under "Writing an adapter" and "Writing a
backend", and carry an explicit stability contract. The README never opens with a
registration call.

### What the README teaches

```r
A <- linop(function(X) ..., adjoint = function(X) ..., dim = c(n, n))

x <- solve(A, b)
eigs(A, k = 5)
svds(A, k = 5)
verify(A)
```

then that the algebra is lazy:

```r
B <- crossprod(A) + 0.1 * linop(diag(n))
```

That is the package. The evidence system, materialisation rules, solve contracts,
provenance envelopes, backend protocol and Hilbert-domain theory are what make those lines
trustworthy. They are not what anyone should have to type.

Design goal: **six functions to use the package, substantial machinery to make those six
functions honest.**

---

## 2. Verified landscape

Checked 2026-07-27. Facts below were read from the source given; anything unchecked is
marked.

### Finite side

| Fact | Source |
|---|---|
| `linop` is free as a CRAN package name (404) | cran.r-project.org/web/packages/linop |
| eigencore 1.0.2, 2026-07-25, Bradley Buchsbaum, `Imports: Matrix, methods`, no `LinkingTo` | CRAN eigencore |
| The `PRIMME` R package was removed from CRAN 2026-06-07, archived because email to the maintainer is undeliverable | CRAN PRIMME |
| PRIMME 3.2-6 checked **OK on all three Windows platforms**. NOTEs elsewhere: C library writing to stdout, unused `Matrix` import | CRAN archived checks, 2026-06-07 |
| PRIMME C library is BSD 3-clause, Make-based, needs BLAS/LAPACK, repo active | github.com/primme/primme |
| `RSpectra::eigs()` accepts `function(x, args)` computing `Ax`, single vector | RSpectra manual |
| `Rlinsolve` 0.3.3, 2025-09-22, Kisung You. Twelve solvers. No MINRES, no LSQR, no LSMR | Rlinsolve index |

### Hilbert side

| Fact | Source |
|---|---|
| ApproxFun operators carry `domainspace` and `rangespace`; `op * f` lands in `rangespace(op)` | ApproxFun docs |
| `Derivative(Chebyshev())` lands in `Ultraspherical(1)`; mismatches repaired by automatic conversion, "if one writes `D + I` it will translate it to `D + C`" | ApproxFun docs |
| Spectral pollution can occur at **any** point in a gap of the extended essential spectrum of a self-adjoint operator | Davies & Plum, second-order relative spectra literature |
| For self-adjoint `A` not bounded above, perturbing a finite basis by arbitrarily small vectors in `D(A)` pushes the finite-section spectrum arbitrarily far right | same |
| No Chebfun or ApproxFun equivalent on CRAN | search, provisional pending a Numerical Mathematics task view read |

**Not verified.** Every claim about eigencore's internals in the strategy documents is
second-hand. Spike S0.4 verifies or refutes each. None may appear in a README, vignette,
benchmark note, issue or paper until then.

### Where the gaps are

Finite: matrix-free eigensolving is not empty, since RSpectra takes a function. Missing is
a composable, typed, capability-carrying operator shared across solver families. Also
missing: operator-native Krylov; MINRES, LSQR and LSMR anywhere; block matvecs in any
matrix-free interface; a preconditioner that is the same kind of thing as what it
preconditions; any conformance suite a third party can run.

Infinite: nothing in R at all.

---

## 3. Architecture: three packages

**`linop`** (core). Pure R. `Imports: methods`, nothing else. No `LinkingTo`, no compiled
code in v0.1. Holds the finite object model, expression graph, algebra, capability and
evidence lattices, adapters, conformance suite, certificates, the provenance envelope, and
all pure-R Krylov solvers. Core does not know Hilbert spaces exist.

**`linop.primme`** (engine backend). Vendored PRIMME C source plus bindings.

**`linop.hilbert`** (Hilbert layer). Spaces, operators, forms, discretisations, lifting,
discretisation certificates. Depends on `linop`; nothing depends on it.

### Why the splits are not optional

CRAN policy, verbatim:

> "Orphaned CRAN packages should not be strict requirements (in the 'Depends', 'Imports'
> or 'LinkingTo' fields, including indirectly)."

That closes depending on a revived-but-orphaned `PRIMME`, so the engine must be vendored,
which means compiled code and a BLAS/LAPACK link. A package whose proposition is "depend
on this abstraction" cannot make dependents compile PRIMME first.

The Hilbert split has a different and equally hard reason. Its cadence is research-paced;
adaptive estimators, finite-section bounds and pollution-free methods are open problems in
places. Core's cadence must be boring, because adopters build on its contract. Coupling a
stable interface to a moving research layer damages the interface.

Rule: **anything that costs an install second, or moves at research pace, goes in a
separate package.**

### Backend protocol

`supports() == TRUE` is insufficient once two backends are installed. Each backend
declares:

```r
linop_register_backend(
  name              = "primme",
  contract_version  = "1.0",
  engine_version    = "3.2.6",
  provides          = c("eigs_hermitian", "eigs_gen", "svds", "eigs_generalized"),
  dtypes            = c("double", "complex"),
  requires          = c("apply", "adjoint"),      # operator capabilities needed
  targets           = c("LM", "SM", "interior"),
  preconditioning   = TRUE,
  warm_start        = TRUE,
  priority          = 100L,
  supports = function(A, k, which, opts) TRUE | "reason for rejection"
)
```

`supports()` returns `TRUE` or a rejection string, never a bare `FALSE`, so every
exclusion is explainable.

Selection order, deterministic and reported:

1. explicit `backend =` argument
2. user-configured priority (option)
3. solver-family default priority
4. declared applicability

No environment-dependent benchmark guessing by default, because reproducibility matters
more than the last few percent.

Selection is reported through `summary(fit)` in normal use. `plan_eigs(A, k = 10)` is a
developer and debugging tool, documented with this protocol rather than in the user
vignettes:

```
primme      eligible    selected
rspectra    eligible    lower configured priority
reference   eligible    reference implementation
torch       rejected    operator is on CPU
```

---

## 4. Three kinds of object

The most important structural decision in the finite layer, and the one an earlier draft
got wrong by treating an iterative inverse as an operator.

```
linop            applies a genuinely linear map
linsolve         computes an approximation to A^{-1} b; may be adaptive
preconditioner   transforms residuals; may be fixed or flexible
```

They share capability and evidence conventions. They are not the same mathematical object
and must not share a class.

### 4.1 The membership test

> An object is a `linop` if and only if its action on `x` is determined without reference
> to `x`.

Equivalently: if the object applies `p(A)` for a polynomial or rational `p`, then `p` must
be fixed before any vector is seen.

| Object | `linop`? | Why |
|---|---|---|
| `precond_polynomial(A, degree, bounds)` | yes | polynomial fixed from spectral bounds a priori |
| Chebyshev semi-iteration, fixed degree | yes | fixed polynomial |
| Neumann series to fixed degree | yes | fixed |
| `k` sweeps of Jacobi, SSOR, IC(0), ILU(0) | yes | fixed splitting, fixed sweep count |
| `k` steps of CG from `x0 = 0` | **no** | see below |
| CG or GMRES to a tolerance | no | iteration count varies with the right-hand side |
| Any restarted or adaptive scheme | no | |

**Why fixed-iteration CG is excluded, and it is not a finite-precision issue.** CG chooses
its polynomial by minimising the A-norm error over `K_k(A, b)`. The minimiser depends on
which eigencomponents `b` has: one step is exact when `b` is an eigenvector and not
otherwise. So the polynomial is a function of `b`, and `S(b1 + b2) != S(b1) + S(b2)` in
exact arithmetic. Chebyshev semi-iteration looks similar and *is* linear, because its
polynomial comes from spectral bounds rather than from `b`. That contrast is the reason
the membership test is stated as a rule rather than as a list.

### 4.2 `linsolve`

```r
S <- linsolve(A, method = "cg", tol = 1e-10, M = , maxit = )
S(b)                      # returns a solution plus a certificate
```

`solve(A, b, ...)` is the one-shot form of the same thing, and is the only form the public
API exposes in v0.1 (section 1.1).

Shift-invert holds a solve object, not an operator. Internally:

```r
S <- linsolve(A - sigma * linop_identity(n), method = "minres", tol = 1e-12)
```

Publicly this is `eigs(A, k, sigma = 0)`, which constructs the solve object itself. The
object surfaces only when a user needs to control the inner solve, which is the point at
which `solver()` becomes public.

**A `linsolve` declares its contract, and compatibility is enforced rather than
documented.** Ordinary Arnoldi must not receive a history-dependent solve merely because
both objects implement `solve()`.

```
fidelity:     exact | linear_approximation | variable_inexact
determinacy:  fixed | input_dependent | history_dependent
randomness:   deterministic | stochastic
adjoint:      available | absent
error_bound:  available | absent
```

| Solve object | fidelity | determinacy |
|---|---|---|
| sparse LU | exact | fixed |
| fixed polynomial inverse | linear_approximation | fixed |
| fixed zero-start stationary iteration | linear_approximation | fixed |
| CG to a tolerance | variable_inexact | input_dependent |
| CG warm-started from previous calls | variable_inexact | history_dependent |
| stochastic inner solver | variable_inexact | input_dependent, stochastic |

Backends and eigensolvers declare which contracts they accept:

```
accepted: exact
          fixed_linear
          variable_inexact with outer-controlled tolerance
```

The third generally requires inexact or flexible rational Krylov, Jacobi-Davidson, or a
relaxation strategy, not ordinary shift-invert Lanczos. `error_bound = available` is what
unlocks relaxation, where the inner tolerance is loosened as the outer iteration
converges; without an inner bound the outer method has nothing to relax against.

**Consequence for conformance.** A `history_dependent` solve returns different output for
the same input by design, so it fails the purity check in section 5.10. `linsolve` objects
therefore get their own conformance suite, which tests the declared contract rather than
assuming purity. Running the `linop` suite against a `linsolve` is itself an error.

### 4.3 `preconditioner`

Two representations are common and mixing them is a known bug source: `M ~ A`, applied by
solving `M z = r`, and `P ~ A^{-1}`, applied by computing `z = P r`. Internally only one
direction ever exists.

```r
preconditioner(
  apply_inverse     = function(R) ...,   # always z = M^{-1} r
  fixed             = TRUE,
  hermitian         = ,
  positive_definite = ,
  side              = "left" | "right" | "split"
)

as_preconditioner(M)          # from an operator M ~ A: applies solve(M, r)
as_preconditioner_inverse(P)  # from an operator P ~ A^{-1}: applies P %*% r
```

Solver requirements, enforced:

| Solver | Requires | Side |
|---|---|---|
| CG | `fixed`, `hermitian`, `positive_definite` | unrestricted |
| MINRES | `fixed`, `hermitian`, `positive_definite` | unrestricted |
| GMRES | `fixed` | unrestricted |
| FGMRES | `fixed` may be `FALSE` | `right` only |
| BiCGSTAB | `fixed` | unrestricted |
| LSQR, LSMR | `fixed`, acting on the normal equations | unrestricted |

`side` is enforced rather than merely recorded. FGMRES forms its update from the
right-preconditioned basis, so `right` is the only side it is defined for, and PETSc's
`KSPFGMRES` states the same restriction. The rows marked unrestricted admit more than one
standard formulation; each needs checking against Saad before it is narrowed, and until
then the table accepts all three sides, which is the conservative reading rather than a
finding.

A `fixed = FALSE` preconditioner passed to CG is an error naming the flag, not a warning.
A `preconditioner` whose `apply_inverse` is a fixed linear map can be converted to a
`linop` with `linop()`; the reverse conversion is always available.

---

## 5. The linop contract

### 5.1 Object model

```r
structure(
  list(
    node       = "product",
    dim        = c(m, n),
    dtype      = "double" | "complex",
    caps       = <capability list, each with value + evidence>,
    args       = list(...),
    cost       = <flops per column, estimated>,
    provenance = NULL
  ),
  class = c("linop_product", "linop")
)
```

Immutable. Every operation returns a new object.

### 5.2 Two-tier apply

Authors implement the easy form; solvers call the fast form.

```r
A <- linop(
  apply   = function(X) M %*% X,          # n x b in, m x b out
  adjoint = function(X) crossprod(M, X),
  dim     = c(m, n)
)
```

`X` is always a matrix; a vector is `b = 1`. Block-first from the start, since every block
Krylov method needs it and retrofitting blocks is a breaking change.

Optional tier 2: `apply_gemm(X, Y, alpha, beta)` computing `Y <- alpha*A%*%X + beta*Y`,
synthesised from `apply` when absent. Solvers only ever call tier 2, so scratch reuse is
uniform. The GEMM signature is the right internal ABI and the wrong thing to demand from
someone wrapping a data format.

### 5.3 Capabilities carry evidence

Random probes can disprove a false claim. They cannot prove one. So a capability is a
value plus a structured description of the argument behind it.

```r
caps$hermitian <- capability(
  value    = TRUE,
  evidence = evidence(source = "construction", guarantee = "identity")
)
```

Values are three-valued: `TRUE`, `FALSE`, `NA`. `NA` means unknown and is never read as
`FALSE`.

**Evidence is not a single strength number.** A ranking such as
`heuristic < estimate < computed < rigorous < exact` would be wrong, because the classes
are not comparable: a probabilistic bound can carry a rigorous confidence statement while
a deterministic estimate carries no guarantee at all; an adapter contract can be stronger
than a hundred random probes; construction proves hermitian structure and says nothing
about positive definiteness. Three independent fields instead:

```
source:      construction | adapter_contract | user_declaration | probe
             | computation | theorem
guarantee:   identity | deterministic_bound | probabilistic_bound | estimate | heuristic
confidence:  1 | a stated probability | NA
```

Solvers state a requirement as a set of acceptable combinations, and satisfaction is a
predicate, never an integer comparison:

```r
requirement(sources    = c("construction", "adapter_contract", "theorem"),
            guarantees = "identity",
            min_confidence = 1)

evidence_satisfies(actual, required)
```

CG's default requirement for `positive_definite` accepts construction, a trusted adapter
contract, or a theorem-backed declaration, and does not accept ten positive Rayleigh
quotients. The user can override explicitly; the override is recorded in the result.

**Construction evidence is conditional or unconditional, and the difference must be
recorded.** `crossprod(A)` is hermitian by construction whatever is or is not known about
`A`, so its evidence has no dependencies. `A + B` is hermitian by construction only if
both `A` and `B` are, so its evidence carries them:

```r
evidence(source = "construction", guarantee = "identity",
         depends_on = list(<A's hermitian evidence>, <B's hermitian evidence>))
```

`evidence_satisfies()` recurses into `depends_on`. Without this, `crossprod(A) + crossprod(B)`
would report bare `construction` even when both inputs were `user_declaration`, which is
exactly the laundering the evidence system exists to prevent.

Examples: `crossprod(A)` hermitian, `source = construction`, no dependencies; a `dsCMatrix`
adapter reports `adapter_contract`; a function-backed operator flagged hermitian by its
author reports `user_declaration`, and after `verify()` gains a second
`source = probe, guarantee = heuristic` record beside it rather than replacing it.

Capability set:

```
hermitian   symmetric   positive_definite   positive_semidefinite
orthogonal/unitary   diagonal   triangular_upper   triangular_lower   full_rank   real
```

`hermitian` and `symmetric` are separate, since `A = A^H` and `A = A^T` differ for complex
operators. Optional capabilities, present or absent: `solve(X)`, `solve_adjoint(X)`,
`diagonal()`, `norm_bound(type)` with evidence, `materialize(target)`.

### 5.4 Transpose, adjoint and conjugate are three different things

Base R distinguishes `t(A)` from `Conj(t(A))` for complex matrices. `linop` preserves that
distinction; an earlier draft made `t()` error on complex operators, which breaks the
matrix contract that is the package's selling point.

```r
t(A)          # transpose,  A^T
Conj(A)       # elementwise conjugate
adjoint(A)    # conjugate transpose, A^H
```

For a real operator, `adjoint(A)` simplifies to `t(A)` and `Conj(A)` to `A`. The concepts
stay separate, since nonsymmetric algorithms and complex bilinear problems need both.

The trap this exposes, worth writing down: for complex `A`, `crossprod(A) = A^H A` is
hermitian, while `t(A) %*% A = A^T A` is symmetric and generally not hermitian. Under a
merged flag these are indistinguishable and a solver takes the wrong path.

### 5.5 dtype lattice

```
promote(double, double) = double
promote(double, complex) = complex
promote(complex, complex) = complex
```

A real operator applied to complex `X` promotes the result. No path silently downcasts.

### 5.6 Flag propagation

| Expression | hermitian | symmetric | positive definite |
|---|---|---|---|
| `A + B` | both | both | both |
| `a * A` | A herm. and `a` real | A sym. | A PD and `a > 0` |
| `A %*% B` | `NA` unless `B` provably `A^H` | `NA` unless `B` provably `A^T` | `NA` |
| `crossprod(A)` (`A^H A`) | always | iff A real | iff full column rank, else `NA` |
| `t(A) %*% A` | iff A real | always | `NA` |
| `adjoint(A)` | preserved | `NA` unless real | preserved |
| `t(A)` | `NA` unless real | preserved | preserved if also hermitian |
| `Conj(A)` | preserved | preserved | preserved |
| `kron(A, B)` | both | both | both |
| `block_diag(...)` | all | all | all |

Every propagated capability records `source = "construction"` with `depends_on` holding the
input capabilities the rule relied on, empty when the rule is unconditional (section 5.3).

### 5.7 Node types

v0.1 surface, deliberately small:

```
leaves:  dense  sparse  fun  identity  diag
unary:   transpose  adjoint  conjugate  scale
n-ary:   sum  product
```

Deferred to v0.2: `lowrank`, `perm`, `zero`, `power`, `inverse`, `hstack`, `vstack`,
`blockdiag`, `kron`.

The reason for deferring is not effort. The extension registry must be exercised by an
**external** adapter before the taxonomy is filled in from inside, otherwise the mechanism
that the whole ecosystem proposition depends on ships untested where it matters. A smaller
complete algebra with a proven registry beats a complete taxonomy with an unproven one.

Registration, unchanged for internal and external node types:

```r
linop_register_node(node = "toeplitz", apply = , adjoint = , transpose = ,
                    caps = , cost = , materialize = )
```

### 5.8 Simplification

Local, applied at construction, never cost-model driven:

```
I %*% A -> A                    A %*% I -> A               0 + A -> A
a * (b * A) -> (a*b) * A        diag(u) %*% diag(v) -> diag(u*v)
transpose(transpose(A)) -> A    adjoint(adjoint(A)) -> A
Conj(Conj(A)) -> A              adjoint(transpose(A)) -> Conj(A)
transpose(adjoint(A)) -> Conj(A)
adjoint(diag(d)) -> diag(Conj(d))     transpose(diag(d)) -> diag(d)
```

For a real operator, `adjoint` and `transpose` have identical action and `Conj` is the
identity, but the node is **not** rewritten at construction. The simplification happens at
apply time. Otherwise `print()` shows `transpose` where the author wrote `adjoint`, the
expression tree stops reflecting what was written, and a real-valued subgraph inside a
complex expression loses the distinction just where it starts to matter.

### 5.9 Materialisation invariant

A literal "never materialise" rule contradicts every Krylov implementation, and so does
the narrower "nothing may be `O(n)`", since every Krylov method holds a basis of `n x m`
and block methods hold `n x b` workspaces. Those are solver state, not operator
materialisation. Three categories:

**User-operator materialisation.** Forbidden without an explicit request. No dense or
sparse representation may be constructed whose row and column extents *both* scale with an
input-space dimension of the operator.

**Tall solver workspace.** Allowed: `n x O(1)`, `n x O(k)`, `n x O(m)`, `n x O(b)`. Krylov
bases, Ritz vectors, residual blocks, scratch.

**Projected problems.** Allowed: `O(m) x O(m)`, `O(k) x O(k)`, `O(b) x O(b)`. Hessenberg,
tridiagonal, bidiagonal, Rayleigh-Ritz, small Gram problems.

Checkable form:

> No internal operation may materialise a representation with both axes proportional to a
> full input-space dimension, except through an explicit materialisation request. Solvers
> report the dimensions and memory cost of every workspace they allocate.

This has one consequence worth stating in the open, because it is an algorithmic decision
and not a formality: for an `m x n` operator, forming the Gram matrix `A^H A` is `n x n`
and therefore prohibited as an automatic step. Normal-equations routes to the SVD may
*apply* `A^H A` as an operator, which is free, and may not silently build it. A route that
needs the explicit Gram matrix has to say so and take the user's request.

Explicit user escapes:

```r
as.matrix(A)                      # dense, size guard the user can raise
as_sparse(A)
collapse(A, expected_applies = 500)
explain(B)
```

`collapse()` returns a transformation report:

```
collapsed:     product of two 80 x 80 dense leaves
retained:      100000 x 80 sparse centring expression
memory added:  51 KB
break-even:    14 applications
```

Cost optimisation becomes an inspectable transformation rather than hidden behaviour, so a
wrong cost model costs a suboptimal plan and never a surprise 40 GB allocation.

### 5.10 Conformance suite

The adoption mechanism. A third party writes `linop.myclass()`, runs one function, and
knows every solver in the ecosystem will work.

```r
verify(A, tol = , n_probe = , seed = )
```

1. **Adjoint consistency.** `|<Ax, y> - <x, A^H y>| <= tol * ||A|| ||x|| ||y||`. Catches
   the most common third-party bug, a wrong or missing conjugation.
2. **Transpose consistency**, separately, for complex operators, plus the three identities
   that tie the nodes together: `adjoint(t(A)) == Conj(A)`, `t(adjoint(A)) == Conj(A)`,
   `adjoint(Conj(A)) == t(A)`. For real operators, that `adjoint` and `transpose` agree in
   action while remaining distinct nodes in the printed tree.
3. **Linearity.** `A(ax + by) == a Ax + b Ay`. Also the test that catches an author who
   wrapped an adaptive solve in a `linop`.
4. **Block consistency.** Column `j` of a block apply equals the single-column apply.
5. **Shape and dtype** for real and complex `X`.
6. **Declared capabilities** probed; contradictions are errors, non-contradictions upgrade
   evidence to `probe_checked` and never further.
7. **Purity.** Two applies to the same `X` give identical output; `X` unmodified.
8. **`apply_gemm` agreement** with the synthesised form.
9. **Materialisation agreement** with column-by-column application.

### 5.11 Provenance is an opaque envelope

Core carries provenance and never inspects it.

```r
provenance = list(provider = "linop.hilbert", payload = <opaque>)
```

The provider registers generics:

```r
provenance_lift(p, x)                 # E_n: C^n -> H
provenance_refine(p, n_new)
provenance_original_residual(p, result)
provenance_summary(p)
```

Hard-coding a schema naming `HilbertOperator`, `Discretisation`, `lift` and `project`
would leak Layer 3 concepts into core and invites retained object graphs, closures
capturing large environments, serialization problems, circular references and unstable
hashed representations. The envelope avoids all of it and keeps the package split honest.

Rules: solvers propagate provenance through results and never depend on it;
`strip_provenance()` exists; `NULL` is the normal value.

This slot goes in during Phase 1. It costs nothing when unused, and fixing the extension
mechanism before adopters exist is much cheaper than after.

---

## 6. Certificates

Two dimensions, not one. A status says whether a check passed; evidence says what kind of
argument produced it. Two `qualified` results can otherwise mean a probabilistic norm
bound with a stated failure probability, empirical stability between `n = 128` and
`n = 256`, or an unsupported heuristic. Those must not flatten together.

```
check                   status        source        guarantee              conf
arithmetic floor        pass          computation   identity               1
residual                pass          computation   identity               1
orthogonality           pass          computation   identity               1
backward error          qualified     computation   probabilistic_bound    0.99
target identity         not_checked   -             -                      -
convergence             pass          computation   identity               1
discretisation error    pass          theorem       deterministic_bound    1
basis tail              qualified     computation   estimate               NA
-------------------------------------------------------------------------------
overall                 qualified     no deterministic bound on: basis tail
```

Status: `pass`, `fail`, `qualified`, `not_checked`. Evidence uses the same three-field
structure as capabilities (section 5.3), for the same reason: the classes are not
comparable on one axis, so a certificate reader needs the source and the guarantee
separately rather than a single strength label. A solve certificate is `qualified` and not
`pass` whenever `forward error` is unchecked, which matrix-free it always is: the forward
bound needs `||A^-1||` and no residual implies one. The example above already has that
shape, through `target identity`; it is the expected outcome rather than a degraded one.

**Every residual and backward-error line carries an arithmetic floor.** S0.6 measured what
happens without one: the a posteriori bound keeps decaying past machine epsilon while the
true error plateaus near 1e-15, so the results that certify as `fail` are exactly the most
converged ones. An earlier version of this table had no roundoff term anywhere.

For a linear solve the floor has a sharper origin than the general `c * ||A|| * eps`.
Higham bounds the computed residual: `fl(b - A x)` differs from `b - A x` by at most a
small multiple of `eps * (||A|| ||x|| + ||b||)`, so no `x`, however exactly it solves the
system, can be shown to have a residual below that level. The Rigal-Gaches identity gives
the normwise relative backward error as `||b - A x|| / (||A|| ||x|| + ||b||)` with the same
denominator, so once the backward error is expressed relatively the floor is a plain
`c * eps` and one quantity serves both lines.

`c` counts terms in an inner product a matrix-free operator never exposes, so it is an
argument rather than a constant, defaulting to 4 from S0.6's measurement. A line that meets
its tolerance only through the floor is `qualified` with an `estimate` guarantee, never a
clean `identity`: the floor is a model of roundoff and this table already has the field to
say so. The `arithmetic floor` row reports `||A||`, the route that produced it and the
resulting floor, and inherits that route's evidence, so an estimated norm is visible in the
summary line rather than buried.

`||A||` comes from `norm2()`, which is exact where the node structure gives it exactly,
exact for a stored matrix small enough to factor, and a power iteration otherwise. Every
route returns a lower bound, which is the safe direction: it shrinks the Rigal-Gaches
denominator, so the reported backward error can only overstate. Section 5.3 lists
`norm_bound(type)` with evidence as an optional operator capability; that protocol can
adopt this estimate later, but the floor needs `||A||` in Phase 2 and cannot wait for it.

Overall status is derived by a documented rule, not case by case: `fail` if any check
fails; else `qualified` if any check is `qualified` or if any check the request implies is
`not_checked`; else `pass`. The summary line names which checks lack a deterministic
bound, since an all-`pass` certificate resting on estimates is a different object from one
resting on theorems, and "weakest evidence" is not well defined once evidence stopped
being a ranking.

### 6.1 Target identity defaults to not_checked

A small residual proves that `(lambda, v)` is approximately an eigenpair. It does not
prove that it is the largest, smallest or nearest requested eigenpair. Solver convergence
history is not evidence for target identity, and treating it as such is the most common
way a diagnostic object gets mistaken for a certificate.

`target_identity` passes only with:

- inertia or Sturm counts, from an `LDL^T` of `A - sigma I`, which genuinely certifies how
  many eigenvalues lie below `sigma` when a sparse factorization is affordable
- interval enclosures
- complete decomposition at small dimension
- a rigorous separation or gap bound
- a theorem tied to the specific solver and operator class

For matrix-free `which = "LM"` there is generally no cheap certificate, so `not_checked`
is the honest and expected answer. This restraint is what separates the certificate from a
generic diagnostics object, and it should be documented as a feature rather than
apologised for.

`verify(fit, A)` is callable on results produced by other packages, given an operator and
claimed eigenpairs. That independence is the point, and it is the second method on the
same generic that checks operators (section 1.1): internally different checks, publicly
one question.

---

## 7. Solvers

### 7.1 Pure-R Krylov, v0.1

| Method | Problem | Notes |
|---|---|---|
| `cg` | SPD square | requires `positive_definite` at declared minimum evidence |
| `minres` | hermitian indefinite | not on CRAN elsewhere; careful recurrence and breakdown tests |
| `gmres(m)` | general square | restarted, MGS with reorthogonalisation |
| `fgmres(m)` | general square, flexible preconditioner | right preconditioning only; stores the preconditioned basis |
| `bicgstab` | general square | |
| `lsqr` | rectangular least squares | not on CRAN elsewhere |
| `lsmr` | rectangular least squares, ill-conditioned | not on CRAN elsewhere; stopping tests are the risk |

Seven methods. An earlier version of this table listed six and omitted FGMRES, which was
the outlier rather than a decision: section 4.3 always carried an FGMRES row, section 16
cites Saad on FGMRES as the source of the `fixed = FALSE` rules, and the enforcement code
and its Gate 1 test both name it. FGMRES is also the only v0.1 consumer of the
`variable_inexact` contract of section 4.2, so without it the inexact half of the fidelity
lattice ships with nothing that reads it until Phase 3. It shares the GMRES implementation
rather than sitting beside it: the Arnoldi loop, the rotations and the residual recurrence
are the same, and the difference is storing the preconditioned basis `Z` and forming the
update from it, at the cost of one further `n x m` array.

```r
solve(A, b, method = "auto", preconditioner = NULL, tol = , maxit = , x0 = , control = )
```

`method = "auto"` picks from declared capabilities and shape, records the choice and the
reason, and never picks a method whose precondition is `NA` or whose evidence is below the
method's declared minimum. It errors and names the capability it would have needed.

**The evidence minimum filters `auto` and does not gate the named method.** These are two
different questions and only the second is a request the caller made. `solve(A, b,
method = "cg")` requires `positive_definite` to be `TRUE` behind any evidence, because
naming the method is the caller asserting their own declaration and refusing it would make
`linop(M, properties = c(positive_definite = TRUE))` unsolvable by the method it was
declared for. `method = "auto"` is the package choosing on the caller's behalf and applies
the minimum. CG's admits `construction`, `adapter_contract`, `theorem` and `computation` at
an `identity` guarantee: an exact check on data the operator already holds is an identity,
not a probe, and `linop_diag()` proves definiteness from the signs of `d`. Probes are
excluded by the guarantee field, which they never satisfy, rather than by the source list.

MINRES and LSMR are the two with real numerical risk, not because of size but because
correct recurrence, breakdown detection, condition estimation and stopping tests take
care. Gate 2 covers this with reference-implementation agreement on ill-conditioned
problems.

### 7.2 Eigen and SVD

v0.1: reference implementations (Lanczos with full reorthogonalisation,
Golub-Kahan-Lanczos), labelled reference rather than production, plus RSpectra delegation
where the operator can present the interface RSpectra wants.

v0.2: `linop.primme`. Hermitian standard and generalized, largest, smallest and interior
singular triplets, real and complex, block methods, preconditioning, warm starts.

---

## 8. The Hilbert layer

### 8.1 The rule

> Infinite-dimensional operators are never solved by pretending they are finite. They are
> discretised into finite linear operators, and the discretisation becomes part of the
> result and part of the certificate.

```
HilbertOperator
      | discretise(scheme, n)
linop  (with provenance)
      | solve
FiniteResult
      | lift + certify against the original operator
HilbertResult
```

### 8.2 The hierarchy

The finite operator is not a sibling of the Hilbert operator. It is its **image under a
discretisation**:

```
A_{m,n} = P_m A E_n,     E_n: C^n -> H,     P_m: K -> C^m
```

Sibling classes cannot express that. In practice `discretise()` returns a `linop` whose
provenance envelope holds whatever `linop.hilbert` needs to lift, refine and compute an
original-operator residual. Without that link, certifying against the original operator is
not implementable and the certificate collapses to certifying `A_n`, which is the thing
the layer exists to improve on.

### 8.3 Spaces are first-class

```
Space
├── FiniteSpace
└── HilbertSpace
    ├── SequenceSpace     (ell^2, weighted ell^2)
    ├── FunctionSpace     (Fourier, Chebyshev, Ultraspherical, Hermite, Laguerre)
    ├── SobolevSpace      (H^1, H^1_0, H^k)
    └── ProductSpace
```

A `HilbertSpace` carries scalar field, domain, inner product, basis family, norm,
coefficient representation, and declared conversions into other spaces.

The domain/range distinction is load-bearing and ApproxFun has proven the design.
`Derivative(Chebyshev())` maps into `Ultraspherical(1)`, not back into Chebyshev.
Pretending every operator maps one coefficient space to itself destroys the banded
structure that makes spectral methods fast. So `linop.hilbert` needs a conversion graph
between spaces, with composition, and arithmetic that inserts conversions.

### 8.4 Bounded and closed operators

```
HilbertOperator
├── BoundedOperator     A: H -> K, with a norm bound
└── ClosedOperator      A: D(A) subset H -> K, closed, densely defined
```

`D(A)` is an essential part of a `ClosedOperator`, governing composition, addition, the
domain of `A + B`, the domain of `A^*`, invertibility, and whether a spectral request is
meaningful.

`adjoint(A)` must produce the adjoint action, the adjoint domain and the boundary form.
For a differential operator, changing boundary conditions changes the adjoint and can
change the spectrum completely.

**Adjoint domain is usually not computable, so unknown is first-class.** Symbolic
integration by parts is tractable for constant or polynomial coefficients with separated
boundary conditions on an interval, and not much further. Outside that class
`adjoint_domain = NULL`, and every solver needing self-adjointness refuses rather than
assuming. Same discipline as section 5.3, applied twice rather than invented twice.

### 8.5 Forms are peers, not subclasses

A sesquilinear form `a(u, v)` maps `H x H -> C`. It is not `H -> K`. It induces an
operator into the dual, with the Riesz map connecting them. Forcing a form into the
operator callback contract means the mass matrix in

```
a(u, v) = lambda m(u, v)  for all v      ->      A_n x = lambda M_n x
```

gets attached by hand instead of falling out of the structure.

```
HilbertOperator     (peer)
SesquilinearForm    (peer)
    induced_operator(a, riesz_map) -> HilbertOperator
    discretise(a, m, scheme)       -> list(A_n, M_n), shared provenance
```

### 8.6 Discretisation is an object

```r
discretisation(trial_space, test_space, basis, project, reconstruct,
               quadrature, boundary_conditions, n, refinement_rule,
               approximation_order)
```

One object covers Galerkin, Petrov-Galerkin, spectral, finite element, collocation and
finite section, since they differ only in the trial and test maps. v1 schemes:
`FiniteSection`, `Galerkin`, `Collocation`.

### 8.7 Compose before discretising

For `L = -D^2 + V`, build `L` continuously, then discretise once. Continuous composition
preserves self-adjointness, differential order, boundary conditions, sparsity in the right
basis, compactness, coercivity and known spectral properties. Discretising first destroys
most of that, silently. Enforced structurally: arithmetic on `HilbertOperator`s stays in
the Hilbert layer, and no operation returns a mixed object.

### 8.8 The three-part certificate

```
total error  <~  algebraic solver error
               + discretisation error
               + representation and quadrature error
```

A tiny residual for `A_n` proves the finite problem was solved accurately and says nothing
about whether `A_n` represents `A`. The `original operator residual` line is the only one
that speaks about the mathematical problem the user posed, and it is computable only
through provenance.

### 8.9 Spectral requests refuse by default

The verified facts make this a refusal rather than a warning. Pollution can occur at any
point in a gap of the extended essential spectrum, and for a self-adjoint operator not
bounded above, arbitrarily small basis perturbations push the finite-section spectrum
arbitrarily far. A polluted eigenvalue is fabricated and passes every check computed on
`A_n`.

`eigs()` on a `ClosedOperator` errors unless the operator or scheme declares justifying
structure: compact; self-adjoint with compact resolvent; coercive generalized problem;
target interval isolated below the essential spectrum; or a scheme declared pollution-free
for the class. Otherwise:

```r
spectral_approximation(A, scheme, n)
```

which returns an object that is not called an eigenvalue and whose print method says so.
Two verbs, because a user reading only the return value must not be able to confuse them.

### 8.10 First rigorous unit

Not banded Toeplitz generally. Two mathematical problems with that choice:

- For a **self-adjoint** banded Laurent operator the spectrum is the range of the symbol
  and is purely **continuous**, with no isolated eigenvalues at all. A flagship for an
  eigensolver that has no eigenvalues is the wrong flagship.
- For **nonnormal** Toeplitz the spectrum is not the symbol range (it is the range
  together with the points of nonzero winding number), and finite sections behave badly.
  The "known answer" is not known in the simple sense the earlier draft assumed.

The first rigorous class instead:

> **Self-adjoint banded Jacobi or Laurent operators on `ell^2`, plus finite-support real
> perturbations producing isolated eigenvalues outside the essential spectrum.**

This gives a known essential spectrum from the limiting coefficients via Weyl, genuine
discrete eigenvalues to approximate, a meaningful isolation gap, computable residuals, a
tractable route to truncation and enclosure bounds, and a natural pollution test inside
the gap.

Canonical instance: the discrete Schrödinger operator on `ell^2(Z)`, in the convention

```
(H u)_j = u_{j-1} + u_{j+1} + V_j u_j
```

so the free operator `V = 0` has purely continuous spectrum `[-2, 2]`. The convention is
fixed here because the tests assert specific numbers.

**Finite support, not general compact `V`, in the first release.** Compact is
mathematically clean and still too broad to yield a computable truncation bound from the
actual object. With `V` of finite support, outside the support the eigenvalue equation is
exactly the free one, so an eigenfunction at `lambda` outside `[-2, 2]` decays exactly
exponentially at a rate obtainable in closed form from `lambda`. The finite-section
truncation error then has a closed form rather than an existence proof. That is a
theorem-sized first release.

The next extension is potentials satisfying a declared decay envelope
`|V_j| <= C exp(-alpha |j|)`, with `C` and `alpha` carried on the discretisation object so
the tail bound stays operational rather than existential. General compact `V` comes after
that, if at all.

The flagship demonstration contains both halves of the philosophy:

1. **A request that rigorously succeeds.** `H` with finite-support real `V`, an isolated
   bound state outside `[-2, 2]`, certified with a closed-form truncation bound and an
   isolation gap.
2. **A request that cleanly refuses.** The free operator `V = 0`, purely continuous
   spectrum on `[-2, 2]` and no isolated eigenvectors: `eigs()` errors,
   `spectral_approximation()` returns something not called an eigenvalue.
3. **A pollution test.** A self-adjoint operator with an essential-spectrum gap where
   naive finite section is known to produce plausible spectral points inside the gap,
   asserting that the refusal fires and that the approximation verb labels its output
   correctly.

Together these show the layer performing mathematical classification rather than feeding
larger matrices to an eigensolver.

### 8.11 Function spaces after that

Fourier on the circle, then Chebyshev and Ultraspherical with the conversion graph, then
Hermite and Laguerre. One-dimensional domains only. This is the larger lift and should not
start until the sequence layer's certificate is rigorous and tested.

---

## 9. Phases and gates

### Phase 0: spikes

Written up in `dev_notes/` before Phase 1. None produce shipped code.

**S0.1 `%*%` dispatch and the R version floor.** R made `%*%`, `crossprod()` and
`tcrossprod()` generic at some version (believed 4.3.0, `matrixOps` group generic for S4).
Verify from R NEWS, build a 20-line test package under S3 and under S4, check dispatch on
the oldest R that must be supported. Also check `t()`, `Conj()` and `solve()` dispatch on
the same footing. Decides S3 vs S4 and the `Depends` floor.

**S0.2 R callback overhead.** Nanoseconds per apply for a `fun` leaf at `n` = 1e3, 1e5,
1e7, block widths 1, 4, 16, 64, against native `dgemm`. Krylov runs 1e2 to 1e5 iterations,
so if the round trip is expensive relative to the matvec, block widths become mandatory
and some solver designs change.

**S0.3 Vendored PRIMME build.** Build the static library under Rtools on Windows, macOS
arm64 and Linux, linking R's BLAS/LAPACK. No bindings.

**S0.4 eigencore source audit.** Full read, not grep. Findings to
`dev_notes/eigencore-audit.md`, nothing quotable without a line number.

**S0.5 Prior art sweep.** Confirm MINRES, LSQR and LSMR are genuinely absent from CRAN.
Read the Numerical Mathematics task view. Confirm no R package already exposes a
function-space or operator-discretisation abstraction worth being compatible with.

**S0.6 Finite-section bound feasibility.** For the discrete Schrödinger operator with real
finite-support `V`, derive the closed-form exponential decay rate outside the support at a
given `lambda`, the resulting finite-section truncation bound, and the isolation gap
bound. Check each is computable from quantities the discretisation object will hold. If
any is not, section 8.10 needs a different first class, and better to learn that now than
in Phase 5.

### Phase 1: core

Object model, expression graph, the v0.1 node set, algebra, capability and evidence
lattices, dtype lattice, adapters for `matrix`, `Matrix` classes and `function`, tree
printing, `verify()`, the provenance envelope, and the `linop` / `linsolve` /
`preconditioner` split internally, even where only operators are implemented. The public
surface stays at the six names of section 1.1 throughout.

**Gate 1 to 2.**
- Conformance suite passes on every node type, real and complex.
- Expression fuzzer: 1e4 generated trees checked against a base R recomputation from
  materialised leaves, agreement to 1e-12 relative.
- Flag propagation table checked by brute force, including the `A^H A` versus `A^T A`
  distinction on complex inputs.
- Evidence propagation checked: every conditional construction rule records `depends_on`,
  and `evidence_satisfies()` recursing through a composite of `user_declaration` inputs
  fails a requirement that the inputs would have failed directly.
- `linsolve` conformance is a separate suite; running the `linop` suite against a
  `linsolve` errors rather than failing purity.
- **API budget enforced.** A test asserts the exported name set against the budget in
  section 1.1. Adding an export is a deliberate act that fails a test until the budget is
  edited, which is the only mechanism that reliably stops a public surface from growing by
  accretion.
- dtype promotion table exhaustive.
- Install cost recorded via `tools::package_dependencies(recursive = TRUE)`.

### Phase 2: solvers, ships as v0.1

Seven Krylov methods, the preconditioner model, certificates with evidence, the `solve()`
method, vignettes, pkgdown. `linsolve()` and `solver()` exist internally and are not
exported.

**Gate 2 to 3.**
- Every solver has a recovery test against closed-form ground truth, run to convergence.
- Certificate coverage: over at least 20 seeds per class, the claimed bound contains the
  true error at least at the nominal rate.
- MINRES and LSMR agree with a trusted reference on ill-conditioned problems, including
  breakdown and near-breakdown cases.
- Every solver refuses a preconditioner whose declared properties it requires and does not
  have, and refuses a side it is not defined for, tested per row of the table in 4.3.
- FGMRES with a fixed right preconditioner reproduces right-preconditioned GMRES to
  machine precision. Two paths that should be identical, so a divergence there is a bug
  and not a tuning question.
- Benchmark harness runs end to end with committed results.

### Phase 3: `linop.primme`, ships as v0.2

Vendored source, bindings using the block matvec callback, backend registration with the
full protocol, warm starts, preconditioning, complex support. The deferred node types
(`kron`, stacks, `lowrank`, `blockdiag`, `perm`, `power`) land here, after at least one
external adapter has exercised the registry.

**Gate 3 to 5.**
- Source install succeeds on Windows, macOS arm64 and Linux unaided.
- The stdout NOTE is fixed at source, not suppressed.
- Every PRIMME path covered by a recovery test against a closed-form spectrum.
- `plan_eigs()` reports correctly with two backends installed.

### Phase 4: adapters

`DelayedArray`, `HDF5Array`, `torch`, memory-mapped. Each a `Suggests` with conditional
`registerS3method()` in `.onLoad()`.

### Phase 5: `linop.hilbert` first unit

Spaces, `BoundedOperator`, `ClosedOperator`, `SesquilinearForm`, `FiniteSection` and
`Galerkin`, lifting, the three-part certificate, and the self-adjoint Jacobi plus compact
perturbation class with proven bounds.

**Gate 5 to 6.**
- Original-operator residual computed from provenance for every supported class.
- Truncation and isolation bounds proven, cited, and shown to hold empirically across
  seeds and truncation sizes.
- The three-part flagship of 8.10 passes: rigorous success, clean refusal, pollution test.
- Adaptive refinement to a requested tolerance converges on the supported class.

### Phase 6: function spaces and pollution-free spectra

Fourier, then Chebyshev and Ultraspherical with the conversion graph, then second-order
relative spectra and dissipative-perturbation methods. Nonsymmetric Krylov-Schur and
harmonic extraction on the finite side land here.

---

## 10. Testing

Shape tests prove plumbing and nothing else. Every solver needs recovery against known
truth.

### Closed-form ground truth, finite layer

- **1-D Dirichlet Laplacian** `tridiag(-1, 2, -1)`: eigenvalues `4 sin^2(k pi/(2(n+1)))`,
  eigenvectors `sin(i j pi/(n+1))`. Exact for any `n`, `k`.
- **Circulant matrices**: eigenvalues are the DFT of the first column. Exact, and complex.
- **Prescribed spectrum** `Q diag(lambda) Q^H`: clustering, exact multiplicity and tiny
  gaps dialled in deliberately.
- **Complex non-hermitian with `A^T A` symmetric but not hermitian**, to test the
  transpose/adjoint split end to end.
- **Hilbert matrices** for LSQR and LSMR stress.
- **KMS matrices** `rho^|i-j|`, known tridiagonal inverse, for `solve` paths.

### Closed-form ground truth, Hilbert layer

- **Discrete Schrödinger** `(Hu)_j = u_{j-1} + u_{j+1} + V_j u_j` on `ell^2(Z)`, `V` real
  with finite support: essential spectrum `[-2, 2]`, isolated bound states outside it,
  eigenfunctions decaying exactly exponentially outside the support at a rate available in
  closed form.
- **The same operator with `V = 0`**: purely continuous spectrum on `[-2, 2]`, no
  eigenvalues. The refusal test.
- **Multiplication operator on Fourier space**: spectrum is the essential range of the
  multiplier. Second refusal test.
- **Harmonic oscillator on Hermite space**: eigenvalues `2k + 1`, self-adjoint with
  compact resolvent. The canonical justified request.
- **An operator with an essential-spectrum gap** where naive finite section pollutes.

### Categories reported separately

Structure and dispatch; algebra correctness; adjoint and transpose consistency; linearity
(which is also the test that catches an adaptive solve wrapped as a `linop`); parameter
recovery; certificate coverage; evidence-propagation correctness; reference-implementation
agreement; convergence-rate tests (CG on SPD with known condition number should converge
near the predicted rate, catching implementations that are wrong but pass a tolerance
test); discretisation convergence order against the scheme's declared
`approximation_order`.

README and paper report the breakdown by category, never a test count.

---

## 11. Benchmark protocol

Primary metric: **time to independently certified answer**.

Reported per run: wall time, matvec count, preconditioner applies, peak memory, achieved
residual, orthogonality loss, target recovery, certificate status and which checks lacked
a deterministic bound. Matvec count is the machine-independent number and gets equal
billing with wall time.

Comparators: RSpectra, irlba, eigencore, vendored PRIMME, Rlinsolve on the linear-solve
side, `Matrix::solve` where applicable. On the Hilbert side there is no R comparator, so
compare against ApproxFun and Chebfun on identical problems, stating that the comparison
crosses languages.

Rules: seeded, scripted, committed, machine spec recorded. No caps on problem size in any
published run. Cache filenames keyed on every knob setting problem extent. Losses
published alongside wins.

---

## 12. The PRIMME decision

**A. Take over the CRAN package.** Policy: "When a new maintainer wishes to take over a
package, this should be accompanied by the written agreement of the previous maintainer
(unless the package has been formally orphaned)." The maintainer's email is undeliverable,
so written agreement is not obtainable by the ordinary route. Whether CRAN has formally
orphaned it needs checking against the orphaned list.

**B. Vendor the C source.** BSD 3-clause permits it. The archived package's clean Windows
checks say it builds under Rtools. Only route giving a strict dependency, since orphaned
packages must not be strict requirements including indirectly.

**C. Both.** Vendor for `linop.primme`, separately offer to adopt the CRAN package as a
service to its reverse dependencies.

Recommendation: **B**, with C as a courtesy afterwards. Commit to no PRIMME interface
language until S0.3 produces a linking build on all three platforms.

---

## 13. Risks

1. **The engine is orphaned upstream at the R level.** Mitigated by vendoring, which moves
   maintenance onto this project. The C library is active, so this is a wrapper problem.
2. **R callback overhead.** S0.2 decides before any solver is written. Contingency is a
   small C fast path for `fun` leaves, which would move core off zero-compiled-code and
   must be a deliberate decision rather than a drift.
3. **Windows source install of a vendored C library** is the main adoption risk for the
   backend and none for core. That asymmetry is the argument for the split.
4. **MINRES and LSMR correctness.** Not size, but recurrence, breakdown detection,
   condition estimation and stopping tests. Gate 2 requires reference agreement including
   near-breakdown cases.
5. **The Hilbert layer is a research programme.** ApproxFun is roughly fifteen years on top
   of Chebfun's prior decade. Control is structural: separate package, separate cadence,
   first unit is one class with provable bounds, nothing in core or v0.1 waits on it.
6. **Spectral pollution is a correctness hazard.** A polluted eigenvalue passes every check
   computed on `A_n`. Mitigated by refusing rather than warning, and the two-verb split.
7. **Competitive pressure on the finite side.** Only a risk if linop is built on the same
   axis as eigencore. Acceptance criteria below are deliberately not benchmark-relative.

---

## 14. Acceptance criteria

**Primary, v0.1:**

> A third party implements `linop.myclass()` for a storage format this project has never
> seen, passes `verify()`, and immediately gets all six Krylov solvers, preconditioners
> and certificates working, with no changes to `linop`.

**Secondary, v0.1:** every solver has recovery and coverage tests; install cost is
`methods` only; MINRES, LSQR and LSMR available operator-natively; every certificate line
carries source, guarantee and confidence; `target_identity` reports `not_checked` wherever
it is not genuinely certified; no `linsolve` reaches an eigensolver whose declared
accepted contracts exclude it.

**v0.2:** `linop.primme` installs from source on all three platforms unaided; complex
matrix-free SVD and preconditioned generalized eigenproblems work end to end;
`summary(fit)` names the backend and the reason it was chosen; published benchmark table
including losses.

**v0.3, the Hilbert unit:**

> For the discrete Schrödinger operator `(Hu)_j = u_{j-1} + u_{j+1} + V_j u_j` on
> `ell^2(Z)` with real finite-support `V`, `linop.hilbert` returns an isolated bound state
> outside `[-2, 2]` whose certificate reports the original-operator residual, a
> closed-form finite-section truncation bound and an isolation gap, with total error
> decomposing into algebraic, discretisation and representation parts that sum within the
> claimed tolerance. On the same operator with `V = 0`, the request refuses.

Falsifiable, and it contains both halves of the philosophy.

**The real win condition**, not a benchmark: eigencore, RSpectra or irlba implementing the
`linop` contract.

---

## 15. Presentation constraints

- No claim about another package's internals without a line-number citation from a full
  read.
- No competitive vocabulary in any public artefact. Describe what linop does, report both
  numbers side by side, let readers compare. Buchsbaum, You, and the RSpectra, irlba,
  ApproxFun and Chebfun authors are peers whose adoption or goodwill is the goal.
- No negative selling. The README opens with what linop does.
- Reach out to Buchsbaum after v0.1 exists and before any benchmark is published, with the
  interface proposal rather than a comparison.
- The Hilbert layer must never be described as certifying more than it certifies. The
  refusal behaviour and the `not_checked` default on target identity are features to
  document prominently, not apologise for.

---

## 16. Theory references

To be verified individually before any appear in a paper. The source documents' citations
were partly mis-attributed: Osborn 1975 is *Math. Comput.* 29:712-725, not the SIAM URL
given, and the Finite Element Exterior Calculus citation for unbounded-operator domain
theory is off-target, since FEEC concerns differential complexes and structure-preserving
discretisation.

- Osborn, "Spectral approximation for compact operators", *Math. Comput.* 29 (1975),
  712-725.
- Chatelin, *Spectral Approximation of Linear Operators*, SIAM Classics.
- Boffi, "Finite element approximation of eigenvalue problems", *Acta Numerica*.
- Davies & Plum, "Spectral pollution", and the second-order relative spectra literature.
- Aurentz & Trefethen, "Block operators and spectral discretizations", *SIAM Review*
  (2017).
- Kirby, "From functional analysis to iterative methods", *SIAM Review*.
- Bottcher & Silbermann on finite sections of Toeplitz operators.
- Teschl, *Jacobi Operators and Completely Integrable Nonlinear Lattices*, for the
  essential spectrum and eigenvalue structure of the first rigorous class.
- Olver & Townsend on the ultraspherical spectral method, for the banded-conversion design.
- Saad on FGMRES and Notay on flexible CG, for the `fixed = FALSE` preconditioner rules.
- Simoncini & Szyld, and van den Eshof & Sleijpen, on inexact Krylov methods and relaxed
  inner tolerances, for the `variable_inexact` transform contract and the role of
  `error_bound = available`.

---

## 17. Open questions

1. Route A, B or C on PRIMME. Recommendation B; confirm before Phase 3.
2. Package names. `linop.primme` and `linop.hilbert` group clearly in listings.
3. Does core stay at zero compiled code if S0.2 says callbacks are slow?
4. ~~Should `linsolve` and `preconditioner` live in core?~~ **Settled: core, through at
   least v0.2.** They are three sides of one contract, not peripheral integrations:
   solvers consume preconditioners, shift-invert consumes solve objects, generalized
   eigenproblems need metric solves, and backend interoperability depends on a shared
   inverse-action contract. Splitting them would make a user install two packages to
   express one iterative solve and force third-party adapters to choose which package owns
   their solve capability. Revisit only if one of these appears: an independent release
   cadence, substantial compiled dependencies, an ecosystem unrelated to `linop`, or
   enough public constructors that they obscure the operator API. Large preconditioner
   families and factorisation backends go elsewhere when they arrive; the protocols stay.
5. JOSS or MEE, and not before demonstrated research use. The Hilbert layer may deserve a
   separate methods paper.
6. Is there value in contacting Sheehan Olver about the space-conversion design before
   building it, rather than reimplementing from documentation?
7. `solve()` on a rectangular operator: recommendation is to error and name the
   least-squares methods, reachable as `solve(A, b, method = "lsmr")`. Base R's `solve()`
   requires square and least squares is a different request, so silently switching meaning
   is worse than one explicit argument. Confirm before Phase 2 freezes the signature.
8. Does `solve()` returning a bare numeric with the certificate in an attribute survive
   contact with real use, given that arithmetic drops it? Revisit after the first
   vignettes are written, when it becomes visible whether users reach for `summary()`
   before or after manipulating the result.
