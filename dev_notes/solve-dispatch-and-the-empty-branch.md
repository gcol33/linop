# solve(), and the branch with a correct rule and no supply

Written 2026-07-28. Section 1.1's verb, over the seven methods of section 7.1. Three findings,
two of which correct the plan.

## 1. It cost no export, and section 1.1's "keeps the verb budget" is literal

`solve` is a base generic, so `solve.linop` reaches it exactly the way `t()`, `crossprod()`
and `%*%` already do. `S3method(solve, linop)` puts nothing into
`getNamespaceExports("linop")`, so `BUDGET` in `test-api-budget.R` did not move to
accommodate the whole of Phase 2's verb surface. The deliberate edit that file exists to
force was to *remove* `"solve"` from `PLANNED_PHASE2` and add the assertion that it works
without being exported, beside the other base generics.

That is the strongest case yet for the section 1.1 rule about reusing R's own generics rather
than coining verbs: the rule is usually argued as familiarity, and here it also turned out to
be free.

`eigs()` and `svds()` will not be free. They are not base generics, they are new names, and
`BUDGET` will have to grow by two when section 7.2 is built.

## 2. The cg branch of auto is correct and has no supply

`method = "auto"` picks from declared capabilities and shape, and never picks a method whose
precondition is `NA` or whose evidence is below that method's declared minimum. CG's minimum
(`CG_PD_REQUIREMENT`) admits `construction`, `adapter_contract`, `theorem` and `computation`
at an `identity` guarantee, and excludes `user_declaration`. That is right, and the plan
argues it well: a caller naming `method = "cg"` is asserting their own declaration, and auto
choosing CG on the strength of that same declaration would be the package treating a claim as
a proof.

Measured across every constructor in the package:

| operator | hermitian | source | positive_definite | source | auto picks |
|---|---|---|---|---|---|
| `linop_eye(n)` | TRUE | construction | TRUE | construction | direct |
| `linop_scaling(d)`, positive `d` | TRUE | computation | TRUE | computation | direct |
| `linop_scaling(d)`, mixed `d` | TRUE | computation | FALSE | computation | direct |
| `linop(M)`, dense SPD | TRUE | computation | NA | - | minres |
| `linop(M)`, dense SPD, declared | TRUE | user_declaration | TRUE | user_declaration | gmres |
| `linop(M)`, dense HPD | TRUE | computation | NA | - | minres |
| `laplacian_1d(n)` | TRUE | user_declaration | TRUE | user_declaration | gmres |
| `shifted_laplacian_1d(n, s)` | TRUE | user_declaration | NA | - | gmres |
| `convdiff_1d(n, mu)` | NA | - | NA | - | gmres |
| `adjoint(X) %*% X` | TRUE | construction | NA | - | minres |
| `as_spd_linop(kms)` | TRUE | user_declaration | TRUE | user_declaration | gmres |
| sum of two declared SPD | TRUE | construction | TRUE | construction | gmres |

**Nothing reaches cg.** The two constructors that establish definiteness at an acceptable
source are `linop_scaling()` and `linop_eye()`, and both are diagonal, so the direct route
takes them first. Every other route to `positive_definite` in the package is a bare user
declaration.

The rule is not the bug. The supply is missing, and there are two ways to close it, both
design decisions rather than fixes:

- a `computation` route for definiteness on a dense hermitian leaf, which means attempting a
  Cholesky factorization at construction. That is O(n^3) where the existing entry checks are
  O(n^2), so it would have to be opt-in or size-gated, and the gate is the design question.
- the product node stamping `positive_definite` on `adjoint(X) %*% X` when `X` is established
  to have full column rank. It already stamps `hermitian` there by construction and
  unconditionally; definiteness is conditional on the rank, which is the same shape as the
  existing rule with one more premise.

Until one of those lands, a caller with a genuinely SPD matrix-free operator gets MINRES from
`auto` if symmetry is computable and GMRES otherwise, and reaches CG by naming it. Both are
correct answers, and the second is slower than it needs to be.

`test-solve-dispatch.R` records this rather than leaving it to be rediscovered: one test
asserts that no constructor reaches the cg branch and says that when a constructor starts
proving definiteness this test is the one that fails and should be replaced rather than
adjusted. The branch itself is covered by an operator whose evidence is stamped directly, so
it is tested rather than assumed.

The last row of the table is the section 5.3 laundering case arriving at dispatch, and it is
worth reading beside the row above it. The sum's `positive_definite` carries a *construction*
source, the same source the branch accepts, and is still refused, because
`evidence_satisfies()` recurses into `depends_on` and finds two user declarations underneath.
Same source, opposite verdict. A test asserts both halves against each other.

## 3. The certificate attribute survives arithmetic, which the plan has backwards

Section 1.1: "`solve()` returns a plain numeric vector or matrix with the certificate in an
attribute. The tradeoff is real and documented: arithmetic on the result drops the
attribute."

Measured, R's arithmetic operators copy attributes from their longer operand, so:

| expression | certificate |
|---|---|
| `x + 0` | kept |
| `2 * x` | kept |
| `x + x` | kept |
| `-x` | kept |
| `t(x)` | kept |
| `X + 0` (matrix) | kept |
| `x[1:3]` | dropped |
| `X[, 1]` | dropped |
| `as.numeric(x)` | dropped |
| `c(x)` | dropped |
| `sum(x)` | dropped |

So the tradeoff is real and it is the other one. Losing the certificate under indexing is the
minor half. The major half is that `2 * x` carries a certificate reporting a residual for a
vector that is no longer the solution: **the certificate can outlive the value it describes.**

That is sharper than what the plan documents, and it does not change the decision, because
the alternative section 1.1 rejects is worse in a way this package cares about more: a
classed result would survive arithmetic *and* break the promise that a solve returns what a
matrix solve returns. But it should be documented as what it is. A test asserts both columns
of the table so the direction stays recorded.

## 4. The direct half of section 1.1, and what it is worth

"If the operator declares a direct solve capability, `solve()` uses it." There is no
capability by that name in `CAPABILITY_NAMES`, and the one that makes a direct solve
available is `diagonal` together with `full_rank`. `diag_solve()` recovers the entries with
one apply of a block of ones rather than reading a leaf's `d`, so the route works for any
operator that can establish diagonality and not only for `linop_scaling()`. It reports zero
iterations and a certificate with the same rows as any other.

It refuses a preconditioner and a starting iterate by name rather than ignoring them: on an
exact solve there is nothing to accelerate and nothing to start from, and silently accepting
either would be the kind of quiet no-op the package is against. Naming an iterative method is
how a caller who wants the preconditioner exercised gets it.

The zero-entry check is on the route rather than on the way in, so it holds whenever a
constructor starts establishing diagonality some other way. It is currently unreachable
through `auto` for the same reason as section 2: a bare declaration of `diagonal` does not
clear the minimum.

## 5. What auto does not choose

**Rectangular is a refusal, not a fallback.** Section 1.1 is explicit and the reason is that
least squares is a different mathematical request rather than the same one on a different
shape. `solve()` on a rectangular operator errors and names `lsqr` and `lsmr`, which are
reachable through the same verb. The roster being total over rectangular shapes and `auto`
being total over them are different statements, and only the first is true.

**GMRES rather than BiCGSTAB as the fallback.** Both require nothing of the operator and both
run without an adjoint, so both would answer. GMRES minimises the residual over the whole
Krylov space, so its residual cannot rise and it has no breakdown; those are the properties
to prefer when the package is choosing rather than the caller. BiCGSTAB is cheaper per step
and in storage and stays reachable by name. The reason is recorded in the certificate, so a
caller who wants the other one can see which was taken and why.

## 6. control

Method-specific knobs go through `control`, and a name the chosen method does not take is an
error rather than being ignored, so a typo is not silent. `control` may not carry anything
`solve()` already has an argument for (`tol`, `maxit`, `x0`, `preconditioner`), so each of
those has one place it is set.

## Files

- `R/solve.R`, `R/certificate.R` (the dispatch line in the print method)
- `tests/testthat/test-solve-dispatch.R`, 22 tests, 106 assertions
- `tests/testthat/test-api-budget.R`, the deliberate edit
- `dev_notes/spikes/auto_reachability.R`
