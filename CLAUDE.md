# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository state

Phase 0 (spikes) and Phase 1 (core) are complete and Gate 1 is met. Phase 2 is under way:
the solve certificate with its arithmetic floor, the `||A||` estimate it rests on, and **all
seven Krylov methods** — CG, MINRES, GMRES, FGMRES, LSQR, LSMR and BiCGSTAB — are in.

`solve()`, `eigs()` and `svds()` stay unexported, and the dispatch they need is fully
determined: `positive_definite` selects CG, `hermitian` selects MINRES, and GMRES requires
nothing, so it is the fallback that makes `method = "auto"` total rather than partial. What
is missing is the dispatch itself and a deliberate edit to `BUDGET` in `test-api-budget.R`.
Promoting them is a decision, not a blocked one.

`eigs()` and `svds()` are a larger gap than the roster suggests: plan section 7.2 is
unstarted, there is no eigensolver code in `R/` at all, and Phase 2 ships as v0.1 with them
in it. The two remaining Gate 2 items are the benchmark harness with committed results and a
solver vignette.

Two documents govern:

- `implementation_plan.md` — the design, in numbered sections. Cite section numbers.
  Several sections record a reversal of an earlier decision and why; before proposing an
  alternative, check whether the plan already rejected it and on what grounds.
- `dev_notes/` — what was actually measured, including **seven corrections to the plan**
  from the spikes. Where the two disagree, `dev_notes/` wins: it has executed evidence and
  the plan does not. `dev_notes/GATE1.md` lists the spike corrections in one table;
  `dev_notes/cg-and-the-arithmetic-floor.md`,
  `dev_notes/fgmres-and-preconditioner-sides.md`,
  `dev_notes/minres-and-the-preconditioned-norm.md`,
  `dev_notes/gmres-and-the-second-pass.md`,
  `dev_notes/lsqr-and-the-least-squares-certificate.md`,
  `dev_notes/lsmr-and-the-monotone-backward-error.md` and
  `dev_notes/bicgstab-and-the-recurrence-with-nothing-to-lose.md` carry the Phase 2 ones.

## Commands

```r
devtools::load_all(".")                  # after any R/ change
roxygen2::roxygenise(".", clean = TRUE)  # regenerate NAMESPACE and man/
devtools::test(".")                      # ~2 min, 11,245 assertions
```

```powershell
R CMD build . ; R CMD check --as-cran --no-manual linop_0.0.0.9000.tar.gz
```

`R CMD check` should end at `Status: 1 NOTE` (new submission / development version number).
Anything else is a regression.

The expression fuzzer defaults to 400 trees so routine checks stay inside CRAN's time
budget. Gate 1 requires 10,000:

```r
Sys.setenv(LINOP_FUZZ_N = "10000"); testthat::test_file("tests/testthat/test-fuzzer.R")
```

Full gate verification, including install cost: `dev_notes/spikes/gate1.R`.

Use the Windows R at `C:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe`. Write `.R` files and
run them; never `Rscript -e` with complex code.

## Settled by the spikes — do not re-litigate

| Decision | Why | Source |
|---|---|---|
| **S3, with a `matrixOps` group method** | One method covers `%*%`, `crossprod`, `tcrossprod`; `.Generic` tells them apart | S0.1 §2 |
| **Not S7** | S7 0.2.2 cannot register a method on `crossprod`: `S7:::group_generics()` hardcodes `matrixOps <- "%*%"`, predating R 4.4.0. Worked around only by bypassing S7's own registrar | S0.1 §3 |
| **`Depends: R (>= 4.4.0)`** | `crossprod`/`tcrossprod` became S3 generic in 4.4.0, not 4.3.0 as the plan assumed | S0.1 §1 |
| **Zero `Imports`, no compiled code** | Callback overhead is 200 ns/apply, under 0.13% of an apply at n >= 1e5. The C fast path contingency is not needed | S0.2, GATE1 |
| **`dim<-`, never `as.matrix()`, on the apply path** | 300 ns vs 2300 ns | S0.2 §1 |

## The correction that bites most often

**Base R's `crossprod(A)` is `A^T A`, not `A^H A`.** For complex `A` it is symmetric and
*not* Hermitian. Plan sections 5.4 and 5.6 have this backwards; the code follows base R,
because the package's whole promise is that a `linop` behaves like a matrix, and the fuzzer
compares against base R.

The Hermitian Gram matrix is `adjoint(A) %*% A`. The product node **recognises**
`adjoint(X) %*% X` structurally (`propagate.R:is_adjoint_of`) and stamps
`hermitian = TRUE, source = "construction"` with an empty `depends_on`, because it holds
unconditionally. `test-propagation.R` asserts both halves.

## Architecture

`R/` in dependency order: `aaa-utils` → `evidence`, `caps`, `dtype` → `node-registry` →
`core` → `leaves`, `nodes` → `propagate` → `simplify` → `linop`, `algebra` →
`materialize`, `print`, `norm` → `certificate` → `verify` → `provenance`, `linsolve`,
`preconditioner` → `solvers-common` → `solvers-cg`, `solvers-minres`, `solvers-gmres`,
`solvers-bicgstab`, `solvers-bidiag` → `solvers-lsqr`, `solvers-lsmr` → `zzz`.

`certificate.R` owns the certificate object: `build_certificate()`, its print method and
`solve_certificate()` all live there, and `verify.R` holds only the operator checks.
`solvers-common.R` owns what every method shares: `KRYLOV_CONDITION_LIMIT`,
`solver_setup()`, which validates the operator, the block, the budget and the starting
iterate and takes `square` as a parameter because a least-squares method needs the same
checks on a different shape, and the two run-time preconditioner refusals
(`split_precond_not_hpd()`, `left_precond_singular()`) that any method implementing all
three sides needs. A new solver adds a `solvers-*.R` file and nothing else.

`solvers-bidiag.R` is the one exception, and it is the FGMRES pattern at a coarser grain:
LSQR and LSMR build the same Golub-Kahan bidiagonalisation and differ only in what they
minimise over it, so the loop, the start and the step live once and each method supplies a
recurrence. `bidiag_solve()` takes that recurrence as an argument. What is *not* shared is
the projection, which is the method.

**Apply has four modes**, BLAS-style: `N` (`A X`), `T` (`A^T X`), `C` (`A^H X`), `R`
(`conj(A) X`), composed through `MODE_COMPOSE` (Klein four-group). This is what keeps
transpose, adjoint and conjugate as three distinct things per section 5.4. A `fun` leaf
supplies only `apply` and `adjoint`; `T` and `R` are derived by conjugation.

**Two-tier apply.** Authors implement tier 1 (`apply(op, X, mode)`); solvers only ever call
tier 2, `linop_gemm(A, X, Y, alpha, beta, mode)`, synthesised from tier 1. Keeps scratch
reuse uniform without demanding a GEMM signature from an adapter author.

**Evidence is three independent fields** (`source`, `guarantee`, `confidence`), never a
strength ranking, and `evidence_satisfies()` recurses into `depends_on`. The case this
exists to prevent has a test: a sum of two `user_declaration` positive-definite operators
must fail CG's requirement at any depth (`test-evidence.R`, "the laundering case").

**Capability values are three-valued.** `NA` means unknown and never reads as `FALSE`.
Propagation is conservative: a rule that cannot prove `TRUE` returns `NA`, never `FALSE`.

**Views are not rewritten at construction.** `t()` and `adjoint()` stay distinct nodes even
for real operators, so `print()` shows what the author wrote. The simplification happens at
apply time. Exception: leaves that absorb a view (`identity`, `diag`) collapse it, which is
itself a section 5.8 rule.

**The sparse leaf is real-only.** Matrix 1.7.5 declares the `zMatrix` virtual classes but
not the concrete ones. A real sparse operator applied to a *complex block* does work, by
splitting into real and imaginary parts (`leaves.R:sparse_apply`).

## Working rules

**The API budget is a test.** `test-api-budget.R` asserts the exported set exactly and
asserts Phase 2 names are absent. Adding an export fails it until `BUDGET` is edited, which
is the point. Never edit the budget as a side effect of adding a function.

**The deferred node types are asserted absent.** `lowrank`, `kron`, stacks, `blockdiag`,
`perm`, `power`, `inverse` are Phase 3, after an external adapter has exercised the
registry. A test fails if one appears. This is the mechanism, not a preference.

**Tests are recovery and contract tests, not shape tests.** 11,245 assertions across 279
tests. The propagation suite alone is 9,892 brute-force soundness checks. When adding a
solver, the bar is parameter recovery against closed-form truth run to convergence, plus
certificate coverage over >= 20 seeds (plan section 10). `helper-linop.R` carries the
section 10 fixtures with their closed forms: `laplacian_1d()` with
`laplacian_1d_eigenvalues()`, `kms_matrix()` with `kms_inverse()`, `spd_prescribed()` /
`hpd_prescribed()` for a dialled-in spectrum, `shifted_laplacian_1d()` with
`shifted_laplacian_solve()` for an indefinite system whose eigenvalues *and* eigenvectors
are both closed form, and `convdiff_1d()` with `convdiff_1d_solve()` for a **nonsymmetric**
one. The last carries a documented accuracy budget: its similarity is `diag(rho^j)`, so the
closed form is only truth while `rho^n` stays below about 1e4, and the table in the helper
says where. It stops being ground truth long before the operator stops being well
conditioned.

The rectangular pair is `lsq_prescribed()`, `A = U diag(sigma) V^H` with the solution, the
residual and `kappa` all closed form and both compatible and incompatible right-hand sides
reachable, and `diff_1d()`, matrix free and rank deficient, whose minimum-norm solution is
closed form. Both are checked against their own definitions before any solver runs on them.

**A test that compares two configurations must give them the same budget, and must check
the knob did what is claimed.** Three drafts have now failed this. The MINRES suite claimed
a preconditioner rescued a solve, on a scalar preconditioner producing bitwise identical
iterates, passing only because the two sides ran to different `maxit`. The GMRES suite
claimed reorthogonalisation improved accuracy; over 12 seeds it is worse on 4 of them, and
the draft passed on the one seed where the ratio is 13.5. The same suite compared a
preconditioner on an operator scaled so badly that the arithmetic floor swallowed the
tolerance and *both* sides certified as met. Before asserting a knob helped: check it
changed the iterates, check the claim survives several seeds, and require `pass` rather than
accepting `qualified`.

A related trap the LSQR floor test walks into: **a test whose threshold has to sit inside a
one-`c eps` window cannot use a constant.** A fixed `tol` lands inside or outside that window
depending on the fixture, so the test measures what the solve achieved and sets `tol` from
it. A constant would have passed or failed for reasons unrelated to the floor.

**A closed-form fixture is code, and wrong closed forms pass plausible tests.**
`convdiff_1d_eigenvalues()` was first written with `+2 sqrt(bc) cos(k pi/(n+1))` where the
sign belongs to `b rho`. The eigenvalue *set* is invariant under that error, because
`cos(k pi/(n+1))` over `k = 1..n` is symmetric about zero, so a check against a sorted
spectrum passed while every eigenvalue was paired with the wrong eigenvector and the solve
was wrong by O(1). Check a decomposition by `A S - S Lambda`, never by comparing sorted
spectra.

**Certificates carry an arithmetic floor, and it is load-bearing.** S0.6 found that the a
posteriori residual bound keeps decaying past machine epsilon while the true error plateaus
at ~1e-15, so a converged result certifies as `fail` without one. `solve_certificate()`
adds `c * eps * (||A|| ||x|| + ||b||)`, which is Higham's bound on the computed residual and
also the Rigal-Gaches denominator, so one term serves both the residual and backward-error
lines. A test asserts that `floor_const = 0` turns a converged solve's certificate from
`qualified` to `fail` on byte-identical iterates. A line that meets its tolerance only
through the floor is `qualified` with an `estimate` guarantee, never a clean `identity`.
The least-squares line reduces to the same `c * eps`, plus `floor_abs / ||r||`, which matters
only where `||r||` has fallen to the residual floor and the compatible reading is in use.

**`||A||` is a lower bound on purpose.** Every route in `norm2()` returns one, which shrinks
the Rigal-Gaches denominator, so a reported backward error can only overstate. Do not
"improve" it into an upper bound without re-reading which direction each certificate line
needs. Structural routes are exact; a structural rule over an estimated child records
`construction <- [computation/estimate]`, so `evidence_satisfies()` still sees the estimate
at the top. That is the laundering case of section 5.3 outside capabilities.

## Landscape — read before writing any comparative copy

`dev_notes/eigencore-audit.md` and `dev_notes/S0.5-prior-art.md` supersede plan section 2.

- **eigencore 1.0.2** (CRAN, 2026-07-25) already has a matrix-free block operator with the
  GEMM apply signature, an `as_operator()` S3 adapter generic, `adjoint()`, complex dense,
  block Krylov, shift-invert and certificates. Several gaps the plan claims are closed.
  What linop genuinely adds, verified by full read: evidence-bearing capabilities
  (eigencore's structure descriptor is one field with two values), preconditioner contract
  enforcement, linear solvers (eigencore has none), and R's matrix generics.
  eigencore's API is the constructor-object DSL plan section 1.1 rejects — state that as a
  deliberate tradeoff, never as a deficiency.
- **`bbuchsbaum/amatrix`** (GitHub, active) is backend-agnostic *stored* matrix classes
  with GPU backends. Complementary, not competing: it is about where a matrix lives, linop
  is about matrices that never exist. Its backend registry/planning/explain is prior art
  for plan section 3 — read it before designing that protocol.
- **`sanic`** (CRAN) is a Krylov linear-solve comparator the plan misses.
- MINRES, LSQR, LSMR confirmed absent from CRAN as user-callable solvers. Phrase it that
  way; "no R implementation exists" is not supported.

**Nothing about eigencore's `certification.R` or `validation.R` is quotable** — they were
not read. Read them in full before writing any certificate copy for a README or paper.

**No competitive vocabulary in any public artefact** (plan section 15). Report both numbers
side by side and let readers compare.

## Phase 2 entry points

`solve()`, `eigs()`, `svds()` are unexported. `linsolve` and `preconditioner` exist
internally with their contracts and enforcement already tested; `solver()` stays private
until inexact shift-invert or a PRIMME warm-start workflow creates a real need (plan
section 1.1), and that defer is only safe while `preconditioner()` is public.

`cg_solve()` (`R/solvers-cg.R`) is the template the other six follow, and `minres_solve()`,
`gmres_solve()`, `lsqr_solve()` and `lsmr_solve()` are the evidence that it generalises: to a
minimal-residual method, to a long recurrence with restarts, to an operator with no declared
capability at all, to a problem that is a minimisation rather than an equation, and to a
second method on that same problem. Three things CG established that the others inherit
rather than re-decide:

- **Several right-hand sides run in lockstep**, not one after another. Each column's
  recurrence is independent, so the iterates are exactly those of per-column CG (asserted
  bitwise) at one block apply per step instead of `k`. This is not block CG.
- **The outer loop measures, the inner loop trusts.** The recurrence residual drifts worst
  where the answer is most converged, so the decision to stop is taken on a recomputed
  `b - A x` and the certificate reports that one.
- **Non-convergence comes back, it is not thrown.** A stalled solve returns a `fail`
  certificate naming the exhausted budget. A *contradicted declaration* does throw, naming
  the capability: `p^H A p <= 0` far from convergence means the operator is not what it
  said it was.

`certificate.R` owns the certificate and `preconditioner.R` owns `precond_applier()` and
`precond_adjoint_applier()`, the single guards every solver applies `M^-1` and `M^-H`
through. A fifth copy of that closure in a new `solvers-*.R` is the thing to not write.

MINRES carried all three and forced two refinements, both in
`dev_notes/minres-and-the-preconditioned-norm.md`:

- **The two loops can measure different quantities.** MINRES minimises `||r||_{M^-1}`, so
  with a preconditioner the recurrence scalar is not the certificate's residual. The outer
  loop converts the target into the recurrence's currency going in and re-measures in the
  caller's coming out. Any method whose natural residual is preconditioned inherits this.
- **A contradiction test has to be on a quantity rounding does not move.** CG thresholds a
  sign. MINRES has to threshold a magnitude, and the identity the recurrence already carries
  (`v_{j-1}^H A v_j = beta_j`) is a test of Lanczos orthogonality, not of symmetry: it
  reaches a relative violation of 1e0 on *correct* hermitian operators with clustered
  spectra. The test that works is the definition, `<x, A y> = <A x, y>`, which needs no extra
  apply and stays at 1e-10. Measure the noise floor on correct operators before choosing any
  such threshold.

GMRES and FGMRES (`R/solvers-gmres.R`) are one implementation, and the four findings are in
`dev_notes/gmres-and-the-second-pass.md`:

- **GMRES requires nothing of the operator, so it has no declaration to contradict.** It is
  the only method that works on an operator supplying no adjoint: every apply is mode `N`.
  That also makes it the fallback that turns `method = "auto"` from partial into total.
- **`side` selects an algorithm here, not a label.** CG's row is unrestricted because the
  three sides *coincide*; GMRES's is unrestricted because all three are implemented and they
  produce different iterates. Right minimises the reported residual; left and split do not
  and inherit MINRES's currency conversion. Split is reachable from `M^-1` alone by MINRES's
  change of variable, carrying `u_j = M v_j` alongside `v_j`. A test asserts the three
  disagree, so "accepts all three sides" cannot decay into one path relabelled.
- **Complex arithmetic stops being optional.** `col_dot()` is `Re(<x,y>)`, which is the whole
  product only for the quantities a hermitian recurrence forms. Arnoldi's `<v_i, A v_j>` is
  genuinely complex, so this file uses `col_cdot()` and a rotation with real `cs` and complex
  `sn`. The three inner products stay separate: routing the hermitian methods through a
  complex product to discard half of it would cost them their real arithmetic.
- **Unconditional reorthogonalisation is worse than none**, by three orders of magnitude on
  the worst fixture, because a second pass over an already-cancelled direction adds noise to
  the Hessenberg column. The conditional (Daniel-Gragg-Kaufman-Stewart) criterion is the one
  to use, and it is per column, so masking it exactly keeps the lockstep identity — the
  reason the first draft rejected it was wrong.
- **A Krylov space keeps admitting directions after they stop meaning anything.** Past that
  point the recurrence residual falls monotonically while the true residual diverges, and
  the projected problem cannot tell. `KRYLOV_CONDITION_LIMIT` stops there.
  `max|R_ii|/min|R_ii|` is a lower bound on `cond(R)`, so it stops later and never earlier,
  and it is bitwise inert on every well-posed fixture tried. LSQR reaches the same `1/eps`
  by a second route: the singular values of its bidiagonal lie inside those of `A`.

LSQR (`R/solvers-lsqr.R`) is the first method whose problem is a minimisation rather than an
equation, and the four findings are in `dev_notes/lsqr-and-the-least-squares-certificate.md`:

- **A least-squares solve has no exact answer, and the certificate had a square-system
  shape.** `b - A x` does not go to zero, so testing it against `tol` reports a converged
  solution as a failure — the S0.6 shape again, structural rather than arithmetic. What goes
  to zero is `A^H r`, and Stewart's perturbation `dA = -(r r^H A)/||r||^2` makes that a
  backward error rather than a heuristic: it is *exhibited*, so `||A^H r||/(||A|| ||r||)` is
  an achievable relative perturbation, and it is also Paige and Saunders' second stopping
  rule. **Both readings are backward errors, so one row carries both** and the certificate
  keeps the same shape whatever method produced it. Which reading applies is decided by
  measurement inside `solve_certificate(least_squares = TRUE)`, never from a solver's flag.
  The floor works out to `c eps` on that line too, and a test asserts `floor_const = 0` turns
  a converged least-squares solve from `qualified` to `fail` on byte-identical iterates.
- **Nothing the LSQR recurrence believes about itself is reportable.** Two arithmetically
  equivalent implementations diverge by three orders of magnitude every two steps, reaching
  O(1) relative around step 12 of a 15-column *well-conditioned* problem, then re-converge.
  The bidiagonal vectors lose orthogonality and rounding noise fills the lost direction. A
  reference test can assert bitwise agreement for CG; here it can assert four steps.
- **Finite termination does not survive floating point.** At step `n` the answer is wrong in
  the second digit; it needs about `1.7 n`. Never assert termination at `n`.
- **The preconditioner rows narrow to `right`, and the contract had a gap.** `M^-1` acts on
  the domain; left and split act on the codomain, which for a rectangular operator is a
  different space, and where the spaces coincide they change the minimiser rather than the
  path. The adjoint of `A M^-1` is `M^-H A^H`, which no earlier method needed, so
  `preconditioner()` takes `apply_inverse_adjoint` and `hermitian = TRUE` supplies it for
  nothing. Neither declared is refused by name, never inferred.

LSMR (`R/solvers-lsmr.R`) shares that bidiagonalisation and minimises a different thing over
it, and the three findings are in `dev_notes/lsmr-and-the-monotone-backward-error.md`:

- **It minimises the number the certificate reports, which is what earns it a method rather
  than a flag.** LSQR minimises `||r||`; LSMR minimises `||A^H r||`, the numerator of the
  least-squares backward error. Measured densely from the iterates over 20 steps, 12 seeds
  and three conditionings, LSMR's reported quantity rose *zero* times and LSQR's rose on
  every seed, by up to a factor of 60 in one step. At a shared budget on the same fixture and
  seed it reports the smaller backward error on 12 of 12 seeds.
- **The estimate a method minimises is still not one it may report.** `|zetabar|` is
  `||A^H r||` exactly for the projected problem, and by step 30 of a well-conditioned
  15-column solve it reads 1.3e-16 against a true 1.1e-13: wrong by 800x, in the direction
  that certifies a solve as better than it is. The advantage is real and the recurrence's
  measurement of it is not, so the outer loop's extra apply is what makes the advantage
  reportable at all. Nothing in the monotonicity claim above was measured from the
  recurrence.
- **The reproducibility ceiling is the method's, not the code's.** LSQR's was measured
  against a second implementation written here; LSMR's was measured against SciPy 1.17.1,
  and diverges on the same schedule (1e-15 at four steps, 1e-8 at eight, O(1) at sixteen,
  re-converging). Four steps is what a reference test can assert.

BiCGSTAB (`R/solvers-bicgstab.R`) is the last of the seven, and the three findings are in
`dev_notes/bicgstab-and-the-recurrence-with-nothing-to-lose.md`:

- **It agrees with the published recurrence bitwise, at every step count tried**, on real,
  nonsymmetric and complex fixtures, where LSQR and LSMR manage four steps. That fact is
  about the other six methods: what drifts in a Krylov method is orthogonality, and
  orthogonality is a property of a stored basis. BiCGSTAB orthogonalises nothing and asks no
  vector to stay orthogonal to any other, so there is no invariant to lose. It is not a
  claim of accuracy; it is the same fact as the next one, read the other way.
- **What it gives up is the minimisation.** GMRES's true residual cannot rise; BiCGSTAB's
  rose on 4 to 10 of 20 steps on the nonsymmetric fixtures, by up to a factor of 41 in one
  step. The strongest case in the package for the outer loop measuring.
- **Breakdown comes back as a certificate.** The cure for a breakdown is re-seeding the
  shadow vector, which is what the outer loop's next round already does, so it costs no new
  machinery. A cure that recovers nothing is reported: on a real skew-symmetric operator,
  where `<A z, z> = 0` makes the alpha breakdown exact at every first step, the result is a
  `fail` certificate naming the breakdown, a finite iterate and no error thrown, while GMRES
  solves the same system to 3.8e-16. The two methods that require nothing of the operator
  are not interchangeable.

BiCGSTAB is the second method that runs without an adjoint, so `method = "auto"` has two
candidates for an operator supplying only a forward action, not one.

`solve_certificate()`'s `forward error: not_checked` got its first live demonstration here.
At `kappa` 1e10 with an incompatible right-hand side, four runs certifying backward errors
between 6e-9 and 2.4e-8 show forward errors from 7.5e-5 to 9.2e-1, because `kappa^2 eps` is
2.2e4 there and no backward error constrains the forward one. Do not read that spread as one
method being more accurate than another.

The evidence minima (`CG_PD_REQUIREMENT`, `MINRES_HERMITIAN_REQUIREMENT`) filter
`method = "auto"` and do not gate `method = "cg"` or `method = "minres"`. Naming a method is
the caller asserting their own declaration; `auto` is the package choosing. Do not collapse
the two. GMRES, BiCGSTAB, LSQR and LSMR have no such minimum because they have no
requirement.

`method = "auto"` is now total over rectangular shapes as well as square ones: LSQR and LSMR
are the two methods that accept one.

For `linop.primme` (Phase 3), S0.3 established that PRIMME builds under Rtools45 and on
macOS arm64, and that Windows R lacks `zheevx_`/`zhegvx_` so a shim over `zheevd_` is
needed. The spike shim is a proof of concept with documented gaps, not shippable.
