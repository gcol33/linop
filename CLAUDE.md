# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository state

Phase 0 (spikes) and Phase 1 (core) are complete and Gate 1 is met. **Phase 2 is complete
and Gate 2 is met**: the solve certificate with its arithmetic floor, the `||A||` estimate it
rests on, **all seven Krylov methods** — CG, MINRES, GMRES, FGMRES, LSQR, LSMR and
BiCGSTAB — and section 7.2's two spectral verbs are in.

`solve()` is an S3 method on the base generic and cost no export. `eigs()` and `svds()` are
new names and cost one each, which is the whole of what Phase 2 spends on the public surface.
The benchmark harness runs end to end with committed results in `dev_notes/bench/`. The
reference-agreement line is met for both methods the gate names: LSMR against SciPy 1.17.1,
MINRES against two references at once over kappa 1e2 to 1e10 with breakdown and
near-breakdown constructed rather than hoped for.

The articles list is seven, `pkgdown::check_pkgdown()` clean. Phase 2 added three:
`solvers` (`solve()`, dispatch, the seven methods, preconditioner sides), `spectral`
(`eigs()`, `svds()`, `ncv`, shift-invert) and `certificates`. The split is deliberate and
the reason is the certificate: it takes **four shapes** (`verify()`'s operator conformance,
the square system, least squares, the eigenpair), and the differences between them are the
content, so one article holds them side by side rather than each verb explaining its own.
eigencore made the same call for the same reason; KrylovKit.jl separates by problem class.
All seven knit in 3.9 s total, so vignette build cost is not a constraint on example size.

Two things section 7.2 lists for v0.1 are deliberately not here, and
`dev_notes/eigs-svds-and-the-third-certificate.md` records why: RSpectra delegation (deferred
to Phase 3 with the section 3 backend registry, since wiring it in first means an ad-hoc
branch the registry then has to absorb) and Arnoldi, so a non-hermitian operator has no
eigensolver. The generalized problem `A x = lambda B x` is refused by name.

`dev_notes/rspectra-and-the-delegation-that-is-not-a-superset.md` measures what that
delegation could carry. RSpectra closes the non-hermitian gap matrix-free and does `svds`
and `which = "SM"`, and it refuses a block apply, a complex dense matrix and `sigma` on a
function, which `eigs()` supplies. **A complex operator behind a callback is not refused:
it is coerced, and the returned spectrum is `Re(A)`'s**, matching `eigen(Re(A))` to
2.842e-14 and wrong by 6.946 against the truth. Any wrapper has to refuse on linop's own
dtype before the call. The note also puts RSpectra ahead of `linop.primme` as the first
backend, since it exercises the registry with no compiled code.

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
  `dev_notes/minres-and-the-breakdown-that-does-not-end.md`,
  `dev_notes/gmres-and-the-second-pass.md`,
  `dev_notes/lsqr-and-the-least-squares-certificate.md`,
  `dev_notes/lsmr-and-the-monotone-backward-error.md`,
  `dev_notes/bicgstab-and-the-recurrence-with-nothing-to-lose.md`,
  `dev_notes/solve-dispatch-and-the-empty-branch.md`,
  `dev_notes/eigs-svds-and-the-third-certificate.md`,
  `dev_notes/benchmarks-and-the-fixture-that-solved-itself.md` and
  `dev_notes/compile-ceiling-and-the-basis-that-was-copied.md` carry the Phase 2 ones.

## Commands

```r
devtools::load_all(".")                  # after any R/ change
roxygen2::roxygenise(".", clean = TRUE)  # regenerate NAMESPACE and man/
devtools::test(".")                      # ~3 min, 11,538 assertions
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
| **Zero `Imports`, no compiled code** | Callback overhead is 200 ns/apply, under 0.13% of an apply at n >= 1e5, so the C fast path contingency is not needed **on the apply path**. What that figure does not cover, and the eigensolver's own ceiling, is `dev_notes/compile-ceiling-and-the-basis-that-was-copied.md` | S0.2, GATE1 |
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
`solvers-bicgstab`, `solvers-bidiag` → `solvers-lsqr`, `solvers-lsmr` → `solve` →
`eigen-common` → `eigs`, `svds` → `zzz`.

`certificate.R` owns the certificate object: `build_certificate()`, its print method and
`solve_certificate()` all live there, and `verify.R` holds only the operator checks.
`solvers-common.R` owns what every method shares: `KRYLOV_CONDITION_LIMIT`, `REORTH_ETA`
(the Daniel-Gragg-Kaufman-Stewart second-pass criterion, used by GMRES and by both
eigensolvers), `solver_setup()`, which validates the operator, the block, the budget and the starting
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

**Tests are recovery and contract tests, not shape tests.** 11,976 assertions across 340
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

**A profiler's self time is not a compile ceiling, because a subscript is billed to the
callee.** R forces a promise at its first use, so `f(V[, seq_len(i), drop = FALSE], W)`
allocates and copies inside `f`, and `Rprof` charges `f`'s own frame for a subscript the
caller wrote. `orth_against()` reads 1.132 s of self time inside `eigs()` on the largest
cell measured, the single largest entry in the run, and **0.0%** profiled on its own, where
its two BLAS products account for 96%. Timed against those products it runs at 1.04x, so
there is nothing in it to compile. Before quoting any cell of
`dev_notes/spikes/results/eigs-compilable-ceiling*.csv`, read
`dev_notes/compile-ceiling-and-the-basis-that-was-copied.md`: the figures are real
measurements of a real cost, and that cost is the basis copy rather than R's interpreter.

The ceiling is **1.01x to 1.12x for the six solvers** and **1.18x to 1.70x, median 1.36x,
for the eigensolver alone**. Never state the second as the package's. Pure R cannot close
it: the one bitwise-identical route (orthogonalise against the whole preallocated basis,
whose unused columns are zero) is measured at 0.98x, because the copy it removes and the
extra BLAS it adds are the same quantity averaged over a round. What closes it is a GEMM
with a leading dimension, and the loop it sits in is the one `linop.primme` deletes rather
than compiles, so the finding argues for Phase 3 rather than for `src/`.

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

`linsolve` and `preconditioner` exist internally with their contracts and enforcement
already tested; `solver()` stays private
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

Gate 2's reference line added two more, in
`dev_notes/minres-and-the-breakdown-that-does-not-end.md`:

- **An exact Lanczos breakdown is one collapsed off-diagonal, not the end of the sequence,
  and stopping there returns the worse answer.** With `b` inside a `d`-dimensional invariant
  subspace, `beta_{d+1}/||A v_d||` collapses to between 8.3e-16 (`d = 1`) and 9.6e-11
  (`d = 8`), and the *next* entry comes back at O(1), because rounding in `A v` re-seeds the
  space. Neither reference stops at `d` either, and SciPy's minres does not. Iterating past
  the breakdown improves the iterate: at `d = 8` the exact-arithmetic termination step is
  wrong by 3.9e-12 and eight more steps bring that to 9.7e-15. So nothing thresholds `beta`
  to detect a breakdown and nothing should be added; the residual test is the stopping rule.
  `spent` in `minres_recurrence()` is not breakdown detection, it separates an exhausted `w`
  from a preconditioner that contradicted its declaration, and it suppresses an error rather
  than a step.
- **The reproducibility ceiling belongs to the method, and here that is provable rather than
  suspected.** The suite carries two references that share no code and no idea: the
  definition (`krylov_argmin`, a dense minimiser over an orthonormalised basis) and the
  published Paige-Saunders recurrence (`reference_minres`, validated against SciPy 1.17.1 in
  the spike, not in the suite, so there is no Python dependency). Run against *each other* at
  kappa 1e6 they agree to 1.4e-15 at four steps, 2.8e-13 at eight, 2.0e-8 at twelve and not
  at all at sixteen, which is the same schedule each of them follows against this
  implementation. Eight steps is what a reference test can assert, where LSQR and LSMR manage
  four. Do not try to extend it by tightening anything.

GMRES and FGMRES (`R/solvers-gmres.R`) are one implementation, and the four findings are in
`dev_notes/gmres-and-the-second-pass.md`:

- **GMRES requires nothing of the operator, so it has no declaration to contradict.** That
  is what makes it the fallback turning `method = "auto"` from partial into total. Every
  apply is mode `N`, so it also runs on an operator supplying no adjoint, but that part is
  not particular to it: measured against an adjoint-less operator, CG, MINRES, FGMRES and
  BiCGSTAB all run too, and **only LSQR and LSMR require an adjoint**. Requiring nothing and
  needing no adjoint are different properties and only the first singles GMRES out.
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

BiCGSTAB is the second method that requires nothing of the operator, so `method = "auto"`
has two candidates for an operator that establishes no capability, not one. Both also run
without an adjoint, which every square method in the roster does.

`solve_certificate()`'s `forward error: not_checked` got its first live demonstration here.
At `kappa` 1e10 with an incompatible right-hand side, four runs certifying backward errors
between 6e-9 and 2.4e-8 show forward errors from 7.5e-5 to 9.2e-1, because `kappa^2 eps` is
2.2e4 there and no backward error constrains the forward one. Do not read that spread as one
method being more accurate than another.

`solve()` (`R/solve.R`) is the one verb over all seven, and the findings are in
`dev_notes/solve-dispatch-and-the-empty-branch.md`:

- **It cost no export.** `solve` is a base generic, so the method reaches it the way `t()`
  and `%*%` do and `BUDGET` did not move. Section 1.1's "keeps the verb budget" is literal.
  `eigs()` and `svds()` will not be free: new names, so the budget grows by two.
- **The `cg` branch of `auto` is correct and has no supply.** Measured over every constructor
  in the package, nothing reaches it: `linop_scaling()` and `linop_eye()` establish
  definiteness at an acceptable source and are diagonal, so the direct route takes them
  first, and every other route to `positive_definite` is a bare user declaration, which
  `auto` will not act on. A dense symmetric leaf gets MINRES because symmetry is *computed*
  and definiteness is not. Closing this is a design decision, not a fix: either a Cholesky
  attempt on a dense hermitian leaf (O(n^3) against the existing O(n^2) entry checks, so it
  needs a gate) or the product node stamping `positive_definite` on `adjoint(X) %*% X` when
  `X` is established full column rank. A test asserts the branch is empty and says to replace
  it, not adjust it, when that changes.
- **The certificate attribute survives arithmetic**, which plan section 1.1 has backwards.
  `x + 0`, `2 * x`, `x + x`, `t(x)` all keep it; indexing, `as.numeric()`, `c()` and
  reductions drop it. So the real cost is not losing the certificate, it is that a
  certificate can outlive the value it describes. The decision stands, since a classed result
  would break the promise that a solve returns what a matrix solve returns, but the cost is
  the sharper one.

The evidence minima (`CG_PD_REQUIREMENT`, `MINRES_HERMITIAN_REQUIREMENT`) filter
`method = "auto"` and do not gate `method = "cg"` or `method = "minres"`. Naming a method is
the caller asserting their own declaration; `auto` is the package choosing. Do not collapse
the two. GMRES, BiCGSTAB, LSQR and LSMR have no such minimum because they have no
requirement.

**`method = "auto"` is total over square shapes and refuses rectangular ones.** That is
section 1.1 and it is deliberate: least squares is a different mathematical request rather
than the same one on a different shape, so `solve()` errors and names `lsqr` and `lsmr`,
which are reachable through the same verb. The roster being total over rectangular shapes and
`auto` being total over them are different statements, and only the first is true.

GMRES rather than BiCGSTAB is `auto`'s fallback. Both require nothing and run without an
adjoint; GMRES's residual cannot rise and it has no breakdown, which are the properties to
prefer when the package is choosing rather than the caller.

## The spectral verbs

`eigs()` (`R/eigs.R`) and `svds()` (`R/svds.R`) are section 7.2, with `R/eigen-common.R`
holding what they share. The findings are in
`dev_notes/eigs-svds-and-the-third-certificate.md`:

- **The eigenpair certificate is a third shape, not a flag on the second.** Both readings of
  `solve_certificate()` treat `b` as data, so Rigal-Gaches divides by `||A|| ||x|| + ||b||`.
  For an eigenpair `theta` and `x` are both outputs, `A` is the only datum, and the
  denominator has no second term. The perturbation is exhibited as Stewart's is: a Ritz value
  is the Rayleigh quotient of its own Ritz vector, so `x^H r = 0` and
  `E = -r x^H - x r^H` is hermitian with `(A + E) x = theta x` and `||E||_2 = ||r||`.
- **`forward error` is a real bound here, the first in the package.** `A + E` is hermitian
  and `theta` is exactly its eigenvalue, so Weyl gives
  `min_j |theta - lambda_j(A)| <= ||r||/||x||`. That makes this the first row carrying
  `guarantee = "deterministic_bound"` and `source = "theorem"`, and it exposed a latent bug:
  `build_certificate()` listed "no deterministic bound on" by testing `guarantee != "identity"`,
  which was indistinguishable from correct while `identity` and `estimate` were the only two
  values in use. It now tests membership in `DETERMINISTIC_GUARANTEES`. The bound is about
  *some* eigenvalue; which one is `target identity` and stays `not_checked`.
- **A certificate row can carry evidence with `depends_on`.** `cert_rows()$add()` takes an
  `evidence()` object instead of the three flat fields, and the forward row records the
  hermitian capability's evidence under it. So a bound resting on a bare `user_declaration`
  fails a requirement the declaration would have failed directly — section 5.3's laundering
  case reaching the certificate. `cert$evidence[["forward error"]]` is where it lives; the
  printed table keeps the flat fields.
- **`eigs()` applies no evidence minimum, and that is the decision rather than an omission.**
  Every other dispatcher that applies one has a fallback; here there is one method and no
  non-hermitian eigensolver until Arnoldi, so a minimum would only refuse to run, on exactly
  the callback operators the package exists for. The value gates the run, the evidence is
  reported. Do not "fix" this by adding a requirement constant.
- **The svds certificate is the eigs certificate, on the augmented operator.**
  `H = [[0, A], [A^H, 0]]` is hermitian for every `A` with nothing to declare, and its
  eigenpairs are `(sigma, [u; v]/sqrt 2)`. So one builder serves both, and a singular value's
  forward bound rests on *no* declaration where an eigenvalue's rests on whatever established
  the operator's symmetry. `||H|| = ||A||` exactly, so the norm is inherited structurally.
  The one weaker line: `orthogonality` measures the augmented basis, where a deviation in `U`
  could cancel one in `V`; `test-svds.R` asserts the two separately.
- **Thick restarting is what makes the methods converge, not a refinement.** Restarting from
  the sum of the unconverged Ritz vectors stalls: on `laplacian_1d(60)` at `ncv = 24` asking
  for four pairs it spends 240 of 300 iterations over 9 rounds and converges none, at a
  backward error of 5.7e-4, where the thick restart reaches 5.1e-13 in 96. Where one round
  suffices the two agree to the last digit, which is the control that confines the difference
  to the swapped block. `dev_notes/spikes/restart-comparison.R` runs both.
- **The Rayleigh quotient is reported, never the Ritz value.** It minimises `||A x - mu x||`
  over `mu` and is measured on `A`, so under shift-invert it corrects what the inner solve got
  wrong instead of inheriting it. That is what lets `sigma` be built out of the package's own
  solvers (`A - sigma I` through MINRES) with nothing downstream trusting the inner tolerance.
- **Reference means storage, not accuracy.** A round holds its whole basis and orthogonalises
  against all of it: O(ncv) vectors, O(ncv^2 n) work. Full reorthogonalisation is block
  classical Gram-Schmidt with the `REORTH_ETA` second pass. The loop is *not* shared with
  GMRES and must not be — GMRES keeps the coefficients because they are its Hessenberg, here
  they are a correction and are discarded, and keeping them would make it Arnoldi. The same
  distinction separates `svds()` from `solvers-bidiag.R`, which builds the same recurrence
  and stores no basis.

## The benchmark harness

`dev_notes/bench/`, run with `Rscript dev_notes/bench/run-bench.R` from the repository root.
2.5 minutes, four CSVs, an environment stamp and a generated `RESULTS.md`, all committed.
`dev_notes/benchmarks-and-the-fixture-that-solved-itself.md` has the findings.

- **A right-hand side never draws from the seed its fixture drew from.** Both prescribed
  generators build their operator from a QR factorisation, and the first draws from a seed are
  the first column of the input, which QR puts in the span of `Q[, 1]`. Sharing the seed makes
  `b` the first eigenvector of `spd_prescribed` (CG converges in **one** iteration) and puts it
  exactly in the range of `lsq_prescribed`, so the incompatible fixture is compatible. Neither
  announces itself: both return the right answer and certify correctly. A fixture whose point
  is a property of `b` now measures that property and carries it in the results.
- **Applies is the primary unit and seconds is second**, because an apply count is the same on
  every machine. It does not price orthogonalisation, so both columns are reported and neither
  is called the cost. Memory is `object.size()` only: `gc()`'s max used counts allocation
  between collections and reads 112 Mb for a solve holding kilobytes.
- **The counted operators declare their capabilities rather than computing them**, since an
  operator with a counter in it is a callback leaf. So the tables do not measure `auto` and
  every row names its method.
- **`ncv` is the binding knob for both spectral verbs.** On `laplacian_1d(400)` at `ncv = 40`,
  five times the budget and a tolerance two orders looser give the same iteration count, the
  same restarts and the same value error; `ncv = 80` converges 6 of 6 to 1.1e-16. The stall
  detector stopping at a quarter of the budget is reporting that, not giving up early.
  `dev_notes/spikes/eigs-ncv-probe.R`.
- **The sparse route wins the size ladder at every size**, and that number is committed. A
  matrix-free method is the route to an operator with no matrix, not a faster route to a banded
  one. No cross-package comparison is in the harness, for the section 15 reasons in the note.

For `linop.primme` (Phase 3), S0.3 established that PRIMME builds under Rtools45 and on
macOS arm64, and that Windows R lacks `zheevx_`/`zhegvx_` so a shim over `zheevd_` is
needed. The spike shim is a proof of concept with documented gaps, not shippable.

## What comes next, and in what order

`dev_notes/hilbert-first-and-the-envelope-that-does-not-dispatch.md` records the decision
and what is already in place. **The Hilbert layer's first unit comes before `linop.primme`
and before the adapters**, which reverses the plan's Phase 5 slot. Neither is a
prerequisite: the unit is self-adjoint, so it wants `eigs()` and the certificate and not
Arnoldi, and it couples through the provenance envelope, which is exported and in `BUDGET`.
It also has a spike behind it where PRIMME has no caller. `dev_notes/S0.6-finite-section-bound.md`
returned all three certificate quantities closed form and verified, and the arithmetic
floor every Phase 2 certificate now carries came out of it.

The first unit is exactly S0.6's class: `FiniteSection` on `ell^2(Z)`, self-adjoint Jacobi
plus finite-support `V`, the three-part certificate, and refusal at `V = 0` on `q - a > 0`.
Everything else in plan section 8 stays out, and that is what the separate package is for.
The package split is not what moved; only the order is.

**The four provenance generics dispatch on `p$payload`, not on the envelope.**
`set_provenance()` stores an unclassed `list(provider =, payload =)`, so dispatching on the
envelope reaches nothing but the defaults: a provider's method was unreachable and all four
`.default` messages named the wrong cause. The payload's class survives and belongs to the
provider, so `UseMethod(generic, p$payload)` fixes it with no signature change, no new
export and a byte-identical `NAMESPACE`. Handing a class to `UseMethod()` is not inspecting
the payload, so section 5.11 holds.

The test that should have caught it was already there and passed:
`test-provenance.R`'s "a provider can register methods and core will route to them"
hand-built its envelope with `structure()` instead of obtaining one from
`set_provenance()`, so it exercised the generics and not core — the half its title names
was the half not run. **A fixture for an object the caller never builds by hand is not a
test of the path the caller takes.** It now goes through `set_provenance()` and covers all
four generics; against the old generics it fails. `dev_notes/spikes/provenance-dispatch-probe.R`.

**The first unit has been built against v0.1 core and the operator path costs no export.**
`dev_notes/hilbert-first-unit-and-the-certificate-a-provider-cannot-build.md` has the run:
a matrix-free `finite_section(V, n)`, `verify()` passing eleven checks, all four provenance
generics routing including through `t(H) %*% H`, and S0.6's table reproduced through
`eigs()` rather than dense `eigen()`. Three things came out of it:

- **The certificate was the one thing a provider could not build, and now it is exported.**
  The Hilbert certificate is a fourth *shape* — truncation, isolation, arithmetic floor, no
  residual row and no backward-error row — but the row table, the evidence fields, the
  `overall` roll-up and the print method are the same object all four shapes share, so the
  alternative was a second copy of the roll-up in every satellite. `build_certificate()` and
  `cert_rows()` are exported; `BUDGET` was edited deliberately and the surface now stands at
  exactly the 32 `test-api-budget.R` allows, so the next export fails that test. Both gained
  the validation a public entry point needs: `CERT_STATUSES` is closed because the roll-up
  reads those four strings exactly and a fifth would be counted as a pass, and `$add()`
  checks `source` and `guarantee` against the same vocabulary `evidence()` uses.
  `test-certificate.R` is new and asserts the contract an outside caller sees.
- **A matrix-free provider cannot say `construction`.** `properties=` stamps
  `ev_declared()` unconditionally, so the finite section — hermitian because it is
  tridiagonal with unit off-diagonals — reports `user_declaration`, while the *dense* leaf
  of the same operator gets `computation` for free. `eigs()` carries that into the
  forward-error row's `depends_on`, so a provable Weyl bound reports as resting on a
  declaration. The smallest route out is a `properties=` form accepting `capability()`
  objects, which needs no new export since `capability()` and `evidence()` are both public.
- **The eigensolver is a second signal for the `V = 0` refusal and not a substitute.** At
  `n = 10` and `40` `eigs()` converges cleanly and returns a value still inside the band; at
  `n = 100` it stalls, and that stall is clustering near the band edge rather than absence of
  an eigenvalue. Only `q - a > 0` separates the two.

**The satellite exists.** `~/dev/linop.hilbert`, its own repository, `linop.hilbert`,
first commit `9e48e98`: `finite_section()`, `discretise()`, `eigenpairs()`,
`decay_rate()`, `sites()`, 319 assertions, `R CMD check` at one WARNING whose two halves
are linop not being on CRAN yet and the GitHub repo not existing yet. It uses no `:::`.
Its own `dev_notes/the-first-unit-and-the-subspace-that-was-too-narrow.md` carries three
findings, and the first is about **this** package:

- **`eigs()`'s default `ncv` is too narrow for a near-band eigenvalue.** On
  `finite_section(0.3)`, whose eigenvalue is 2.0224 against a band edge of 2, `ncv = 21`
  (the default at `k = 1`) converges *nothing* at `n = 80` and lands 7.3e-04 from the
  truth, where `n = 40` reaches 4.2e-07 — a wider block giving a worse answer, with only
  `nconv` saying so. `ncv = 40` reaches 2.7e-12 and `ncv = 80` and `160` change no digit.
  This is the benchmark note's "`ncv` is the binding knob" on a fixture that is not
  `laplacian_1d`, and it argues the default deserves a second look here rather than only
  in the satellite.
- The satellite's certificate needed a fifth row for that, as a **qualification and never
  a failure**, because `dist(theta, sigma(H)) <= ||(H - theta)u||/||u||` holds for every
  `theta` and every nonzero `u`. A stalled inner solve gives a true bound, a worse one.
- **The arithmetic floor is load-bearing only once the inner solve is good enough.** At
  `ncv = 21` the truncation residual plateaus three orders above the true error and the
  crossing never happens; at the satellite's default it happens at `n = 80`, 7.734e-17
  against a true error of 8.882e-16, which is S0.6's row exactly. A test that fixed `n`
  would be measuring the eigensolver rather than the floor.

Two repos and the name `linop.hilbert` are settled, with the reasoning and the one
constraint the name carries (**never a class called `hilbert`**, since `linop()` is a
generic whose adapter convention is `linop.<class>()`) in
`dev_notes/hilbert-first-and-the-envelope-that-does-not-dispatch.md`. **`linop` goes to
CRAN before the satellite**, or the satellite needs `Additional_repositories:` and a live
`PACKAGES` file. Core must never gain a `Suggests` on the layer; the dependency runs one
direction.

Backend order is RSpectra then PRIMME, for the reasons in
`dev_notes/rspectra-and-the-delegation-that-is-not-a-superset.md`.

Undecided and not blocking: one repo with two package directories against two repos (a
tooling question, not a design one), and the name `linop.hilbert` (plan open question 2).
