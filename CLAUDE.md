# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository state

Phase 0 (spikes) and Phase 1 (core) are complete and Gate 1 is met. Phase 2 is under way:
the solve certificate with its arithmetic floor, the `||A||` estimate it rests on, and CG,
the first of the seven Krylov methods, are in. MINRES, GMRES/FGMRES, BiCGSTAB, LSQR and
LSMR are not, and `solve()`, `eigs()` and `svds()` stay unexported until there is more than
one method for `method = "auto"` to choose between.

Two documents govern:

- `implementation_plan.md` — the design, in numbered sections. Cite section numbers.
  Several sections record a reversal of an earlier decision and why; before proposing an
  alternative, check whether the plan already rejected it and on what grounds.
- `dev_notes/` — what was actually measured, including **seven corrections to the plan**
  from the spikes. Where the two disagree, `dev_notes/` wins: it has executed evidence and
  the plan does not. `dev_notes/GATE1.md` lists the spike corrections in one table;
  `dev_notes/cg-and-the-arithmetic-floor.md` and
  `dev_notes/fgmres-and-preconditioner-sides.md` carry the Phase 2 ones.

## Commands

```r
devtools::load_all(".")                  # after any R/ change
roxygen2::roxygenise(".", clean = TRUE)  # regenerate NAMESPACE and man/
devtools::test(".")                      # ~30 s, 10,403 assertions
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
`preconditioner` → `solvers-cg` → `zzz`.

`certificate.R` owns the certificate object: `build_certificate()`, its print method and
`solve_certificate()` all live there, and `verify.R` holds only the operator checks. A
second solver adds a `solvers-*.R` file and nothing else.

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

**Tests are recovery and contract tests, not shape tests.** 10,403 assertions across 120
tests. The propagation suite alone is 9,892 brute-force soundness checks. When adding a
solver, the bar is parameter recovery against closed-form truth run to convergence, plus
certificate coverage over >= 20 seeds (plan section 10). `helper-linop.R` carries the
section 10 fixtures with their closed forms: `laplacian_1d()` with
`laplacian_1d_eigenvalues()`, `kms_matrix()` with `kms_inverse()`, and `spd_prescribed()` /
`hpd_prescribed()` for a dialled-in spectrum.

**Certificates carry an arithmetic floor, and it is load-bearing.** S0.6 found that the a
posteriori residual bound keeps decaying past machine epsilon while the true error plateaus
at ~1e-15, so a converged result certifies as `fail` without one. `solve_certificate()`
adds `c * eps * (||A|| ||x|| + ||b||)`, which is Higham's bound on the computed residual and
also the Rigal-Gaches denominator, so one term serves both the residual and backward-error
lines. A test asserts that `floor_const = 0` turns a converged solve's certificate from
`qualified` to `fail` on byte-identical iterates. A line that meets its tolerance only
through the floor is `qualified` with an `estimate` guarantee, never a clean `identity`.

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

`cg_solve()` (`R/solvers-cg.R`) is the template for the remaining six. Three things it
established that the others inherit rather than re-decide:

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

The evidence minimum (`CG_PD_REQUIREMENT`) filters `method = "auto"` and does not gate
`method = "cg"`. Naming a method is the caller asserting their own declaration; `auto` is
the package choosing. Do not collapse the two.

For `linop.primme` (Phase 3), S0.3 established that PRIMME builds under Rtools45 and on
macOS arm64, and that Windows R lacks `zheevx_`/`zhegvx_` so a shim over `zheevd_` is
needed. The spike shim is a proof of concept with documented gaps, not shippable.
