# FGMRES is the seventh method, and `side` is now enforced

Decided 2026-07-27, entering Phase 2. Two changes to the preconditioner model and one
correction to the plan.

## 1. FGMRES ships in v0.1

Plan section 7.1 listed six Krylov methods and omitted FGMRES. That table was the outlier,
not a decision. FGMRES was already committed in five places:

| Location | What it says |
|---|---|
| `implementation_plan.md` section 4.3 | an FGMRES row in the enforcement table, `fixed` may be `FALSE` |
| `implementation_plan.md` section 16 | Saad on FGMRES cited as the source of the `fixed = FALSE` rules |
| `vignettes/adapters.Rmd` | shipped copy: fgmres is "the method designed to accept it", and the table "covers `cg`, `minres`, `gmres`, `fgmres`, `bicgstab`, `lsqr` and `lsmr`" |
| `R/preconditioner.R` | the `fgmres` row, and CG's refusal message pointing at it |
| `test-linsolve-preconditioner.R` | two Gate 1 assertions |

Removing it instead would have meant editing shipped vignette copy, deleting a plan row,
orphaning a cited reference and dropping two passing tests.

Two further reasons, either of which would carry it on its own:

**It is the only v0.1 consumer of `variable_inexact`.** `accepts_contract()`
(`R/linsolve.R:44`) has three branches. Shift-invert `eigs()` accepts `exact` and
`fixed_linear`. Without FGMRES the inexact half of the section 4.2 fidelity lattice ships
in v0.1 with nothing that reads it until PRIMME in Phase 3.

**It is a flag on GMRES, not a second solver.** The Arnoldi loop, the Givens rotations and
the residual recurrence are unchanged, because `A Z_m = V_{m+1} H_bar_m` holds by
construction. The difference is storing the preconditioned basis `Z` and forming the update
as `x = x0 + Z_m y`, at the cost of one further `n x m` array. The recurrence is to be
written against Saad section 9.4 rather than from recall.

## 2. `side` was stored and never read

`check_preconditioner()` looped over required boolean flag names only. `side` was set by
the constructor, printed by `print.preconditioner()`, and consulted by nothing. A
`preconditioner(fixed = FALSE, side = "left")` handed to FGMRES would have passed.

`PRECOND_REQUIREMENTS` now carries `flags` and `sides` per method, and
`check_preconditioner()` checks both. Changed while the table has no callers; after seven
solvers consume it, the same change is a migration.

## 3. Which side rows are sourced

Only FGMRES is restricted:

- **FGMRES, `right` only.** It forms its update from the right-preconditioned basis. PETSc
  `KSPFGMRES` states the same restriction: "Only right preconditioning is supported."
  (read 2026-07-27, `https://petsc.org/release/manualpages/KSP/KSPFGMRES/`).
- **CG, MINRES, GMRES, BiCGSTAB, LSQR, LSMR: unrestricted, pending a sourced pass.** Each
  admits more than one standard formulation (left in the `M`-inner product, split via
  `M = L L^H`, right), and the correct per-method set is a question for Saad, *Iterative
  Methods for Sparse Linear Systems*, chapter 9, not for recall. Accepting all three sides
  is the conservative reading: it can only fail to refuse, never wrongly refuse. A test
  asserts this explicitly so that narrowing a row is a deliberate edit that changes a
  passing test.

This is the open item. Narrowing those six rows needs the book, and the rows should not be
guessed at from memory in the meantime.

## 4. `linsolve` -> `preconditioner`

An inner solve run to a loose tolerance is the canonical flexible preconditioner, and there
was no path from a `linsolve` into section 4.3: `as_preconditioner()` took a `linop` plus a
solver function.

`as_preconditioner()` is now S3, with methods for `linop`, `linsolve` and a `default` that
refuses. The exported name count is unchanged, so the API budget is untouched.

The `linsolve` method reads `fixed` off the contract on all three axes
(`fidelity`, `determinacy`, `randomness`), mirroring `verify.linsolve()`, which holds a
solve to repeatability under exactly that condition. `side` defaults to `right` for a
flexible solve, since FGMRES is its only consumer, and to `left` otherwise.

`hermitian` and `positive_definite` stay `NA`. A linsolve declares a fidelity and a
determinacy, never a symmetry, and an exact solve of an SPD operator is not the same claim
as a preconditioner declaring itself SPD. CG refusing the result by name is the
three-valued rule working rather than a gap; a caller who can justify the claim passes it
through explicitly.

## Corrections to the plan

| Plan location | Correction |
|---|---|
| 7.1 | seven Krylov methods, not six; FGMRES was missing from this table only |
| 4.3 | `side` is enforced, and the table carries a side column; FGMRES is `right` only |
| 9, Gate 2 | the per-row refusal test covers sides; FGMRES with a fixed right preconditioner must reproduce right-preconditioned GMRES to machine precision |
