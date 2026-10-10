# PRIMME, and the basis that is bounded rather than removed

Issue #1 (2026-07-30) asked for `linop.primme` against the section 3 registry, on
the premise that PRIMME's implicit restart "never forms the 610 MB basis and never
runs the orthogonalisation", so that the 1.42x of
`compile-ceiling-and-the-basis-that-was-copied.md` would be the floor of what the
backend gives. The engine is built, linked against R's own LAPACK on Windows,
routed through an internal registry, and certified exactly as the reference is.
The premise is half right: the basis is bounded, the orthogonalisation is not
removed, and on the fixture that converges the reference is as fast or faster.

Probes, all single process on the i9-14900K with R 4.6.1's reference BLAS:
`dev_notes/spikes/primme-restart-probe.R`, `primme-variants-probe.R`,
`primme-block-probe.R`, `primme-agreement-probe.R`, results beside them.

## 1. Where it lives

Not on CRAN as `PRIMME`: that package was archived on 2026-06-07 (maintainer
address undeliverable), so it can be neither imported nor suggested from a
mainstream repository. The plan's route B stands: vendor the C library.

The one-package decision (`one-package-and-the-abstractions-that-were-holding-a-boundary.md`)
says capability never leaves `linop` and only a compiled engine is external. So
`linop.primme` is an engine and nothing else: a sibling package at
`~/Documents/dev/linop.primme`, private repository `gcol33/linop.primme`, with
three exports (`primme_eigsh()`, `primme_svds()`, `primme_engine()`) that only
`linop` calls. Users write `eigs(A, method = "primme")`. Core stays at zero
compiled code; the registry entry names the package as data and loads it with
`requireNamespace()` when it is asked about.

`linop`'s DESCRIPTION does not list it yet. A `Suggests` entry needs a repository
serving the package (`Additional_repositories`, the gcol33 r-universe) and the
repository is private, so listing it now would break the public CI's dependency
resolution. That is the step when it is published, not before.

## 2. The LAPACK gap closes without the expert drivers

S0.3 found Windows R lacks `zheevx_` and `zhegvx_` and proved a shim over
`zheevd_` with documented gaps. Reading PRIMME's `blaslapack.c` before patching,
as S0.3 asked, settles it differently. PRIMME has compile switches `USE_XHEEV` and
`USE_XHEGV` that make its dense kernels call the simple drivers `zheev`/`dsyev`
and `zhegv`/`dsygv`. Against R 4.6.1's `Rlapack.dll` (846 exports) every one of
those is present except `zhegv`, and its building block `zhegst` is missing too.

So `src/lapack_ext.c` supplies one routine, `zhegv` for `itype = 1`, by the
reference reduction: `zpotrf`, two `ztrsm` for `U^-H A U^-1`, `zheev`, one `ztrsm`
back. `itype = 1` is the only type PRIMME requests (both call sites pass the
constant 1, and the Hermitian projected problem is the only caller). The routine
is linked under a private name, `-Dzhegv_=lp_zhegv`, so it never interposes on a
full LAPACK's `zhegv` elsewhere. Tested against `solve(B, A)` for both triangles
with the unused triangle poisoned with `NaN`, B-orthonormality of the vectors, and
`info = n + i` for an indefinite `B`. It is exercised end to end because a block
size above one makes PRIMME orthogonalise explicitly and solve its projected
problem against `V^H V`, which is the generalized driver.

The S0.3 shim is retired, gaps included: nothing calls `zheevx` or `zhegvx`.

`R CMD check --as-cran` on the engine: no compiled-code NOTE (no stdout, no
`abort` from PRIMME). One trap on the way: `roxygen2::roxygenise()` compiles
through pkgbuild with `-UNDEBUG`, and a following `R CMD INSTALL` reused those
objects, so the installed DLL carried PRIMME's `assert()`s and an operator error
became an abort. The CLAUDE.md rule on stale objects applies to flags as well as
headers: delete the objects after any `load_all()` before an install that is
measured or shipped.

## 3. The registry

`R/backends.R` is section 3's protocol, internal. Each entry carries `name`,
`family` (`eigs_hermitian`, `svds`), `priority`, `contract_version`,
`engine_version`, `package`, `dtypes`, `requires`, `targets`, `supports()` and
`run()`. `supports()` returns `TRUE` or a reason, never `FALSE`. Selection: a
named `method` takes that engine or errors with its reason; `"auto"` takes the
eligible engine of highest priority, where `options(linop.backend_priority =)`
overrides. An external engine is eligible only when installed and reporting the
same contract version. Every result's `certificate$dispatch$candidates` lists
each engine with `selected` / `eligible` / `rejected` and the reason, which is
what section 3 asked `plan_eigs()` to show; no export was spent.

PRIMME's `supports()` rejects, with the reason recorded and `auto` falling back
to the reference:

- `tol` below machine epsilon. PRIMME refuses it (code -12), and reads `tol = 0`
  as a request for its own default of 1e4 eps, which is not what `tol = 0` means
  here.
- `svds()` with `v0`: PRIMME starts a singular value problem from both sides.
- an `ncv` too small for a block of `k` (section 5).

## 4. The certificate survives the swap

What crosses from the engine is vectors. The reported eigenvalue is the Rayleigh
quotient on `A`, measured in the adapter; under `sigma` PRIMME iterates on
`(A - sigma I)^-1` and its own value is never used. A singular value is
`Re(u^H A v)` of the unit vectors, the Rayleigh quotient of the augmented vector.
Residuals, convergence and the forward bound are the certificate's, from the same
`eigen_certificate()`, so the forward row's `depends_on` still carries the
hermitian evidence: a bare declaration fails a strict requirement through PRIMME
exactly as through Lanczos (`test-backends.R`).

PRIMME's own stopping test is the backward error alone, `||r|| < eps ||A||`, and
a pair stopping exactly there fails the certificate's residual row whenever
`|value| < ||A||`: 4.0e-10 against 1e-10 on the well. The engine's
`convTestFun` is replaced with the certificate's two rows, each with its floor
`c eps (||A|| + |value|)`; without the floor the relative row is unreachable for
a small eigenvalue and the run exhausts its budget.

Under `sigma` the replacement test runs on the inverse, which bounds the residual
on `A` by `tol ||A - sigma I||` rather than by `tol |mu|` (1.13e-10 against
1e-10 measured). And the reference has the mirror problem: it locks on the
backward error alone and fails the residual row on `laplacian_1d(80)`'s small end
(7.7e-10). Both are #2; fixing the reference moves measured baselines.

A budget that converges nothing in `svds()` returns no triplet: PRIMME's singular
value solver packs only converged triplets. The adapter used to read `k` columns
at a fixed offset and returned garbage there; it now unpacks `nconv`, and the
certificate has an empty shape whose convergence row fails "0 of 4". Non-convergence
still comes back rather than being thrown.

## 5. Two defaults the measurement chose

**The preset is GD+k, not PRIMME's default.** `PRIMME_DYNAMIC` switches between
GD+k and JDQMR on measured timings: on the well at n = 1e5 two identical calls
took 180 and 149 applies (`primme-variants-probe.csv`), and six identical calls
on `laplacian_1d(2000)` returned six different sets of vectors in 8869 to 10000
applies. Section 3 puts reproducibility ahead of the last few percent, so the
engine defaults to `PRIMME_DEFAULT_MIN_MATVECS`, which is bitwise reproducible
(asserted in the engine's tests).

**The block is k.** A method grown from one vector reaches one vector of each
eigenspace in exact arithmetic and only rounding supplies the rest. Asked for the
`k` copies of a `k`-fold eigenvalue (orders 30, 60 and 200, multiplicities 2, 3,
5 and 8, ten seeds each, `primme-block-probe.R`):

| | 2- and 3-fold | 5-fold | 8-fold, order 30 | 8-fold, order 60 | 8-fold, order 200 |
|---|---|---|---|---|---|
| PRIMME GD+k, block k | all | all | 8 | 8 | 8 |
| PRIMME GD+k, block 1 | all | all | 6 or 7 | 7 | 8 |
| reference Lanczos | all | 4 at order 30, else all | 3 to 5 | 6 | 8 |

A block of `k` is the only configuration that returned every copy every time. The
reference's random restart is not the guarantee its test reads as: it fails at
order 30 from five copies up. The certificate never claimed otherwise, since
`target identity` is `not_checked`, but `test-eigs.R`'s triple eigenvalue passes
because three is below where it starts to miss. Under PRIMME's DYNAMIC preset a
block of one missed a triple on six of six seeds in an earlier run of the same
probe (its log was overwritten when the preset was fixed), which
is the other reason the preset is fixed.

The price, on the well at n = 1e6: 32.6 s and 166 applies against 15.1 s and 121
for a block of one (`primme-variants-probe.R`). Correctness was taken over the
2.2x; it is the knob to revisit if a caller can promise simple eigenvalues.

A block of `k` needs a basis of at least `2k + 1` (4 for one eigenpair at an
interior target), measured for `k` = 1 to 8 at every target of both solvers;
below it PRIMME's own parameter check fails, so `supports()` rejects it.

## 6. What the restart removes

`primme-restart-probe.R`, through `eigs()`, certificate included, `k = 4`.
`budget`: the compile-ceiling stencil at equal work, 400 applies, nothing
converges. `well`: the stencil minus a well, four bound states to 1e-10.

| fixture | n | engine | seconds | applies | peak MB |
|---|---|---|---|---|---|
| budget | 1e6 | lanczos, ncv 20 | 104.9 | 400 | 1191 |
| budget | 1e6 | lanczos, ncv 80 | 235.9 | 400 | 3529 |
| budget | 1e6 | primme, default (basis 20) | 71.5 | 400 | 712 |
| budget | 1e6 | primme, ncv 80 | 186.3 | 400 | 1636 |
| well | 1e6 | lanczos, ncv 20 | 16.0 | 69 | 1175 |
| well | 1e6 | lanczos, ncv 80 | 28.8 | 80 | 1708 |
| well | 1e6 | primme, default | 33.5 | 166 | 720 |
| well | 1e6 | primme, ncv 80 | 79.8 | 162 | 1636 |

n = 1e5 shows the same ordering (CSV).

**The basis is bounded.** PRIMME's default holds 15 to 20 vectors whatever the
problem (20 with a block of 4), and peak memory at n = 1e6 is 712 to 720 MB
against 1175 MB at `ncv = 20` and 1708 to 3529 MB at `ncv = 80`. That half of the
premise holds.

**The orthogonalisation is not removed.** PRIMME orthogonalises every new vector
against its basis, 35 to 50 inner products per apply at its default size and 148
to 167 at a basis of 80, and on top of that recomputes Ritz vectors, `A` times
them and residuals every iteration rather than at a restart: 58 `n` dense flops
per apply for GD+k at a block of one and 125 to 129 `n` at a block of four
(`primme-variants-probe.R`). Orthogonalisation is 23 to 30% of PRIMME's own time
at its default basis and 34 to 43% at 80; the dense updates are 33 to 51% of it
for GD+k. At equal work and equal basis
(ncv 80) PRIMME is 1.27x the reference; at each engine's default, 1.47x against
`ncv = 20` and 3.3x against `ncv = 80`. The compiled `lda` the issue asked for is
inside PRIMME, and it is one of its costs rather than the removal of one.

**Time to an answer is the reference's.** On the well the reference converges in
69 applies at `ncv = 20`; GD+k needs 121 at a block of one and 166 at a block of
four. The fastest deterministic PRIMME configuration, timed on the engine alone,
takes 15.1 s against the reference's 16.0 s through `eigs()` with its certificate
(about 0.5 s of that is the `||A||` probe), so call it a tie; the shipped one,
block `k`, takes 2.1x as long. JDQMR, the other deterministic preset, takes
17.5 s at a block of one and 24.6 s at four.

So the 1.42x is not the floor of what the backend gives, and "deletes the loop
the copy is in" was wrong. What PRIMME gives at these sizes is the bounded basis
(memory), a basis size that does not have to be tuned (`ncv = 80` was needed by
the reference on hard spectra, `SECTION_NCV_MIN` exists for the same reason, and
on `finite_section(linop_jacobi(diagonal = 0.3), 80)` at a basis of 21 PRIMME's
eigensolve meets its target, 2.7e-12 from the eigenvalue of the operator, where
the reference's stalls at 9.8e-3; the certificate still fails both, correctly,
on the truncation term at that width),
determinism under its GD+k preset, and the engine for what the reference cannot
do yet (preconditioning, the generalized problem, interior targets without
shift-invert). It is not a speedup on a cheap matrix-free operator with R's
reference BLAS, and nothing here measures an optimised BLAS, where PRIMME's dense
kernels would gain and the reference's `%*%` with them.

## 7. Agreement with the reference

With PRIMME's Arnoldi preset, a block of one, a minimum restart size of one, a
basis wider than the run and its arithmetic floor at zero, its basis is the
Krylov space of the starting vector, which is what the reference's first round
builds. Started from the same unit vector, the Rayleigh quotients after `m`
applies agree to 4.7e-16 relative at every `m` from 1 to 40, on four fixtures
(real kappa 1e2 and 1e6, complex kappa 1e4, `laplacian_1d(400)`) and five
seeds, with no drift in `m`. That is BiCGSTAB's row of the roster, not LSQR's
four steps, and for the reason the BiCGSTAB note gives: both engines keep a fully
orthogonalised basis, so there is no orthogonality to lose. The test asserts
1e-14 at 1, 2, 4, 8, 16 and 32 steps and that PRIMME spent exactly `m` applies.

Under the shipped preset the iterates are a different method's and only the
answers agree, which `test-backends.R` asserts against closed-form spectra for
both engines.

## 8. The compile-ceiling spike, re-run

`eigs-compilable-ceiling.R` resolves to the reference (it passes `tol = 0`, which
PRIMME's `supports()` rejects, and the test helper's priority option is loaded by
`load_all()` in any case). Re-run 2026-10-10: 1.14x to 1.49x, median 1.30x,
against 1.18x to 1.70x, median 1.36x, with the stencil n = 1e6, ncv = 80 cell at
1.47x against 1.57x. The loop is unchanged, so this is run-to-run spread on a
shared machine, and the 1.42x projection for a compiled GEMM in the reference
still sits below its cell's ceiling. What is corrected is what the projection
was said to bound, section 6.

## 9. Corrections

| Where | Was | Is |
|---|---|---|
| Issue #1, compile-ceiling note section 5 | PRIMME "never forms the 610 MB basis and never runs the orthogonalisation"; the 1.42x is the floor of what it gives | It bounds the basis (712 MB peak at n = 1e6 against 1175 to 3529) and keeps the orthogonalisation, plus dense updates every iteration. 1.27x to 3.3x at equal work; to an answer 0.48x as shipped and a tie at best |
| S0.3, plan 13.3 | Windows needs `zheevx`/`zhegvx` shimmed over `zheevd` | `USE_XHEEV`/`USE_XHEGV` route PRIMME to the simple drivers; only `zhegv` is missing, supplied for `itype = 1`, reference-exact by construction |
| CLAUDE.md | PRIMME's default method | Not reproducible run to run; the engine uses GD+k |
| `print.linop_eigen` | printed "section n = <nconv>" on every non-section result (`x$n` partial-matched `nconv`) | fixed, `39e5763` |
| `SECTION_NCV_MIN` | a property of finite sections | a property of the reference's eigensolve: PRIMME at a basis of 21 meets the eigensolve target the reference misses by 9.8e-3 |

## 10. Not done

- macOS arm64 and Linux builds of `linop.primme`. S0.3 built PRIMME on macOS,
  not this package; Gate 3 needs all three.
- `Suggests: linop.primme` with `Additional_repositories`, at publication.
- #2, the stopping rules against the residual row.
- Preconditioning and `B` through PRIMME. The registry entry declares neither.
