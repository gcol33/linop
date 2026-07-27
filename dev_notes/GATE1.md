# Gate 1 (Phase 1 -> Phase 2)

Verified 2026-07-27 on R 4.6.0 (2026-04-24 ucrt), x86_64-w64-mingw32.
Reproduce with `dev_notes/spikes/gate1.R`.

## Status: met

| Gate 1 requirement (plan section 9) | Result |
|---|---|
| Conformance suite passes on every node type, real and complex | **pass**, with the amendment below |
| Expression fuzzer, 1e4 trees vs base R, 1e-12 relative | **pass**, 10,000 trees, 0 failures, 6.3 s |
| Flag propagation checked by brute force, incl. `A^H A` vs `A^T A` on complex | **pass**, 9,892 assertions |
| Evidence propagation: `depends_on` recorded, `evidence_satisfies()` recurses | **pass** |
| `linsolve` conformance separate; `linop` suite against a `linsolve` errors | **pass** |
| API budget enforced by a test | **pass** |
| dtype promotion table exhaustive | **pass**, all 4 combinations |
| Install cost recorded | **zero hard dependencies** |
| `R CMD check --as-cran` | 1 NOTE, expected (see below) |

## Numbers

```
files: 10  tests: 76  assertions: 10214  failed: 0  errors: 0  skipped: 0  time: 26.6s

                           file passed failed error
    test-algebra-and-simplify.R     51      0     0
              test-api-budget.R     13      0     0
             test-conformance.R    114      0     0
                   test-dtype.R     24      0     0
                test-evidence.R     30      0     0
                  test-fuzzer.R      8      0     0
 test-linsolve-preconditioner.R     37      0     0
             test-propagation.R   9892      0     0
              test-provenance.R     22      0     0
   test-registry-and-adapters.R     23      0     0
```

Fuzzer at the full Gate 1 count: **10,000 trees, 0 failures, 6.3 s.** The committed
test defaults to 400 trees so routine `R CMD check` stays inside CRAN's time budget;
`LINOP_FUZZ_N=10000` runs the gate count and is what the number above reports. Only the
coverage of the space changes with the count, not what is checked.

## Install cost

```
Depends:   R (>= 4.4.0)
Imports:   (none)
LinkingTo: (none)
Suggests:  Matrix, testthat (>= 3.0.0), knitr, rmarkdown
recursive hard deps: 0
compiled code: none
```

This is stronger than the plan's target. Section 3 specifies `Imports: methods`; with S3
throughout and no S4 classes, `methods` is not needed either. Section 14's secondary
criterion "install cost is `methods` only" is met with room to spare.

Matrix is in `Suggests`, not `Imports`, so the sparse leaf is available when Matrix is
present and costs nothing when it is not.

## R CMD check

`Status: 1 NOTE`:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Gilles Colling <gilles.colling051@gmail.com>'
New submission
Version contains large components (0.0.0.9000)
```

Expected for a development version, and this is not a submission; Phase 2 ships v0.1.
Everything else is OK, including examples, tests, S3 registration, Rd, and cross-references.

## Amendment carried from S0.1

"Conformance passes on every node type, **real and complex**" cannot hold literally:
Matrix 1.7.5 declares the `zMatrix` virtual classes but not the concrete ones, so a complex
sparse operator cannot be constructed at all. Restated and satisfied as: *every node type
real, and complex on every node whose backing storage admits it, which is all but
`sparse`.* A complex operator with sparse structure goes through a `fun` leaf.

A real sparse operator applied to a **complex block** does work: `sparse_apply()` splits
into real and imaginary parts and recombines rather than downcasting, so section 5.5's "no
path silently downcasts" holds there too. This was found by the conformance suite, not by
inspection.

## What the tests actually assert

Not shape and plumbing. The suite is:

- **Conformance (114).** All nine section 5.10 checks against every node type; plus
  negative controls that a wrong adjoint, a nonlinear operator, an inconsistent block
  apply, an impure operator and a contradicted declaration are each *caught*.
- **Propagation (9,892).** 300 random operator/expression combinations over 12 expression
  shapes, each asserting that no `TRUE` flag contradicts the materialised matrix. Includes
  the corrected `crossprod(A) = A^T A` rows and the structural `adjoint(A) %*% A`
  recognition.
- **Evidence (30).** Including the laundering case: a sum of two `user_declaration`
  positive-definite operators must fail CG's requirement, at every depth.
- **Fuzzer (10,000 trees).** Materialisation *and* action compared against base R.
- **API budget (13).** Exported set asserted exactly; Phase 2 names asserted absent.
- **Adapters (23).** A storage format the package has never seen (a circulant held as its
  first column) registered from outside, passing `verify()` and the whole algebra with no
  change to linop. This is the plan's primary acceptance criterion, exercised.

## Corrections to the plan made during Phase 1

All traceable to Phase 0 spikes. Full detail in the individual spike notes.

| Plan location | Correction | Source |
|---|---|---|
| 9, S0.1 | floor is R >= 4.4.0, not 4.3.0 | S0.1 section 1 |
| 5.4, 5.6 | `crossprod(A)` is `A^T A`: symmetric always, hermitian iff real. The Hermitian Gram is `adjoint(A) %*% A` | S0.1 section 6 |
| 1.1 | `nrow`/`ncol` are not generic; they come free from `dim` | S0.1 section 5 |
| 3 | core needs no `Imports` at all, not even `methods` | this gate |
| 5.2 | normalise vectors with `dim<-`, not `as.matrix()` | S0.2 section 1 |
| 13.2, 17.3 | callback overhead is not a risk; core stays pure R | S0.2 sections 2-3 |
| 9, Gate 1 | complex conformance cannot cover the sparse node | S0.1 section 7 |

## Not done in Phase 1, by design

`solve()`, `eigs()`, `svds()`, the six Krylov methods, and the certificate for *results*
(as opposed to operators) are Phase 2. `verify()` currently has methods for `linop` and
`linsolve`; the `verify(fit, A)` method arrives with the solvers.

The deferred node types (`lowrank`, `kron`, stacks, `blockdiag`, `perm`, `power`,
`inverse`) are Phase 3, and a test asserts they are absent, so they cannot be added from
inside before an external adapter has exercised the registry.
