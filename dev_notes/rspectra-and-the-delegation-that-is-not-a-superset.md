# RSpectra, and the delegation that is not a superset

Plan section 7.2 lists RSpectra delegation for v0.1;
`dev_notes/eigs-svds-and-the-third-certificate.md:188-192` defers it to Phase 3 with the
section 3 backend registry, and describes the interface it would delegate through as
"single vector and real". The probe is `dev_notes/spikes/rspectra-delegation-probe.R`,
with results in `results/rspectra-delegation.csv` and `.txt`. RSpectra 0.16.2, R 4.6.0,
x86_64-w64-mingw32.

Both halves of that description hold. What the probe adds is that they fail differently,
and that the delegation runs in both directions: three capabilities RSpectra has and
`eigs()` does not, and three `eigs()` has and RSpectra does not.

## 1. What it has that core does not

| case | result |
|---|---|
| real non-symmetric, matrix-free, n = 60, k = 3 | max abs diff 3.682e-12 against dense `eigen()` |
| `which = "SM"`, matrix-free | works, smallest 0.00833893 |
| `svds` matrix-free with `Atrans` | max abs diff 1.243e-14 against `svd()` |

The first is the gap `R/eigs.R:111` names by hand: `eigs()` refuses a non-hermitian
operator and says Arnoldi arrives with the backend. RSpectra is ARPACK, it takes a
`function(x, args)`, and it closes that gap through `Suggests` at zero install cost. It
also restarts implicitly, so it does not build the basis
`dev_notes/compile-ceiling-and-the-basis-that-was-copied.md` prices at 610 MB per call at
n = 1e6, ncv = 80.

The `svds` route wants the second action as a separate `Atrans` argument. linop has it as
mode `C` of one operator, so supplying it is a closure over `linop_gemm()` and costs
nothing structural.

## 2. What core has that it does not

| case | result |
|---|---|
| block apply | `60x1` on all 37 callback invocations of a run |
| complex dense | error: `REAL() can only be applied to a 'numeric', not a 'complex'` |
| `sigma` on a function | error: `unsupported matrix type` |

The third is the one that was not anticipated anywhere. `eigs()` does matrix-free
shift-invert, built on the package's own MINRES against `A - sigma I`, and reports the
Rayleigh quotient so nothing downstream trusts the inner tolerance. RSpectra takes `sigma`
and it works on a dense matrix; on a function it is refused outright. So delegation is not
a superset in either direction, which is what makes the registry's `supports()` returning a
rejection string rather than a bare `FALSE` load-bearing rather than tidy. There is a real
request `plan_eigs()` has to route back to the reference implementation after a backend is
installed.

## 3. The complex path answers a different operator

The interesting failure is the one that does not fail. Given a hermitian complex `A` of
order 40, built as `Ac + t(Conj(Ac))` and checked exactly symmetric, the callback path
returns three real numbers:

```
returned            : -16.560944   15.872478  -15.506388
eigenvalues of A    : -23.507308  -22.238910   22.172162
eigenvalues of Re(A): -16.560944   15.872478  -15.506388
```

Agreement with the spectrum of `Re(A)` is 2.842e-14. Disagreement with the spectrum of `A`
is 6.946. The callback returns a complex vector, RSpectra coerces it to double, the
imaginary part goes, and what comes back is the correct answer to a question nobody asked.
It raises 54 coercion warnings over the run and then returns a converged result with no
error.

"Cannot carry a complex operator" is therefore the right conclusion and the wrong
mechanism. A complex *dense* matrix is refused, cleanly and at the door. A complex operator
behind a callback is accepted, and the wrapper is the only thing that could tell. If linop
delegated on dtype without checking, `svds()` and `eigs()` would hand back a certificate
whose backward error is honestly computed against `Re(A)`, and every line of it would be
true of an operator the caller did not supply. The refusal has to be linop's, on the dtype
it already tracks, before the call rather than after it.

## 4. What this decides

**RSpectra is the first backend and PRIMME is the second.** RSpectra needs no vendoring, no
`Makevars`, no LAPACK shim and no `LinkingTo`; it is a CRAN package in `Suggests` behind a
`requireNamespace()` gate. It is a genuine external consumer of the registry, which is what
Phase 3 gates the deferred node types on, and it exercises the protocol at a fraction of
the cost of finding out through PRIMME whether the protocol is right.

It does not make PRIMME redundant. The three rows in section 2 are close to a description
of what PRIMME is for: block methods, complex arithmetic, preconditioned and shift-invert
spectral solves. `dev_notes/S0.3-primme-build.md` has it building under Rtools45 and on
macOS arm64 with the `zheevx_`/`zhegvx_` shim, so the order is a sequencing decision rather
than a reprioritisation.

If `Suggests: RSpectra` is added while `linop.primme` is not yet on CRAN, DESCRIPTION needs
`Additional_repositories:` for the second one. RSpectra itself is on CRAN and needs
nothing.

## 5. Corrections

| Where | Was | Is |
|---|---|---|
| `dev_notes/eigs-svds-and-the-third-certificate.md:189` | "single vector and real, so it can carry neither a complex operator nor a block apply" | Correct on both, by two different mechanisms. Block: never asked for, `60x1` throughout. Complex: refused when dense, silently answered against `Re(A)` when matrix-free |
| Delegation as a strict gain | A backend adds capability | Not here. Matrix-free shift-invert is refused by RSpectra and supplied by `eigs()`, so the registry has to route in both directions |
| Phase 3 backend order | `linop.primme` implied first, being the named backend | RSpectra first: no compiled code, CRAN `Suggests`, and it exercises the registry before PRIMME's build cost is spent |

## 6. Not done

- No delegation wrapper written. This probe establishes what one would have to refuse; it
  does not implement one.
- Accuracy compared against dense `eigen()`/`svd()` only, at n = 60 and m = 40. Nothing
  here measures RSpectra at the sizes where its implicit restart is the point.
- `opts` was left at its default throughout, so nothing is known about how its tolerance
  and `ncv` controls would map onto the certificate's.
