# Bringing the Hilbert layer forward, and the envelope that does not dispatch

Two things, recorded so a fresh session can start on the Hilbert unit without re-deriving
them.

**Sections 1, 2, 4 and 5 are a decision and an inventory, not a measurement.** They record
what was chosen and what already exists, and nothing in them was run. Section 3 is the only
part with a probe behind it, `dev_notes/spikes/provenance-dispatch-probe.R`, with output in
`results/provenance-dispatch.txt`. `dev_notes/` outranks the plan because it carries
executed evidence, so a note that mixes the two has to say which is which.

## 1. The decision

The Hilbert layer is wanted, and its first unit moves ahead of the plan's Phase 5 slot,
ahead of `linop.primme` and ahead of the adapters.

The package split itself is unchanged and is not the thing being reconsidered. Section 3's
rule stands: anything that costs an install second, or moves at research pace, goes in a
separate package. What moves is the order.

## 2. Why the first unit does not have to wait

**It has a spike behind it and PRIMME does not have a caller.**
`dev_notes/S0.6-finite-section-bound.md` already ran the plan's own risk control for this
layer, which risk 5 states as "first unit is one class with provable bounds". All three
quantities the certificate needs came back closed form and computable, verified
numerically: the decay rate `rho(lambda)` to 3.3e-15 against measured ratios, the
truncation bound `eta` from two entries of the computed eigenvector with `eta / rho^n`
constant in `n` to four digits, and a Kato-Temple bracket four orders tighter than the
linear bound at `n = 20`. The refusal case at `V = 0` behaves as section 8.9 requires, and
the isolation condition and the refusal criterion turn out to be the same inequality.

That spike has already paid for itself once: the arithmetic floor now carried by every
Phase 2 certificate came out of its last table.

**The core surface it couples through exists and is exported.** `set_provenance()`,
`strip_provenance()`, `provenance()` and the four generics `provenance_lift()`,
`provenance_refine()`, `provenance_original_residual()` and `provenance_summary()` are all
in `NAMESPACE` and all in `test-api-budget.R`'s `BUDGET`. Section 8.2's link between the
Hilbert operator and its image under a discretisation is what those seven names are for,
and building the layer costs core no export.

**The first unit is self-adjoint, so it needs nothing from Phase 3 or Phase 4.** It wants
`eigs()` and the certificate, both of which are in. It does not want Arnoldi, PRIMME,
RSpectra or an adapter. The v0.3 acceptance criterion at plan section 14 is reachable
against v0.1 core as it stands.

## 3. The envelope does not dispatch

The four generics dispatch on the envelope. `set_provenance()` stores the envelope as a
bare `list(provider =, payload =)` with no class, so the dispatch sees `"list"`:

| call | result |
|---|---|
| `provenance_lift(p, x)` on the envelope as stored | `no provenance_lift() method registered by provider 'linop.hilbert'` |
| same envelope, `structure(p, class = "hilbert_provenance")` | provider method reached |
| `class(payload)` after `set_provenance()` | `hilbert_payload`, preserved |
| `class(envelope)` after `set_provenance()` | `list` |

So a provider cannot reach its own method. Three of the four generics land on a `.default`
that errors naming the provider, and `provenance_summary()` lands on one that returns
`"provenance from '<provider>' (no summary method registered)"`. Every one of those
messages is accurate about the symptom and misdescribes the cause: the method may be
registered and is unreachable.

The class a provider puts on its *payload* does survive. Dispatching on it works:

```r
prov_lift_v2 <- function(p, x, ...) UseMethod("prov_lift_v2", p$payload)
#> provider method reached via payload class
```

Three routes, and the last is the recommendation:

- **Provider classes the envelope itself.** Needs `set_provenance()` to take and keep a
  class, so core grows an argument and the provider has to remember to use it.
- **Core derives the class from `provider`**, e.g. `c(provider, "linop_provenance")`. Uses
  a string core already holds, and makes a package name into a class name, which is a
  coupling section 5.11 would rather not have.
- **Dispatch on the payload.** `UseMethod("provenance_lift", p$payload)` in all four
  generics. Core still never looks inside the payload, it hands the class to `UseMethod`
  and the provider owns both. One token per generic, no signature change, no new export,
  no budget movement.

Whichever is chosen, `test-provenance.R` currently asserts round-trip and propagation and
does not assert that a registered method is reachable. That test is the gap that let this
through, and it is the one to add first.

## 4. What the first unit is, and is not

In scope, and it is exactly S0.6's class:

- `FiniteSection` on `ell^2(Z)`, self-adjoint Jacobi plus finite-support `V`.
- `discretise()` returning a `linop` carrying the provenance envelope.
- The three-part certificate: original-operator residual, truncation bound, isolation gap,
  with the arithmetic floor S0.6 forced.
- Refusal at `V = 0`, on `q - a > 0` rather than on a separate check.

Out of scope, and this is the part with the research cadence that the package split exists
for: the space hierarchy and its conversion graph (8.3), `ClosedOperator` domains and
adjoint domains (8.4, which already makes unknown first-class because the symbolic route
runs out past constant or polynomial coefficients), forms as peers (8.5), and `Galerkin`
and `Collocation` (8.6). Plan section 15's constraint applies to all of it: the layer must
never be described as certifying more than it certifies.

## 5. Still open

- **Repo layout.** Two packages is settled by install cost. One repo with two package
  directories against two repos is a tooling question, not a design one, and is undecided.
  Separate repos suit r-lib workflows, pkgdown and r-universe topic discovery; a monorepo
  costs `working-directory` overrides and buys a single tracker.
- **The name.** Plan open question 2 still has `linop.hilbert` as a proposal.
- **Whether the first unit needs any new core export.** Unknown until it is built. `BUDGET`
  is a test, and editing it as a side effect of adding a function is the thing that test
  exists to catch.
- **`Additional_repositories:`** if core ever gains a `Suggests` on the layer before the
  layer is on CRAN.

## 6. Corrections

| Where | Was | Is |
|---|---|---|
| Assessment given in session before reading S0.6 | The Hilbert layer is the weakest of the three splits and should not be a committed deliverable | The split reasoning is unchanged, but the first unit is the one part of section 8 with executed evidence behind it. Risk 5's control was run as a spike and returned feasible |
| Plan Phase 5 ordering | Hilbert unit follows `linop.primme` and the adapters | Neither is a prerequisite. The unit is self-adjoint and couples through provenance, which is in v0.1 |
| `R/provenance.R` | Four generics a provider registers methods on | As stored the envelope is unclassed, so the methods are unreachable and the `.default` messages name the wrong cause |
