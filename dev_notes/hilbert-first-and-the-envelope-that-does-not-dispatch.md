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

**Fixed.** All four generics now dispatch on `p$payload`, the third route below, and the
test was rebuilt to obtain its envelope from `set_provenance()`. The counterfactual was
run rather than assumed: against the previous generics the rebuilt test fails 2 of its
assertions with `no provenance_lift() method registered by provider 'demo'` while the
method is registered, which is the misdescription this section is about. `NAMESPACE` is
byte-identical afterwards. The rest of the section is what was measured before the fix.

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

Whichever is chosen, the test that should have caught this is already there and does not.
`test-provenance.R:46` is named "a provider can register methods and core will route to
them", and it builds its envelope by hand:

```r
p <- structure(list(provider = "demo", payload = list(n = 256)), class = "demo_prov")
```

That `structure()` is the one step `set_provenance()` does not do. So the test asserts the
generics dispatch, which they do, and never asserts that what core hands back is
dispatchable. The half the title names is the half not exercised. Constructing the fixture
the way a caller would would have failed it. That is the edit to make first, and it is an
edit rather than an addition.

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

- ~~**Repo layout.**~~ **Two repos.** The case for a monorepo is atomic cross-package
  commits, and it does not apply: the satellite depends on a *released* `linop`, not on
  `linop@main`, so the coupling is versioned and asynchronous and a monorepo would only
  make it look otherwise. Against it: r-universe discovery is driven by a per-repository
  GitHub topic, every install instruction grows a `subdir`, every CI step needs a
  `working-directory` override plus path filters, and two pkgdown sites need two deploy
  workflows. The single tracker is a real win and does not cover four frictions.
  `linop.primme` makes it three packages, and three uniform repos beat a monorepo of two
  plus an outlier.
- ~~**The name.**~~ **`linop.hilbert`**, and the dot is forced rather than chosen. R package
  names take only letters, digits and dots, so a family is marked by a dot or by
  concatenation and nothing else; concatenation breaks on this stem (`linophilbert` reads
  as "lino-philbert", `linopprimme` has a double seam). **One constraint comes with it:
  the satellite must never name a class `hilbert`.** `linop()` is an S3 generic whose
  adapter convention is `linop.<class>()`, so a class of that name would produce a method
  whose name is the package's. Nothing breaks technically, since the two live in different
  namespaces, but it is the kind of collision that surfaces years later as a confusing
  `R CMD check` message.
- **Sequencing.** `linop` goes to CRAN first. While it is not there, the satellite's
  `Imports: linop` needs `Additional_repositories: https://gcol33.r-universe.dev` with a
  live `PACKAGES` file, which is the failure mode taxify 0.2.5 was rejected for. With core
  on CRAN the field is not needed at all, and it must stay that way: the dependency runs
  one direction only and core must never gain a `Suggests` on the layer.
- ~~**Whether the first unit needs any new core export.**~~ Answered by building it:
  `dev_notes/hilbert-first-unit-and-the-certificate-a-provider-cannot-build.md`. The
  operator path needs none. The certificate object needed `cert_rows()` and
  `build_certificate()`, and they are exported, deliberately and with the budget edited to
  match. The surface now sits at exactly the 32 `test-api-budget.R` allows.
- **`Additional_repositories:`** if core ever gains a `Suggests` on the layer before the
  layer is on CRAN.

## 6. Corrections

| Where | Was | Is |
|---|---|---|
| Assessment given in session before reading S0.6 | The Hilbert layer is the weakest of the three splits and should not be a committed deliverable | The split reasoning is unchanged, but the first unit is the one part of section 8 with executed evidence behind it. Risk 5's control was run as a spike and returned feasible |
| Plan Phase 5 ordering | Hilbert unit follows `linop.primme` and the adapters | Neither is a prerequisite. The unit is self-adjoint and couples through provenance, which is in v0.1 |
| `R/provenance.R` | Four generics a provider registers methods on | As stored the envelope is unclassed, so the methods are unreachable and the `.default` messages name the wrong cause |
| This note, first draft | `test-provenance.R` does not assert a registered method is reachable | It does, at line 46, and passes: it hand-classes the envelope instead of obtaining one from `set_provenance()`. The assertion is present and the caller's path is not |
