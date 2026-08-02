# One package, and the abstractions that were only holding a boundary

Decision, 2026-08-03. This supersedes plan section 3's package-split rule, the
three-package model, and `dev_notes/hilbert-first-and-the-envelope-that-does-not-dispatch.md`
section 5. Nothing here is measured; it is a design and a work order.

## The premise

`linop` owns the whole mathematical system and the complete user-facing API. A user
installs and learns one package. RSpectra and PRIMME remain external numerical
libraries, and their integration lives inside `linop` behind `method =`. Whether a
library is an `Imports`, a `Suggests`, or an optional compiled component is an
installation detail, and installation details do not fragment an API.

The target:

```r
library(linop)
A <- linop_jacobi(diagonal = v_at_origin)   # an operator on l^2(Z)
F <- finite_section(A, n = 20)              # a linop, and a node that knows its parent
eigs(F)                                     # native
eigs(F, method = "rspectra")                # same verb, external engine
eigs(F, method = "primme")                  # same verb, compiled engine
solve(A2, b, method = "minres")
certificate(F)
```

## The reassessment

Each mechanism below existed to hold a package boundary. With no boundary, each has
to earn its place again on its own.

### Opaque provenance envelopes and provider-owned payloads: **delete**

Section 5.11's rationale was that hard-coding a schema naming a discretisation would
"leak Layer 3 concepts into core". Layer 3 *is* core now. There is no leak because
there is no boundary, and opacity buys exactly one thing -- core not needing to know
the payload's structure -- which is worth nothing when core owns the structure.

What replaces it is stronger, not weaker. **A finite section is a node.** `linop`
already has a typed node registry with propagation, printing, materialisation and
composition; a truncation is a shape-changing node over one child, in the same family
as the deferred `stacks` and `blockdiag`. The child *is* the original operator, so
"this finite thing approximates that infinite thing" is structural rather than
annotated, and `explain()` and `print()` show it with no extra mechanism.

This removes seven exports (`provenance`, `set_provenance`, `strip_provenance`, and
the four generics), one payload class, the `UseMethod(generic, p$payload)` dispatch
rule, and the bug that rule was introduced to fix. `db7e289` fixed a real defect in a
mechanism that should not exist.

Touches `core.R`, `nodes.R`, `print.R`, `provenance.R`.

### Cross-package generics: **delete**

They are the four above. Nothing else crossed.

### Exported certificate builders: **back to internal**

`build_certificate()` and `cert_rows()` were internal for all of Phase 2 and became
public yesterday for exactly one consumer, in another package, written the same day.
The validation added with them is kept -- a closed status vocabulary and a checked
evidence vocabulary are worth having on an internal constructor too -- but the names
come out of the public surface.

What is public instead is an **accessor**: `certificate(x)`. Users read certificates;
they do not build them.

### Node registration machinery: **internal**

`linop_register_node()` and `linop_nodes()` are how a package outside this one would
add a composition type. There is no package outside this one. The extension axis that
survives for a *user* is `linop.<class>()`, plain S3, which needs no registry.

The deferred node types are linop's own work and register internally.

### Finite sections treated as a provider: **gone with the envelope**

They become `R/section-*.R`, a node type and a certificate shape, indistinguishable
in kind from the product node or the eigenpair certificate.

### What survives, and why

- **Evidence and capabilities.** Never about packaging. `evidence()`, `capability()`,
  `requirement()`, `evidence_satisfies()`, `cap()` stay.
- **The certificate object and its shapes.** The fifth shape is the finite section's.
- **`verify()`.** It checks an operator against the contract, which a user writing a
  callback needs whoever wrote the callback.
- **`method =`.** Already the algorithm selector; now the engine selector too.

## The consequence that is not packaging

**Arnoldi has to be native.** `eigs()` currently refuses a non-hermitian operator and
the plan defers Arnoldi until a backend arrives -- a capability gap in `linop` filled
by a satellite. Under this premise that is backwards: a package owning eigensolvers
cannot be missing the non-hermitian one, and Arnoldi is pure R with no dependency.
The generalized problem `A x = lambda B x`, currently refused by name, is the same
argument.

This largely deletes the case for RSpectra as a source of capability.
`dev_notes/rspectra-and-the-delegation-that-is-not-a-superset.md` measured that it
refuses a block apply, a complex dense matrix and `sigma` on a function, and that it
silently coerces a complex operator and returns `Re(A)`'s spectrum, wrong by 6.946
against the truth. With Arnoldi native, `method = "rspectra"` is an optional
accelerator for a narrow case, and it must refuse on `linop`'s own dtype before
delegating or it inherits that coercion.

PRIMME stays external for one reason: it is compiled, and `linop` installs with no
compilation today. `Suggests` plus a gated call, not a second package.

## The public surface this produces

| | |
|---|---|
| removed | `provenance`, `set_provenance`, `strip_provenance`, `provenance_lift`, `provenance_refine`, `provenance_original_residual`, `provenance_summary`, `build_certificate`, `cert_rows`, `linop_register_node`, `linop_nodes` |
| added | `finite_section`, `certificate`, `decay_rate`, `linop_jacobi` |

Thirty-two becomes twenty-five, covering strictly more mathematics. The budget cap
comes *down* to match, which is the first time it has moved in that direction.

`linop_jacobi()` is the one place this deviates from the sketch, which wrote
`A <- linop(...)`. `linop()` is a generic dispatching on its first argument, and an
operator given by coefficient sequences has no `x` to dispatch on; it belongs with
`linop_eye()` and `linop_scaling()` in the constructor family.

## Work order

Steps 1 to 3 are done. Step 3's own findings, including two corrections to this
list, are in
`dev_notes/the-section-and-the-certificate-about-an-operator-that-is-not-there.md`.

1. Delete provenance. Un-export the certificate builders and the node registry. Move
   the budget down. This is subtraction and it lands first, so nothing later is built
   against machinery that is going away.
2. `dim` admits `Inf`. Audit every consumer: `verify()` cannot probe an infinite
   operator, `as.matrix()` must refuse it, the algebra must still compose it. This is
   the real engineering and the rest depends on it.
3. `linop_jacobi()` and the `section` node, with the five-row certificate.
4. `eigs()` generic, `certificate()` accessor.
5. Arnoldi, closing the non-hermitian gap natively.
6. `method = "rspectra"` and `method = "primme"` as gated optional engines.

## Corrections

| Where | Was | Is |
|---|---|---|
| Plan section 3 | anything at research pace or costing an install second goes in a separate package | Capability never leaves `linop`. Only a compiled engine is external, and only because it is compiled |
| This session, twice | `linop.hilbert` as a satellite, then merged while keeping the provider architecture internally | Neither. The mechanisms that held the boundary go with it |
| Section 5.11 | core carries provenance and never inspects it | Deleted. A truncation is a node with a child, which is structural rather than annotated |
| `db7e289` | the provenance generics must dispatch on the payload | A real fix to a mechanism that should not exist |
