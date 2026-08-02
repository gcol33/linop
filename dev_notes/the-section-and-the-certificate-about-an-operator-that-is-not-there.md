# The finite section, and a certificate about an operator that is not there

Step 3 of `dev_notes/one-package-and-the-abstractions-that-were-holding-a-boundary.md`,
2026-08-03. `linop_jacobi()`, `finite_section()` and `decay_rate()` are in, with the
certificate shape S0.6's mathematics was derived for. Every number below is from
`dev_notes/spikes/section-ncv-probe.R` or from `tests/testthat/test-section.R`, both
of which run from a clean load.

## What landed

`linop_jacobi(diagonal, offdiagonal, from, diagonal_limit, offdiagonal_limit)` is a
self-adjoint operator on `l^2(Z)`, tridiagonal in the standard basis, with both
coefficient sequences real and eventually constant. `finite_section(A, n)` is
`P_n A P_n` as a **node holding the operator it truncates**, so the printed tree
shows an `Inf x Inf` operator underneath a finite one and `explain()` walks to it
with no extra mechanism. That is what step 1 promised when it deleted the
provenance envelope, and it costs nothing: `node_children()` gained one word.

`eigs()` on a section returns the fifth certificate shape. Which builder runs is the
node type's to say -- `linop_register_node()` gained a `certify` field and `eigs()`
gained one line -- rather than `eigs()` testing for a node name.

The limits are arguments rather than fixed at `0` and `1`, which costs two arguments
and generalises the class from a discrete Schrodinger operator to any Jacobi operator
with an eventually constant coefficient pair. It is one affine change of variable:
with `mu = (lambda - a_inf)/b_inf` every closed form below is the one S0.6 derived.
`test-jacobi.R` asserts the recovery on `H = c I + s H_0`.

## The identity the whole certificate rests on

Let `u` solve the finite problem on `|j| <= n` and let `u~` be `u` extended by zero.
Then `(H - q) u~` agrees with `(H_n - q) u` at every `|j| <= n`, because the couplings
the truncation removed multiply entries that are zero, and it vanishes at every
`|j| > n + 1`. It is nonzero at exactly two sites, `b_n u_n` and `b_{-(n+1)} u_{-n}`,
and both coefficients are `b_inf` because `n` exceeds the window. So

```
||(H - q) u~||^2 = ||(H_n - q) u||^2 + b_inf^2 (u_n^2 + u_{-n}^2)
```

exactly. `test-section.R` checks it by embedding the computed pair in a section four
times wider and measuring there; it agrees to 1e-12.

**S0.6 has a term this does not.** Its table bounds the error by `eta` alone, because
it computed the finite eigenpair with dense `eigen()`, where the residual is at
rounding level and invisible. Run through `eigs()` at a tolerance the residual is not
negligible, and at `v = 1, n = 80` it is what dominates: residual 1.6e-10 against a
truncation term of 1.8e-11. So the certificate reports **two terms, not one**, and
they are separately actionable -- the first answers to `tol` and `ncv`, the second to
`n`. Folding them into one number would report the size of the error without saying
which knob moves it.

## Eight rows, where the plan sketched five

| row | what it is |
|---|---|
| arithmetic floor | `c eps ||H||`, with `||H||` the largest row sum, closed form |
| finite residual | `||H_n u - q u|| / ||u||`, what the eigensolve left behind |
| orthogonality | shared with the eigenpair certificate |
| truncation bound | `|b_inf| sqrt(u_{-n}^2 + u_n^2) / ||u||`, what the truncation costs |
| isolation | whether the value has cleared `sigma_ess`, known exactly by Weyl |
| enclosure | Temple's bracket against the band edge |
| convergence | shared with the eigenpair certificate |
| forward error | `dist(q, sigma(H)) <= sqrt(residual^2 + truncation^2) + floor` |

The sketch's five was written before the residual term was noticed. `orthogonality`
and `convergence` are the eigenpair certificate's own rows and are now shared code
(`cert_add_orthogonality()`, `cert_add_convergence()`, `cert_level()` in
`certificate.R`) rather than a second spelling of each.

## Three findings

### The forward-error bound is the first in the package that rests on nothing

`eigs()` on an ordinary operator records the hermitian capability under its Weyl
bound, so a bound resting on a bare `user_declaration` fails any requirement the
declaration would have failed directly. Here the class is self-adjoint **by
construction** and the residual identity is arithmetic, so the row's evidence has an
empty `depends_on` and satisfies `requirement(guarantees = "deterministic_bound")`
outright. It is also the only bound in the package about an object that was never in
memory.

That is the resolution of the hilbert spike's section 3, which found a matrix-free
provider forced to report `user_declaration` for a property it knew by construction,
and called it "the laundering case running the wrong way". Inside one package there
is a constructor to say `construction`, and no new export was needed for it.

### Temple's bracket covers one value per side, and it is not a deterministic bound

The bracket needs no spectrum of `H` between the band edge and the eigenvalue. That
is available for exactly one value on each side of the band -- the one adjacent to it
-- because for any other the near end of the interval is a neighbouring eigenvalue no
residual locates. Even for that one the hypothesis is not provable from a residual.

So the row carries `source = "theorem"` and `guarantee = "estimate"`: the theorem is
the source, and what it yields with an unverified hypothesis is not a deterministic
bound. `without_deterministic_bound` names it, which is the honest printed answer,
and the three evidence fields being independent is what lets one row say both. It is
worth 6.3e-9 against the linear bound's 3.9e-5 at `n = 20`, so it earns its place; it
does not earn the stronger guarantee.

The far end is free and needs no hypothesis at all: the Rayleigh quotient of `u~`
against `H` **is** the Rayleigh quotient of `u` against `H_n`, exactly, so the
variational principle puts the extreme eigenvalue beyond `q`. In exact arithmetic. A
test caught it failing at the last bit, so the reported end carries the floor -- the
S0.6 correction reaching a third row.

### The default subspace was wrong for a section, and the failure inverts with n

Confirmed and quantified. At `v = 0.3`, whose eigenvalue sits 2.2e-2 outside the
band, the ordinary default of `ncv = 21`:

| n | err at ncv 21 | err at ncv 40 |
|---|---|---|
| 20 | 1.707e-04 | 1.707e-04 |
| 40 | 4.241e-07 | 4.241e-07 |
| 80 | **9.803e-03** | 2.724e-12 |
| 140 | 4.441e-16 | 0 |

A **wider section giving the worse answer**, by four orders. The cause is that
enlarging `n` densifies the discretised continuum the wanted value has to be told
apart from rather than adding isolated eigenvalues, so a subspace adequate at one
size is inadequate at the next. More budget does not recover it, which is what
separates this from an ordinary shortage of iterations.

`ncv = 40` was then measured over `n = 80` to `800` and two potentials, and **it does
not have to grow with n**: the failure is a threshold rather than a trend. So
`default_ncv()` applies a floor rather than a formula, and a request for many pairs
still widens the subspace by the ordinary rule. The floor was measured at `k = 1`; at
`k = 3` on a well carrying three eigenvalues above the band, `ncv = 21` already
converges all three, because those values are far from the edge. The binding case is
proximity to the band, not the number of pairs.

**The narrow run is wrong and says so.** Residual 4.5e-2, `nconv` 0, certificate
`fail`. This was a bad default, never a silent wrong answer, and the certificate is
why. The archived note recorded this as "a defect in this package, found from outside
it, still open"; it is closed.

## What the truncation term stops measuring, and when

`eta / rho^n` is constant to four digits while the tail is resolved -- 0.2719 at
`v = 0.5`, 0.5845 at `v = 1`, 0.9852 at `v = 2`, 1.2629 at `v = 4`, reproducing S0.6
exactly through `eigs()` rather than dense `eigen()`. Past that it breaks, and the
level it breaks at is **the accuracy of the computed eigenvector, not the arithmetic
floor**: at `v = 1, n = 80` the analytic value is 1.1e-17, the arithmetic floor is
2.7e-15, and the measured term is 8.2e-13.

The bound holds throughout, because the identity above is about the vector actually
computed rather than about the exact eigenvector. Only the *interpretation* of the
term as a truncation measurement degrades. The row says so where the term has reached
the arithmetic floor, which under-fires: it is silent in the `v = 1, n = 80` case
above, where the term has plateaued four orders higher. A test that fired on the
eigenvector's own noise level would need a separation estimate, which is the same
thing `target identity` refuses to guess at elsewhere, so nothing was added.

## Corrections

| Where | Was | Is |
|---|---|---|
| One-package note, work order step 3 | a five-row certificate | Eight. The finite residual is a term the sketch omitted, because S0.6 measured it with dense `eigen()` where it is invisible |
| S0.6, the truncation bound | `dist <= eta + c eps ||H||` | `dist <= sqrt(res^2 + eta^2) + c eps ||H||`. The residual term dominates whenever the finite eigensolve runs to a tolerance rather than to machine precision |
| S0.6 section (c) | the Kato-Temple bracket against the known band edge | Correct, and it applies to one value per side. For any other the near end is a neighbouring eigenvalue, and the hypothesis is unverifiable in every case, so the guarantee is `estimate` |
| The hilbert spike, section 3 | a matrix-free provider cannot say `construction` | It can, once the constructor is in the package. No new export, and the forward-error bound is unconditional because of it |
| The archived first-unit note | `eigs()`'s default `ncv` at `k = 1` is a defect in this package, still open | Closed. A section gets a floor of 40, measured over `n = 80` to `800` |
