# The first Hilbert unit against v0.1 core, and the certificate a provider cannot build

Run 2026-08-02. Script `dev_notes/spikes/hilbert-first-unit-probe.R`, output in
`results/hilbert-first-unit.txt`, `results/hilbert-first-unit.csv` and
`results/hilbert-first-unit-surface.csv`. Every number below is from that run.

The open question blocking the package decision was "does the first unit need any new core
export?", and `dev_notes/hilbert-first-and-the-envelope-that-does-not-dispatch.md` records
it as unknown until the unit is built. So the unit was built, as a provider would build it,
using only names a separate package could call.

## Verdict: one export question, and it is the certificate

The operator path needs nothing. The unit constructs, verifies, carries provenance and goes
through `eigs()` with **zero** names that are not already in `NAMESPACE`. What it cannot do
is produce a `linop_certificate`.

## 1. What runs, and against what

`finite_section(V, n)` is a matrix-free `linop` of dimension `2n+1` whose apply is three
shifted adds, carrying an envelope with a `finite_section` payload class. The four
provenance generics are registered against that class and all four route, including through
`t(H) %*% H` where the envelope is the propagated one. `verify(H)` returns `pass` on eleven
checks.

S0.6's table reproduced through `eigs()` rather than through dense `eigen()`, `V = 1
delta_0`, exact `lambda = sqrt(5)`:

| n | lambda_n | true error | eta | floor | error <= eta + floor | eta / rho^n |
|---|---|---|---|---|---|---|
| 5 | 2.233209543058 | 2.858e-03 | 5.407e-02 | 2.665e-15 | yes | 0.5997 |
| 10 | 2.236045382463 | 2.260e-05 | 4.754e-03 | 2.665e-15 | yes | 0.5847 |
| 20 | 2.236067976007 | 1.493e-09 | 3.864e-05 | 2.665e-15 | yes | 0.5845 |
| 40 | 2.236067977500 | 4.441e-16 | 2.554e-09 | 2.665e-15 | yes | 0.5844 |
| 80 | 2.236067977500 | 8.882e-16 | 7.734e-17 | 2.665e-15 | yes | 4.0499 |

`eta` is two entries of the computed eigenvector and nothing else, and the floor is the
provider's own `4 (2 + max|V|) eps`, closed form for this class with no norm estimate
anywhere. The n = 80 row is S0.6's failing row and it now holds, which is the floor doing
the work it was added for.

**The last column reads the floor from the other side.** `eta / rho^n` is constant to four
digits at 0.584 through n = 40 and then breaks to 4.05, because at n = 80 the analytic
`eta` is 1.1e-17 and the measured one is 7.7e-17: `eta` itself has reached the arithmetic
plateau. The truncation bound stops being a truncation measurement before the certificate
stops being issued, which is why the floor is a term and not a tolerance.

Kato-Temple against the band edge `a = 2`, known from the symbol rather than estimated:

| n | eta (linear) | KT half-width | contains sqrt(5) |
|---|---|---|---|
| 5 | 5.407e-02 | 1.254e-02 | yes |
| 10 | 4.754e-03 | 9.576e-05 | yes |
| 20 | 3.864e-05 | 6.325e-09 | yes |
| 40 | 2.554e-09 | 2.763e-17 | yes |

## 2. The certificate is the one thing a provider cannot build

`build_certificate()`, `cert_rows()`, `arithmetic_floor()` and `norm2()` are all internal.
A satellite package can compute every number in its own certificate -- it just did, above --
and has no way to put them in the object the rest of the package reports.

The Hilbert certificate is a **different shape**, the fourth: truncation bound, isolation
gap and arithmetic floor, with no residual row and no backward-error row, because a
finite-section eigenvalue is not the answer to an equation with a right-hand side. That is
the same argument that made the eigenpair certificate a third shape rather than a flag on
the second. But shape is the rows; the *object* is the row table, the three evidence fields,
the `overall` roll-up and the print method, and those are identical across all four.

So the fork is:

- **Export the builder.** `cert_rows()` and `build_certificate()`, two names, and the
  satellite's certificate prints and rolls up like every other one. `arithmetic_floor()` is
  a third if the floor is to be spelled the same way rather than re-derived; `norm2()` is
  not needed here, since `||H|| <= 2 + max|V|` is closed form for this class.
- **Reimplement it in the satellite.** No budget movement, and two copies of the roll-up
  rule, the evidence fields and the print method in two packages, drifting apart at the
  first change to either.

The second is the facade-over-duplicate outcome, so the first is the recommendation, but it
is a `BUDGET` edit and therefore the user's call rather than a side effect of adding a
function.

Nothing else was reached and not exported. `cap(A, name)$value` is the public read; the
internal `cap_value()` is a convenience with an unknown-default and the unit does not need
it.

## 3. A matrix-free provider cannot say `construction`

The finite section is hermitian because it is tridiagonal with unit off-diagonals and a real
diagonal. The provider knows this the way `linop_eye()` knows it about the identity. It
cannot say so:

| leaf | value | source | guarantee |
|---|---|---|---|
| matrix-free, `properties = c("hermitian")` | TRUE | `user_declaration` | `identity` |
| the same operator as a dense leaf, no properties | TRUE | `computation` | `identity` |

`normalise_properties()` stamps `ev_declared()` unconditionally, and `ev_construction()` is
internal with no `properties=` form reaching it. So the route that materialises the operator
gets the better evidence and the matrix-free route -- the one the layer exists for -- gets
the weakest source in the package.

This is not cosmetic here. `eigs()` records the hermitian capability's evidence under
`depends_on` on the forward-error row, precisely so that a Weyl bound resting on a bare
declaration fails a requirement the declaration would have failed directly. For the finite
section that bound is provable and reports as declared, which is the laundering case running
the wrong way: a true `construction` argument downgraded rather than a weak one promoted.

Three routes, none costed yet: a `properties=` form accepting `capability()` objects
(no new export, `capability()` and `evidence()` are both already public), an exported
`ev_construction()`, or an `adapter_contract` convention for adapter-authored leaves. The
first looks smallest and reuses two exports that exist for exactly this vocabulary.

## 4. The refusal, and what the eigensolver does not tell you

`V = 0`, where `H_free` has no eigenvalues at all:

| n | q (top Ritz) | q - a | certificate issued | eigs nconv | eigs certificate |
|---|---|---|---|---|---|
| 10 | 1.979642883762 | -2.036e-02 | no | 1 of 1 | qualified |
| 40 | 1.998532362102 | -1.468e-03 | no | 1 of 1 | qualified |
| 100 | 1.999602918148 | -3.971e-04 | no | 0 of 1 | fail |

**The eigensolver is a second signal and not a substitute for the first.** At n = 10 and 40
it converges cleanly and returns a value that is still inside the band, so a layer trusting
`nconv` would issue a certificate for a discretisation of continuous spectrum. At n = 100 it
stalls, and the stall is about clustering near the band edge rather than about the absence
of an eigenvalue -- the same stall a genuinely hard isolated eigenvalue would produce. Only
`q - a > 0` separates "no eigenvalue here" from "a hard eigenvalue here", which is section
8.9 and section 8.10 turning out to be one inequality for the second time.

## 5. Corrections

| Where | Was | Is |
|---|---|---|
| `hilbert-first-and-the-envelope...`, section 5 | whether the first unit needs any new core export is unknown until it is built | Built. The operator path needs none; the certificate object needs `cert_rows()` and `build_certificate()`, or a second copy of them in the satellite |
| Implicit in the provenance surface | a provider states its capabilities through `properties=` | It states values through `properties=` and cannot state evidence. A construction argument reports as `user_declaration`, and `eigs()` carries that into the forward-error row's `depends_on` |
| S0.6 section (d) | at `V = 0` the finite-section eigenvalues lie inside the band, so `eigs()` must refuse | True, and `eigs()` does not refuse -- it converges and returns one. The refusal is the layer's, on `q - a > 0` |
