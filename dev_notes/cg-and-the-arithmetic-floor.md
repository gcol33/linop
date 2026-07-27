# CG, the norm estimate, and the arithmetic floor

Written 2026-07-27, entering Phase 2. Three pieces landed together because none of them is
testable alone: the certificate needs `||A||`, the floor needs the certificate, and a
certificate with no solver behind it cannot meet its own coverage gate.

All numbers below are reproduced by `dev_notes/spikes/cg-floor.R`.

## 1. The floor, in the linear-solve setting

S0.6 measured the correction and stated it for an eigenvalue bound:

```
dist(lambda_n, sigma(H))  <=  eta  +  c * ||H|| * eps
```

For a linear solve the same term has a sharper origin. Higham bounds the *computed*
residual: `fl(b - A x)` differs from `b - A x` by at most a small multiple of
`eps * (||A|| ||x|| + ||b||)`, so no `x`, however exactly it solves the system, can be
shown to have a residual below that level. The Rigal-Gaches identity gives the normwise
relative backward error of `x` as

```
omega(x) = ||b - A x|| / (||A|| ||x|| + ||b||)
```

with the *same* denominator, so once the backward error is expressed relatively the floor
is a plain `c * eps`. The residual line and the backward-error line share one quantity
rather than needing two separate roundoff models.

`c` counts terms in an inner product a matrix-free operator never exposes, so it is not
knowable here. It is an argument, `floor_const`, not a constant. The default is 4, from
S0.6's own measurement: a plateau near 2.2e-15 against `||H||` between 4 and 6 puts `c`
between 2 and 4.

**Measured, `tol = 0` on a 30 x 30 operator with spectrum in [1, 4].** Asking for an exact
solve is the cleanest way to reach the plateau, because nothing in floating point delivers
one.

| quantity | value |
|---|---|
| relative residual achieved | 1.789e-16 |
| relative floor at `c = 4` | 2.529e-15 |
| certificate at `c = 4` | residual `qualified`, overall `qualified` |
| certificate at `c = 0` | residual `fail`, overall `fail` |
| iterates in the two runs | bitwise identical |

The last row is the point. The floor changes nothing about the answer and everything about
whether the package calls a converged answer a failure.

**A line that meets its tolerance only through the floor is reported as `qualified` with an
`estimate` guarantee, never as a clean `identity`.** The floor is a model of roundoff, not
a theorem, and section 6's two-dimensional certificate already has the field to say so.

## 2. `||A||` and why a lower bound is the safe error

The floor scales linearly in `||A||`, so the norm has to come from somewhere. Section 5.3
lists `norm_bound(type)` with evidence among the optional operator capabilities, but that
protocol is Phase 3 and the floor is Phase 2, so `norm2()` computes it now and carries
evidence for the same reason a capability does.

| operator | route | value | exact |
|---|---|---|---|
| `linop_eye(9)` | structure | 1.000000 | yes |
| `linop_scaling(c(2, -5, 3))` | structure | 5.000000 | yes |
| `3 * linop_scaling(c(2, -5, 3))` | structure | 15.000000 | yes |
| dense 20 x 20 | svd | 7.853595 | yes |
| Laplacian, n = 200 | power | 3.854263 | no (truth 3.999756) |

Every route returns a **lower** bound: `||A v||` for a unit `v` never exceeds `||A||_2`,
and the structural rules are equalities. That direction is the safe one. A smaller `||A||`
shrinks the Rigal-Gaches denominator, so the reported backward error can only grow, and the
certificate errs towards overstating the error rather than understating it. A test asserts
this over 20 seeds by forcing the iteration route and comparing against the exact norm.

Two consequences worth recording:

- **The guarantee is `estimate`, not `deterministic_bound`,** even though the power
  iteration is deterministic and its output is genuinely a bound. "Bound" on a norm reads
  as an upper bound, and this is the other side; labelling it `deterministic_bound` would
  be accurate and would still mislead.
- **A structural rule over an estimate stays an estimate.** `norm2(4 * A)` multiplies an
  exact rule by an estimated child, and the composite records
  `construction <- [computation/estimate]`, so `evidence_satisfies(guarantees = "identity")`
  returns `FALSE` at the top. This is the laundering case of section 5.3 appearing
  somewhere other than capabilities, and the same recursion catches it.

CG also hands the certificate `max ||A p|| / ||p||` observed over its own iteration, at no
extra apply. It is another lower bound, so the larger of the two is used whenever the norm
is not exact.

## 3. CG

**Recovery against closed-form truth** (section 10's 1-D Dirichlet Laplacian, eigenvalues
`4 sin^2(k pi / (2(n+1)))`), `tol = 1e-12`:

| n | kappa | iterations | relative residual | max abs error |
|---|---|---|---|---|
| 20 | 178 | 20 | 9.50e-16 | 3.33e-15 |
| 60 | 1.51e3 | 60 | 1.54e-15 | 1.11e-14 |
| 200 | 1.64e4 | 200 | 5.85e-15 | 1.32e-13 |

Iterations equal `n` exactly in all three, which is finite termination holding in floating
point at these condition numbers.

**The sharper Krylov statement**, which a tolerance test cannot make: with `d` distinct
eigenvalues the CG polynomial of degree `d` annihilates the error, so the run must finish
in `d` steps and not merely fast.

| distinct eigenvalues | iterations |
|---|---|
| 3 | 3 |
| 5 | 5 |
| 10 | 10 |

**Reference agreement.** A textbook recurrence on the dense matrix reaches the same
iteration count and a **bitwise identical** iterate on all five seeds tried, at n = 50 with
a spectrum spanning [0.2, 20]:

| seed | this implementation | reference | max gap |
|---|---|---|---|
| 1 | 66 | 66 | 0 |
| 2 | 68 | 68 | 0 |
| 3 | 68 | 68 | 0 |
| 4 | 68 | 68 | 0 |
| 5 | 68 | 68 | 0 |

This is what licenses the two structural departures below.

### Several right-hand sides run in lockstep

Each column's recurrence is independent of the others, so running them together with
per-column scalars gives exactly the iterates of per-column CG while spending one block
apply per step instead of `k`. Converged columns leave the active set, so the block narrows
as the solve proceeds. Measured on 8 right-hand sides, n = 40:

- **8 of 8 columns bitwise identical** to their own separate solves
- **55 block applies** against 440 per-column applies

This is not block CG, which couples the columns and is a different method with different
iterates. It is the two-tier apply of section 5.2 paying for itself.

### The outer loop measures, the inner loop trusts

The recurrence residual drifts, and it drifts worst exactly where the answer is most
converged. The inner loop iterates on the recurrence; the outer loop recomputes `b - A x`
and decides. A restart that does not reduce the true residual stops rather than looping.
The certificate always reports the recomputed residual, at the cost of one apply, which is
the whole difference between a certificate and a progress report.

**What this buys is honesty, not accuracy.** Measured at `kappa = 1e6`, n = 50,
`tol = 1e-12`, against the same textbook recurrence trusting itself:

| | iterations | residual it reports | true relative residual | verdict it issues |
|---|---|---|---|---|
| textbook, self-trusting | 399 | 6.912e-13 | 3.297e-11 | tolerance met, and it was not |
| this implementation | 500 (1 restart) | 3.950e-11 | 3.950e-11 | the number it measured, qualified against the floor |

Neither run reached `tol = 1e-12` in the true residual. The restart spent 100 further
iterations and did not improve the answer; across six seeds it landed better than the
reference on three and worse on three. What changed is that a solver trusting its own
recurrence reports meeting a tolerance it missed by a factor of 50, and this one reports the
number it recomputed wherever that lands. The verdict it then issues is the next section,
and it is not "no".

### One number could not have told this story

Same system. The certificate splits into two lines that disagree, and both are right:

```
check                          status       source           guarantee             conf
------------------------------------------------------------------------------------------
arithmetic floor               pass         computation      identity              1
residual                       qualified    computation      estimate              -
backward error                 pass         computation      identity              1
convergence                    pass         computation      identity              1
forward error                  not_checked  computation      identity              -
------------------------------------------------------------------------------------------
overall                        qualified   no deterministic bound on: residual, forward error
```

with `||A|| = 1e6` exactly by the svd route, a relative floor of 2.019e-10, a residual of
3.950e-11 and a backward error of 1.738e-16.

`||A|| ||x||` dwarfs `||b||` here, so the residual sits above the requested tolerance and
*below* the arithmetic floor, which is `qualified` with an `estimate` guarantee. The
backward error is at machine epsilon: `x` is the exact solution of a system
indistinguishable from the one asked about, a clean `pass` with an `identity` guarantee.
Collapsing these into one status would have had to pick a story. A test pins both.

This is also what `converged = TRUE` means on this run. The tolerance was met once the
arithmetic floor is allowed for, which is the strongest statement floating point supports
here; reading `qualified` as failure would reinstate exactly the S0.6 behaviour the floor
exists to remove.

### Sides

Section 4.3 leaves the CG row unrestricted, and that is a statement rather than a gap: for
hermitian positive definite `M` the left, split and right preconditioned forms generate the
same iterates (Saad, section 9.2). One implementation therefore serves all three, and since
the reported residual is always the true `b - A x` rather than a preconditioned one,
nothing downstream depends on which side was declared. The test asserts the contract, that
none of the three is refused; the mathematical coincidence is Saad's and is cited, not
re-derived by a test.

### Refusals, each naming what it needed

- `positive_definite` unknown: refused, with `Unknown is not false` and the declaration
  that would settle it. `hermitian` likewise.
- A rectangular operator: refused, naming least squares as a different request.
- A false `positive_definite` declaration: `linop_scaling(c(-1, 1, 1, 1))` declared
  positive definite gives `p^H A p = -1` on the first step with `b = e_1`, so the refusal is
  structural rather than a matter of which seed was drawn. Near convergence `p` can
  underflow, so a column whose residual has already fallen by `sqrt(eps)` is frozen instead;
  anywhere else the declaration is wrong and CG says so.
- A false declaration on the preconditioner, by the same route through `<r, M^-1 r>`.
- A preconditioner that returns the wrong shape.

### Non-convergence is returned, not thrown

Jacobi on a badly scaled system, `kappa = 2.91e8`, `tol = 1e-12`:

| run | iterations | relative residual | certificate |
|---|---|---|---|
| unpreconditioned | 600 (budget) | 3.017e-02 | `fail` |
| Jacobi | 40 | 4.419e-13 | `qualified` |

Unpreconditioned CG loses orthogonality two digits short and stalls. The result comes back
with a failing certificate naming the exhausted budget rather than an error, which is what
"the result records what happened" in section 1.1 has to mean when what happened was not
success.

## 4. Evidence: naming the method is not the same as `auto`

Plan 7.1 says CG requires `positive_definite` at declared minimum evidence, and 7.1 also
says `method = "auto"` never picks a method whose evidence is below its declared minimum.
Those are two different gates and only the second is a filter:

- **`method = "cg"`** requires the capability to be `TRUE` with any evidence. Naming the
  method is the caller asserting their own declaration, and refusing it would make
  `linop(M, properties = c(positive_definite = TRUE))` unsolvable by the method it was
  declared for.
- **`method = "auto"`** applies `CG_PD_REQUIREMENT`, which admits `construction`,
  `adapter_contract`, `theorem` and `computation` at an `identity` guarantee and confidence
  1.

`computation` is in that list because an exact check on data the operator already holds is
an identity, not a probe: `linop_diag()` proves definiteness from the signs of `d`. Probes
are excluded by the **guarantee** field, which `ev_probe()` never satisfies, rather than by
the source list. The test-local `cg_requirement()` in `test-evidence.R` omits `computation`;
it is an illustration of the machinery and its assertions hold either way.

## 5. Corrections to the plan

| Plan location | Correction |
|---|---|
| 6, certificate table | an `arithmetic floor` row, and the roundoff term on `residual` and `backward error`. Without it the results that certify as `fail` are exactly the most converged ones |
| 6 | a solve certificate is `qualified`, not `pass`, whenever `forward error` is unchecked, which matrix-free it always is. The plan's own example table already has this shape; it was not stated as the expected outcome |
| 7.1 | CG's evidence requirement is a filter for `auto`, not a gate on the named method |
| 5.3 | `norm_bound(type)` is listed as an optional capability for Phase 3, but the floor needs `||A||` in Phase 2. It ships now as an internal estimate carrying its own evidence, and the capability protocol can adopt it later |

## 6. Open

- **`forward error` stays `not_checked`.** It needs `||A^-1||`, which no residual implies.
  A Hager-Higham estimate of `||A^-1||_1` is reachable matrix-free, since it needs only
  solves against `A` and `A^H`, and CG supplies both. That is a real follow-up rather than
  a permanent gap, and until it exists section 6.1's restraint applies: the honest answer is
  that it was not checked.
- **The six unrestricted side rows** in `PRECOND_REQUIREMENTS` are still the conservative
  reading. Narrowing them needs Saad chapter 9 and should not be guessed at.
- **`solve()` is not exported.** `method = "auto"` cannot mean anything with one method to
  choose from, and the API budget asserts Phase 2 names absent until they are deliberately
  promoted. It lands when there is a choice to make.
