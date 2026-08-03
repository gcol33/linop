# The width that is a division, and the certificate that is read

Step 4 of `dev_notes/one-package-and-the-abstractions-that-were-holding-a-boundary.md`.
`eigs()` on the operator itself, choosing the truncation, and `certificate()`,
the last name under the 25 cap. Everything below is measured on
R 4.6.0, Windows 11, i9-14900K, unless it says otherwise.

Ground truth throughout is S0.6's single-site model: the potential `v` at the
origin of the free Jacobi operator, whose eigenvalue is `sign(v) sqrt(v^2 + 4)`
and whose eigenvector decays at `rho = (-|v| + sqrt(v^2 + 4)) / 2`.

## 1. The width is a division, not a search

There is no width available before the solve. The tail falls at
`decay_rate(H, lambda)` and `lambda` is what the solve produces; all that is
known in advance is that `lambda` lies outside the band and inside `||H||`, which
puts `rho` anywhere in `(0, 1)`. So the width has to be measured.

One measurement is enough, and the reason is the same one that made this class
first: outside the window the eigenvector is *exactly* geometric rather than
asymptotically so. The certificate's truncation term

    eta(n) = |b_inf| sqrt(u_{-n}^2 + u_n^2) / ||u||

therefore satisfies `eta(n') = eta(n) rho^(n' - n)`, and the width that brings it
to a target is

    n' = n + log(target / eta(n)) / log(rho)

which is a division. Measured against the law itself, over `n = 40` to `60`:

| `v` | `rho` | `eta(60)/eta(40)` | `rho^20` | ratio |
|---|---|---|---|---|
| 0.3 | 0.861187 | 5.0340e-02 | 5.0345e-02 | 0.9999 |
| 0.5 | 0.780776 | 7.0882e-03 | 7.0882e-03 | 1.0000 |
| 1 | 0.618034 | 6.1753e-05 | 6.6107e-05 | 0.9341 |
| 2 | 0.414214 | 2.3807e-01 | 2.2105e-08 | 1.08e7 |

The law holds to four digits while `eta` is above the level the computed
eigenvector stores its own tail at, degrades as it approaches it (`v = 1`, where
`eta(60)` is 1.6e-13), and is meaningless below it (`v = 2`, where both readings
are the 1e-16 plateau). That is the plateau of the step 3 note seen from the
other side, and it is why the search has to detect it rather than trust the
prediction.

What the division buys, in whole solves. Both columns start from the same first
width and stop at the same acceptance level; `applies` is the sum of
`fit$iterations` over the ladder.

| `v` | closed form | applies | doubling | applies |
|---|---|---|---|---|
| 0.2 | 41 -> 219 | 340 | 41 -> 82 -> 164 -> 328 | 760 |
| 0.3 | 41 -> 150 | 260 | 41 -> 82 -> 164 | 420 |
| 0.5 | 41 -> 93 | 200 | 41 -> 82 -> 164 | 300 |
| 1 | 41 -> 49 | 120 | 41 -> 82 | 120 |

Two widths at every rate, against two to four. The row where they tie is the one
where a doubling happens to land on the first try.

The landing is tight and it is tight for a reason. The search accepts at half the
tolerance and each prediction aims at an eighth of it, so the measured `eta` at
the chosen width should sit at about a quarter of the acceptance level. It sits
at 0.2264, 0.2186, 0.2206 and 0.2236 for `v` of 0.2, 0.3, 0.5 and 1. The margin
exists because `eta` is predicted from a vector that is itself only accurate to
its own residual; aiming at the acceptance level would miss about half the time
and spend a whole extra solve doing it.

## 2. The first width is a free tail, not a width

`n_start` counts from the edge of the window, not from the origin, so the first
section leaves the same free region whatever the potential. That is not a
convenience: `eta` is read at the section's edge and the rate is
`decay_rate(P, q)` with `q` the value that section produced, so a section which
barely contains the window has nowhere for the eigenvector to have become
geometric and predicts from a `q` that is not yet the right value.

Measured on a 121-site well of depth `-1`, radius 61:

| `n_start` reading | first width | widths |
|---|---|---|
| absolute, floored at `radius + 1` | 62 | 62 -> 80 -> 96 |
| free tail beyond the window | 101 | 101 |

Three widths against one. On a single-site potential the two readings differ by
one index (40 against 41), so nothing is paid for it where the window is small.

## 3. Three stops, and the common one is not the interesting one

The search stops on the target, on the plateau, or on `n_max`, and each is
reported in the certificate's dispatch line rather than retried.

The plateau branch is real and it is rare. Sweeping `tol` on `v = 1`:

| `tol` | `n` | widths | `eta` | why it stopped | certificate |
|---|---|---|---|---|---|
| 1e-12 | 59 | 41 -> 59 | 3.35e-13 | target | pass |
| 1e-14 | 68 | 41 -> 68 | 3.73e-15 | target | pass |
| 1e-15 | 73 | 41 -> 73 | 2.50e-16 | target | fail |
| 1e-16 | 78 | 41 -> 78 | 1.23e-16 | target | fail |
| 1e-17 | 114 | 41 -> 83 -> 93 -> 100 -> 107 -> 114 | 1.79e-16 | plateau | fail |
| 1e-18 | 97 | 41 -> 88 -> 97 | 7.43e-17 | plateau | fail |

At 1e-15 and 1e-16 the search still meets its own target: the truncation term can
be driven to 1e-16, and what fails the certificate is the *residual*, at 3.1e-15
and 5.1e-15. That is the outcome the split exists for. The search owns `n` and
reports it met its target, and the certificate points the reader at `tol` and
`ncv` instead. Only below about 1e-17 does the tail itself become the wall.

The 1e-17 row is where the closed form stops paying: near the plateau the ratio
law no longer holds, so each prediction asks for a small increment and six widths
go by before `eta` fails to improve. It is bounded and it is reported, and it is
the regime where the answer is already at the arithmetic floor.

**Correction to my own first draft.** At the plateau the `truncation bound` row
reads `qualified`, not `fail`. At `tol = 1e-18` an `eta` of 7.4e-17 is inside the
arithmetic floor `c eps ||H||` of 2.7e-15, so `cert_level()` returns `qualified`
and the detail says the tail is at the arithmetic floor. The test asserting
`fail` was wrong and the code was right; this is the S0.6 floor doing what it is
for, on a row it was not written for.

## 4. `certificate()` is a reader over results, and the sketch had it on an operator

The one-package sketch wrote

```r
F <- finite_section(A, n = 20)
eigs(F)
certificate(F)
```

The third line is refused. A certificate belongs to a result and not to an
operator: every claim a certificate makes is about a computation, and the one
claim an operator can make alone is that it satisfies the contract, which is
`verify()` and returns a certificate of its own. So `certificate.linop()` errors
and names `verify()`. That is the design the whole of Phase 2 already followed
(`README.md` calls it "the certificate that belongs to a result rather than to an
operator"); the sketch predates the shapes being settled.

What the accessor is for is that each verb puts its certificate where its own
return value has room, and those places are not alike: an attribute on a solution
so it stays the matrix a matrix solve returns, a field of a classed result for
the spectral verbs, the certificate itself from `verify()`, and a field of the
plain list `solve(details = TRUE)` returns. One reader over four conventions.

**A second correction, caught by a test rather than by reading.** The first
draft's error message said arithmetic on a solution drops the attribute. It does
not: `x + 0` and `2 * x` keep it, and indexing, `as.numeric()`, `c()` and the
reductions drop it. That is already in
`dev_notes/solve-dispatch-and-the-empty-branch.md` as a correction to plan
section 1.1, and it is where a user meets it, so the message names the four
operations that actually drop it.

## 5. What it cost

Nothing on the eigensolver's surface. `eigs()` on an operator with no matrix is
the same name and the same arguments; the only addition is `section =`, a list of
two knobs in the shape `inner =` and `norm_control =` already have. A fixed width
is not among them, because `eigs(finite_section(A, n))` is already how a caller
fixes it and a knob with two places to set it has one place too many.

The node registry gained `spectrum`, alongside `certify` from step 3: a node type
on a sequence space says how it becomes computable, and `eigs()` carries no
sequence-space vocabulary. `svds()` and `solve()` were deliberately not given the
hook. A singular value of a self-adjoint operator on `l^2(Z)` is not a second
question this layer answers, and a solve against one has no right-hand side to be
truncated with; both still refuse by name.

`certificate` is the twenty-fifth export and the cap is now exactly met. The
budget test asserts equality as well as the bound, so the next name in either
direction is a decision with a note behind it.

## Corrections table

| Where | Was | Is |
|---|---|---|
| One-package sketch | `certificate(F)` on a finite section | Refused. A certificate belongs to a result; an operator's own claim is `verify()` |
| First draft of the plateau test | `truncation bound` reads `fail` | `qualified`. The arithmetic floor covers a tail at 7.4e-17 |
| First draft of the accessor's error | arithmetic drops the attribute | It survives arithmetic; indexing, `as.numeric()`, `c()` and reductions drop it |
| First draft of `n_start` | an absolute first width | A free tail beyond the window. Three widths against one on a 121-site well |
