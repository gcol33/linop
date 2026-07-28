# The benchmark harness, and the fixture that solved itself

Written 2026-07-28. Gate 2's last line: a benchmark harness that runs end to end with
committed results. It lives in `dev_notes/bench/`, runs in 2.5 minutes, and writes four CSVs,
an environment stamp and `RESULTS.md`, all committed. Six findings. Two of them are about the
harness catching itself, which is the reason to write one down rather than run it once.

```
Rscript dev_notes/bench/run-bench.R      # from the repository root
```

## 1. Two fixtures had drawn their right-hand side from their own seed, and both were degenerate

The first version of `bench-solvers.R` built each case by calling the section 10 generator with
a seed and then seeding the right-hand side with the same number, which reads as tidy and is
the whole bug. Both generators build their operator out of a QR factorisation of a random
matrix, `matrix(rnorm(n * n), n, n)` fills column-major, and `qr()` gives `A = QR` with
`A[, 1] = Q[, 1] R[1,1]`. So the first draws from a seed are the first column of the input,
and that column is a multiple of the first column of `Q`.

What that did to each fixture:

| fixture | seed for `b` | result |
|---|---|---|
| `spd_prescribed(300, kappa 1e4)` | 5, the fixture's own | `b` is the eigenvector of `lambda_1`; **CG converges in 1 iteration** |
| `spd_prescribed(300, kappa 1e4)` | 105 | 685 iterations |
| `lsq_prescribed(600 x 200)` | 7, the fixture's own | `\|\|(I - U U^H) b\|\| / \|\|b\|\| = 0.0000`; the incompatible case is compatible |
| `lsq_prescribed(600 x 200)` | 107 | 0.7939 outside the range |

Neither failure announces itself. The one-iteration CG returns the right answer to 1.5e-13 and
certifies `qualified`; the compatible least-squares problem returns the right answer, and the
certificate correctly reports the *compatible* reading, because that reading was correct. The
table would have been a page of plausible numbers describing nothing, and the only visible
symptom was an iteration count that nobody had a prior for.

Two rules came out of it, both now in `bench-solvers.R`:

- A right-hand side never draws from the seed its fixture drew from.
- A fixture whose point is a property of `b` measures that property and reports it. The
  rectangular case computes `||(I - U U^H) b|| / ||b||` and carries it in the `note` column of
  every one of its rows, so the table itself says the right-hand side is outside the range.

The second rule is the one that generalises. The suite already has the discipline for knobs:
before asserting a knob helped, check it changed the iterates. This is the same discipline one
step earlier, applied to the fixture instead of the knob.

## 2. ncv is the binding knob for both spectral verbs, and maxit and tol are not

`eigs(laplacian_1d(400), k = 6)` at `ncv = 40` stops after 240 iterations with nothing
converged, a quarter of the way into a budget of 1000. The obvious reading is that the stall
detector gave up early. It did not:

| which | ncv | maxit | tol | iterations | nconv | value error |
|---|---|---|---|---|---|---|
| largest | 20 | 1000 | 1e-10 | 260 | 0/6 | 9.99e-04 |
| largest | 40 | 1000 | 1e-10 | 240 | 0/6 | 7.44e-04 |
| largest | 40 | 5000 | 1e-10 | 240 | 0/6 | 7.44e-04 |
| largest | 40 | 5000 | 1e-08 | 240 | 0/6 | 7.44e-04 |
| largest | 80 | 5000 | 1e-10 | 477 | 6/6 | 1.11e-16 |
| smallest | 40 | 5000 | 1e-10 | 160 | 0/6 | 4.40e+00 |
| smallest | 80 | 5000 | 1e-10 | 480 | 6/6 | 5.30e-15 |

Five times the budget buys the same iteration count, the same restart count and the same value
error to every digit printed, and so does loosening the tolerance by two orders. The subspace
is the only knob that moves anything, and the stall detector is reporting a fact about it:
there was nothing left in that subspace to find. `dev_notes/spikes/eigs-ncv-probe.R` runs the
whole table above; the harness carries the finding by running `ncv = 20, 40, 80` on both verbs,
so the transition is in the committed results.

Both ends of this spectrum are equally hard and for the same reason. The relative gap is
4.6e-5 at the top and 4.6e-5 at the bottom, which is `laplacian_1d` being symmetric about its
midpoint, and the two ends behave identically at every setting above.

The `smallest` rows are the ones to read for what an unconverged answer costs. At `ncv = 40`
the value error is 4.4: the returned values are wrong by a factor, not by a digit, because the
bottom of this spectrum is `6.1e-05, 2.5e-04, 5.5e-04` and a Ritz value that has not converged
lands among the wrong ones. `nconv` is the column that says so, and it is why `nconv` is in
`RESULTS.md` ahead of the error.

## 3. The sparse column wins at every size on this ladder, and printing it is the point

`tridiag(-1, 3, -1)` through its action, through a dense matrix, and through `Matrix`:

| n | matrix free, cg | sparse | dense, form | dense, LU |
|---|---|---|---|---|
| 500 | 0.0026 s | 0.00045 s | 0.0003 s | 0.011 s |
| 4000 | 0.0084 s | 0.0015 s | 0.072 s | 5.31 s |
| 8000 | 0.013 s | 0.0024 s | 0.128 s | 43.3 s |
| 16000 | 0.022 s | 0.0043 s | not run, 1.9 GB | - |
| 100000 | 0.111 s | 0.028 s | not run, 74.5 GB | - |

A tridiagonal matrix has an O(n) direct solve and a sparse library that implements it, so the
sparse route is faster than the iterative one at every size here, by a factor of four to five.
That number belongs in the committed table, under section 15's rule about reporting both and
letting readers compare, and it is also the honest statement of what the package claims: a
matrix-free method is not a faster route to a banded system, it is the route to an operator
that has no matrix. The dense column is where that starts to be the difference between 43
seconds and 13 milliseconds, and then between 74.5 GB and 13 milliseconds.

The shift is load-bearing in the fixture. The plain Laplacian's condition number grows like
`n^2`, so a size ladder on it measures conditioning while claiming to measure size; shifted,
the spectrum sits in `(1, 5)` at every `n` and CG takes 20 iterations at every rung of the
ladder, from 500 to 100000. The flat iteration count is in the committed table, so the fixture
is checked rather than asserted.

## 4. gc()'s max used is churn, and reporting it as peak memory would have been wrong by four orders

The first draft carried a `peak_mb` column from `gc(reset = TRUE)` around each solve. It reads
112 Mb for a CG solve at `n = 1000` whose entire working set is a handful of 8 kb vectors,
because R's max-used counter accumulates allocation between collections rather than tracking
what is live. A reader would have concluded that a matrix-free solve needs a hundred megabytes.

Memory in the harness is `object.size()` of what a route holds, which is exact: 0.048 Mb for
the operator at every `n` on the ladder against 488 Mb for the dense matrix at `n = 8000`. The
column that could not be measured honestly was removed rather than qualified.

## 5. The lockstep identity holds outside the test suite, and the two costs move opposite ways

| k | lockstep applies | separate applies | lockstep seconds | separate seconds | agreement |
|---|---|---|---|---|---|
| 1 | 698 | 698 | 0.066 | 0.069 | 0 |
| 4 | 698 | 2775 | 0.124 | 0.263 | 0 |
| 16 | 708 | 11155 | 0.366 | 1.05 | 0 |

On `spd_prescribed(300)`. The apply count is flat in `k` and the answers agree to 0.0e+00, at
every width, on both the stored and the matrix-free fixture: the same bitwise identity the CG
suite asserts, reached here by a different route. The `k = 1` row is the control, where the two
routes are the same computation and the numbers agree.

What the block does not save is column-applies, and the harness reports both because they move
in opposite directions: 10975 columns through the block against 11155 through the separate
solves at `k = 16`, because the block carries every column until the last one converges. The
saving is in calls, not in arithmetic, which is why the stored fixture (where a block apply is
one GEMM instead of `k` GEMVs) and the matrix-free one (where it is one R closure call instead
of `k`) are both in the table.

## 6. Rows worth knowing about, none of which are new claims

- Restarted GMRES(30) exhausts its budget on both stiff SPD fixtures, `laplacian_1d(200)` and
  `spd_prescribed(300)`, and certifies `fail` at 4.6e-08 where CG certifies at 3.0e-16 in a
  tenth of the applies. Restarted GMRES stagnating on an SPD system is textbook, and the
  harness reproducing it is a check on the harness.
- LSQR and LSMR on a *square* SPD system square the condition number, and both fail at
  `spd_prescribed(300, kappa 1e4)` and `hpd_prescribed(200, kappa 1e3)`. They are on those
  rosters because naming a method is allowed, and their failing there is the shape of the cost.
- At the shared budget on the rectangular fixture, LSMR converges at 1977 iterations and LSQR
  does not converge at 2000, with backward errors 9.95e-09 and 1.07e-07. Consistent with
  `dev_notes/lsmr-and-the-monotone-backward-error.md`, measured here on one seed rather than
  twelve, so it is an illustration of that finding and not a second one.
- `norm2()` on `convdiff_1d`, which supplies no adjoint, takes the probe route at **1 apply**
  where the power route costs 10 to 12. Every route returns a lower bound, so the certificate
  stays sound and its floor is simply looser.
- The certificate picks the compatible reading on `diff_1d(400)`, a rectangular operator,
  because `D` is onto and every right-hand side really is compatible. Which reading applies is
  a fact about `b` and the range of `A`, decided by measurement, and the table shows both
  readings occurring on rectangular fixtures.

## What the harness deliberately does not do

No cross-package comparison. `sanic`, `eigencore`, `RSpectra` and `irlba` all solve overlapping
problems, and a committed table putting linop beside them is a different artefact with
different rules: section 15 governs it, `dev_notes/eigencore-audit.md` records that nothing in
eigencore's `certification.R` or `validation.R` has been read, and the plan asks that
Buchsbaum be approached with an interface proposal before any comparison is published. The
routes the harness does compare against are the ones already in the R a user has open: base
`solve()`, `eigen()`, `svd()`, and `Matrix`.

The counted operators declare their capabilities rather than establishing them by computation,
because an operator with a counter in it has to be a callback leaf. So `method = "auto"` is not
what these tables measure and every row names its method. What `auto` chooses is
`dev_notes/solve-dispatch-and-the-empty-branch.md`.
