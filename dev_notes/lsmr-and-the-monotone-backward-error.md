# LSMR, and the estimate that cannot be reported because the method minimises it

Written 2026-07-28. The sixth of the seven Krylov methods, and the second on the rectangular
least-squares problem. It shares the Golub-Kahan bidiagonalisation with LSQR and differs only
in what it minimises over the space that bidiagonalisation spans, so the question this note
answers is whether that difference earns a method or is a flag on the last one.

It earns a method, and the reason is a property of the certificate rather than of the
recurrence.

## 1. The certificate reports the quantity this method minimises

`solve_certificate(least_squares = TRUE)` reports `||A^H r|| / (||A|| ||r||)` as the
least-squares backward error, because Stewart's perturbation exhibits a `dA` of exactly that
relative size for which `x` is the exact minimiser. LSQR minimises `||r||` over the Krylov
space. LSMR minimises `||A^H r||`, which is the numerator of the number the certificate will
print.

So LSMR's progress is monotone in the reported quantity and LSQR's is not. Measured densely
from the iterates rather than read off either recurrence, over 20 steps, 12 seeds and three
conditionings:

| kappa | LSMR rises, worst seed | worst relative rise | LSQR rises, min and max over seeds |
|---|---|---|---|
| 1e2 | 0 | 0 | 5 to 9 |
| 1e4 | 0 | 0 | 6 to 9 |
| 1e6 | 0 | 0 | 7 to 10 |

Not "rarely rises": exactly zero rises in 36 runs of 20 steps, against a LSQR sequence that
rose on every seed, by up to a factor of 60 in one step. The worst single-step rises measured
for LSQR were +1076%, +1014%, +1274%, +4348%, +3807% and +5993% across six seeds.

The practical form of that is what a shared budget buys. Same operator, same right-hand side,
same `maxit`, same `tol = 0`, compared on the certificate's own number over 12 seeds:

| kappa | budget | LSMR smaller on | ratio LSQR/LSMR, min | median |
|---|---|---|---|---|
| 1e4 | 20 | 12 of 12 | 1.99 | 4.29 |
| 1e4 | 40 | 12 of 12 | 1.21 | 10.41 |
| 1e4 | 80 | 11 of 12 | 0.46 | 5.14 |
| 1e6 | 40 | 12 of 12 | 1.09 | 1.55 |
| 1e6 | 80 | 12 of 12 | 49.05 | 1768.92 |
| 1e8 | 40 | 12 of 12 | 8.94 | 122.63 |
| 1e8 | 80 | 12 of 12 | 105.91 | 2401.83 |

The one row where LSQR came out ahead on a seed is the one to keep: the claim the test
asserts is that LSMR's number is smaller, not that it is smaller by any particular factor,
and the test uses `kappa = 1e8, budget = 40` where the margin held on every seed.

Run to convergence at the same tolerance the two land in the same place at roughly the same
cost, so this is a statement about a truncated budget rather than about the answer.

## 2. The estimate the method minimises is the one that cannot be reported

`|zetabar_{k+1}|` is `||A^H r_k||` exactly, for the projected problem. It is what LSMR's
stopping rule reads, it costs nothing, and it is the sharpest available quantity right up to
the point where it stops being true. On a 15-column problem with `kappa` about 40:

| step | estimated `\|\|r\|\|` | true `\|\|r\|\|` | estimated `\|\|A^H r\|\|` | true `\|\|A^H r\|\|` |
|---|---|---|---|---|
| 12 | 5.22001e+00 | 5.22001e+00 | 3.25589e-01 | 3.25589e-01 |
| 18 | 5.20614e+00 | 5.20614e+00 | 2.62265e-03 | 2.62265e-03 |
| 20 | 5.20614e+00 | 5.20614e+00 | 6.50361e-07 | 6.50361e-07 |
| 25 | 5.20614e+00 | 5.20614e+00 | 2.80861e-13 | 2.92858e-13 |
| 30 | 5.20614e+00 | 5.20614e+00 | 1.31427e-16 | 1.09916e-13 |

At step 30 the recurrence believes the backward error is 800 times smaller than it is, and
believes it in the direction that would certify a solve as better than it is.

This is the CG rule at its sharpest so far. CG's version was that the recurrence residual
drifts where the answer is most converged; LSQR's was that two implementations of the same
recurrence disagree by O(1) mid-solve. Here the drifting quantity is the one the certificate
prints, and the method's whole advantage is that it minimises it. **The advantage is real and
the recurrence's measurement of it is not**, so the outer loop's extra apply per round is
what makes the advantage reportable at all. Nothing in section 1 above was measured from
`zetabar`.

## 3. The reproducibility ceiling, against a third party this time

The LSQR note recorded two arithmetically equivalent implementations diverging three orders
of magnitude every two steps. That was two implementations written here. LSMR was checked
against SciPy 1.17.1's `scipy.sparse.linalg.lsmr`, which is Fong and Saunders' own
transcription, at a fixed step count with every stopping test switched off, so what is
compared is the recurrence:

| steps | kappa 1e2 | kappa 1e6 | kappa 1e10 | 200 x 30, kappa 1e4 |
|---|---|---|---|---|
| 1 | 8.3e-15 | 1.9e-15 | 1.4e-15 | 1.2e-15 |
| 4 | 6.9e-15 | 2.0e-14 | 6.1e-12 | 7.5e-15 |
| 8 | 2.6e-13 | 2.6e-08 | 4.7e-08 | 6.9e-11 |
| 16 | 1.3e-04 | 8.8e-01 | 7.1e-07 | 2.0e-07 |
| 32 | 2.6e-06 | 4.1e-03 | 2.9e-08 | 1.3e-01 |
| 64 | 1.2e-14 | 2.7e-05 | 1.4e-02 | 2.8e-02 |

Same schedule, and against an implementation nobody here wrote. The ceiling is the method's,
not the code's. A reference test can assert four steps.

Against an R transcription of the paper rather than SciPy the four-step gap is tighter and
tracks the conditioning, which is what the ill-conditioned Gate 2 test's tolerances are set
from: 1.2e-15 at `kappa` 1e2, 3.8e-14 at 1e4, 1.2e-13 at 1e6, 6.5e-13 at 1e8, 1.9e-11 at
1e10, 2.1e-10 at 1e12.

Finite termination goes the same way as LSQR's. At step 15 of a 15-column problem the
published recurrence is wrong in the second digit (4.7e-2); it reaches 1.9e-6 at step 18,
1.9e-10 at 20 and 6.4e-15 at 25. About `1.7 n`, as before.

## 4. A converged solve, and where the forward error goes

Running to convergence at `tol = 1e-12` against the closed form, with SciPy on the same
fixtures at the same tolerance:

| fixture | this LSMR | SciPy LSMR | this LSQR | SciPy LSQR |
|---|---|---|---|---|
| kappa 1e2 | 39 it, fwd 7.7e-13 | 38 it, fwd 3.1e-12 | 39 it, fwd 7.7e-13 | 38 it, fwd 3.0e-12 |
| kappa 1e6 | 152 it, fwd 8.5e-09 | 131 it, fwd 6.4e-09 | 157 it, fwd 1.6e-10 | 131 it, fwd 6.4e-09 |
| kappa 1e10 | 448 it, fwd 9.8e-03 | 291 it, fwd 9.1e-01 | 457 it, fwd 7.5e-05 | 291 it, fwd 9.2e-01 |
| 200 x 30 | 178 it, fwd 1.7e-11 | 163 it, fwd 2.1e-09 | 178 it, fwd 1.7e-11 | 163 it, fwd 2.1e-09 |

Two things to read out of the last row of numbers, and one not to.

**The iteration counts differ because the stopping tests read different things.** SciPy stops
when its recurrence estimates say the tolerance is met, at 291 steps on the `kappa` 1e10
fixture. The loop in `solvers-bidiag.R` re-measures `b - A x` and `A^H r` once per round and
keeps going to 448. Section 2 is why: at that point the recurrence's own estimate has fallen
below the truth. This is the same design decision stated three times now, and it is the first
time its cost and its benefit are both visible in one table.

**The forward errors at `kappa` 1e10 are not a difference in accuracy delivered.** All four
runs certify a backward error between 6e-9 and 2.4e-8. For a least-squares problem with an
incompatible right-hand side the forward-error bound carries `kappa^2` on the residual term,
and `kappa^2 eps` is 2.2e+04 there: no backward error either method achieves constrains the
forward error at all. The spread is which nearby problem each landed on. That is exactly why
the certificate's forward-error row says `not_checked` and names `||A^-1||` as what is
missing, and this table is the first live demonstration that the restraint is right rather
than merely cautious.

**What not to read out of it** is that LSMR is less accurate than LSQR. Its forward error was
larger on the `kappa` 1e10 fixtures here and equal or smaller elsewhere, and SciPy's two
methods agree with each other to within 1% on the same fixture. At a conditioning where
`kappa^2 eps` exceeds 1 the forward error is not a property either method controls.

## 5. What was shared, and what was not

`R/solvers-bidiag.R` holds three things: `bidiag_solve()`, the loop, which is identical for
both methods and takes the recurrence as an argument; `bidiag_start()`, which forms
`beta_1 u_1 = r_0` and `alpha_1 v_1 = M^-H A^H u_1` and drops the columns where `alpha_1`
is zero; and `bidiag_step()`, one bidiagonalisation step including the two lower bounds on
`||A||` that fall out of applies the step was making anyway.

`bidiag_step()` returns `v_{j+1}` unnormalised, because the two recurrences reach the point
where they need it at different places and both need its norm before they need it.

What is not shared is the projection, which is the method. LSQR carries one triangular
factor, one direction `w` and `phibar`; LSMR carries two factors, two directions `h` and
`hbar`, and a third rotation sequence whose only purpose is to recover `||r_k||`, which its
own projected problem does not carry. The `||r||` estimate is Fong and Saunders section 3.2
and touches no vector, so it costs no apply.

The condition test is the one LSQR already had, on the extreme diagonals of the first
triangular factor. `rho` is the same sequence in both methods, and `max|rho| / min|rho|` is a
lower bound on `cond` of the QR factor of `B_k`, whose singular values lie inside those of
`A`. Fong and Saunders estimate `cond(A)` from the *second* factor instead; that estimate was
not used, because the first one is shared with LSQR, is already justified in
`solvers-common.R`, and bounds `cond(A)` directly.

Refactoring LSQR onto the shared file was checked for bitwise neutrality before LSMR was
written: all 30 tests and 120 assertions of `test-solve-lsqr.R` unchanged.

## 6. Gate 2

The plan's Gate 2 asks that MINRES and LSMR agree with a trusted reference on ill-conditioned
problems including breakdown and near-breakdown. Three tests cover it:

- agreement with an R transcription of Fong and Saunders at four steps, at `kappa` 1e4, 1e8
  and 1e10, five seeds each, with the tolerance tracking the conditioning rather than being
  one constant;
- the SciPy cross-check above, which is not in the suite because it would be a Python
  dependency, and lives in `dev_notes/spikes/lsmr_scipy.py` with its driver;
- exact breakdown as its own case rather than as a perturbation of one: `b = A v_1` for a
  single right singular vector makes `alpha_2` exactly zero, the Krylov space closes at step
  one, and the solve terminates there with a relative error of 1.3e-15.

## 7. Not included

`damp`, the Tikhonov parameter, for the reason it was left out of LSQR: it changes the
problem being solved, so the certificate would have to certify the damped problem and say so,
and the row in section 4.3 would need re-reading against a different normal equation. It is
a deliberate absence, not an oversight, and the same absence in both methods.

## Files

- `R/solvers-bidiag.R`, `R/solvers-lsmr.R`, `R/solvers-lsqr.R`
- `tests/testthat/test-solve-lsmr.R`, 29 tests, 141 assertions
- `dev_notes/spikes/lsmr_reference.R`, `lsmr_reference_fn.R`, `lsmr_illcond.R`,
  `lsmr_thresholds.R`, `lsmr_export.R`, `lsmr_scipy.py`, `lsmr_converged.R`,
  `lsmr_converged.py`
