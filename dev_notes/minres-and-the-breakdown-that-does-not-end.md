# MINRES against a reference, and the breakdown that does not end

Written 2026-07-29. Gate 2 asks that MINRES and LSMR agree with a trusted reference on
ill-conditioned problems, including breakdown and near-breakdown cases. LSMR cleared that
line against SciPy 1.17.1 in `lsmr-and-the-monotone-backward-error.md`. This is the MINRES
half. It is measured by `dev_notes/spikes/minres_reference.R` and
`dev_notes/spikes/minres_export.R` with `minres_scipy.py`, and asserted by the
`Gate 2` and `breakdown, and past it` sections of `tests/testthat/test-solve-minres.R`.

Four things came out of it.

## 1. Two references, sharing no code, and the ceiling belongs to neither

The suite carries two references and they are different kinds of object.

`krylov_argmin()` is the definition: an orthonormal basis for `K_m(A, b)` by modified
Gram-Schmidt with a second pass, then the minimiser of `||b - A x||` over that basis by
truncated SVD. It says what MINRES is *for* and nothing about how MINRES works.

`reference_minres()` is the published recurrence, transcribed from the Paige and Saunders
form with the preconditioner set to the identity and no shift. It says how MINRES works and
nothing about what it is for. It is validated against SciPy 1.17.1's `minres` by
`minres_export.R`, which is what separates a reference from a second copy of the same
beliefs: over four conditionings and both breakdown fixtures it tracks SciPy exactly as
closely as this implementation does.

| steps | kappa 1e2 | kappa 1e4 | kappa 1e6 | kappa 1e10 |
|---|---|---|---|---|
| 4  | 5.8e-15 / 5.3e-15 | 1.0e-14 / 1.0e-14 | 5.8e-15 / 4.4e-15 | 3.2e-15 / 2.9e-15 |
| 8  | 6.6e-15 / 6.9e-15 | 1.4e-14 / 1.6e-14 | 5.5e-13 / 5.1e-13 | 3.6e-12 / 4.3e-12 |
| 12 | 1.0e-14 / 1.3e-14 | 6.2e-13 / 1.8e-13 | 1.0e-08 / 9.3e-09 | 3.8e-01 / 1.5e-01 |

Each cell is `this implementation / the transcription`, both against SciPy, relative in the
iterate.

The first reading of that table is that MINRES is reproducible to about eight steps, where
LSQR and LSMR manage four. The second reading is the one that matters. Run the two
references against *each other*, with this package out of the picture entirely:

| kappa | 2 steps | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|
| 1e2  | 1.4e-15 | 4.7e-15 | 2.2e-15 | 3.8e-14 | 9.0e-12 |
| 1e6  | 1.1e-15 | 1.4e-15 | 2.8e-13 | 2.0e-08 | 1.5e+00 |
| 1e10 | 8.4e-15 | 7.7e-15 | 2.1e-11 | 1.3e+00 | 1.3e+00 |

A variational characterisation and a short recurrence, no shared code, no shared idea, and
they diverge on the same schedule that each of them diverges from this implementation. So
the ceiling is the method's. A reference test cannot be extended past it by tightening
anything, because past it there is nothing to tighten towards: neither object is ground
truth for the other. That is a stronger statement than LSQR or LSMR could make, where the
second implementation was written here and the schedule could in principle have been a
property of the pair.

What survives the ceiling is the answer. Run to convergence rather than to a step count, the
forward error sits between 0.5 and 22 times `kappa * eps` across kappa 1e2 to 1e8, which is
what the conditioning allows and what the test asserts.

## 2. An exact breakdown is one collapsed off-diagonal, not the end of the sequence

`beta_{j+1} = 0` is the Lanczos breakdown: the Krylov space is exhausted, the iterate is the
exact solution, and in exact arithmetic there is no step `j+1`. Putting `b` inside a
`d`-dimensional invariant subspace of the shifted Laplacian, whose eigenvectors are closed
form, constructs that exactly. The sequence `beta_{j+1} / ||A v_j||`, from a dense Lanczos
with no reorthogonalisation:

```
d = 1:  8.3e-16  3.6e-01  8.1e-01  6.5e-01  6.0e-01
d = 2:  4.7e-02  8.9e-14  1.5e-01  5.3e-01  5.7e-01
d = 3:  1.3e-01  1.5e-01  7.3e-13  1.4e-01  3.2e-01
d = 5:  2.7e-01  4.0e-01  2.9e-01  2.7e-01  3.8e-11  7.3e-02
d = 8:  7.6e-01  5.7e-01  6.3e-01  7.8e-01  4.6e-01  2.9e-01  5.7e-01  9.6e-11  4.2e-02
```

The collapse is there, at exactly the right step, and then the next entry comes back at
O(1). Rounding in `A v` puts a component outside the invariant subspace, the iteration
amplifies it, and the space re-seeds itself. The depth of the collapse degrades with `d`,
from 8.3e-16 at `d = 1` to 9.6e-11 at `d = 8`, so a fixed threshold on `beta` tuned on a
small case misses a larger one.

Neither reference stops at `d` either. `krylov_basis()` runs to 60 columns on the `d = 5`
and `d = 8` fixtures, and SciPy's minres runs 9 iterations when asked for 16 on the `d = 3`
one. Three independent programs, three failures to notice a breakdown that is
mathematically exact.

## 3. Iterating past the breakdown improves the answer

Which makes the obvious repair the wrong one. Forward error at `maxit = d, d+1, ... d+8`,
with `tol = 0` so the steps are actually taken:

```
d = 1:  1.7e-14  1.8e-14  1.8e-14  1.8e-14 ...
d = 3:  2.1e-14  1.9e-14  1.9e-14  1.9e-14 ...
d = 5:  5.9e-13  6.1e-14  3.0e-14  1.3e-14  1.2e-14 ...
d = 8:  3.9e-12  1.7e-13  1.1e-14  9.8e-15  8.3e-15 ...
```

The iterate at the breakdown step is the exact solution in exact arithmetic and is wrong by
3.9e-12 at `d = 8`; the steps taken past it bring that to 9.7e-15. A solver that detected
the breakdown and stopped would return the worse answer, by a factor of 400.

So nothing in `solvers-minres.R` thresholds `beta` to detect a breakdown, and nothing should
be added. What stops the solve is the residual, which is already at rounding level at step
`d` and stays there. The `spent` guard in `minres_recurrence()` is not breakdown detection:
it exists to separate an exhausted `w` from a preconditioner that has contradicted its
positive-definiteness declaration, and it suppresses an error rather than a step.

This is the same shape as LSQR's finite-termination finding read from the other side. There,
the iterate at step `n` was wrong in the second digit and the solve needed about `1.7 n`
steps. Here the iterate at step `d` is wrong in the twelfth and the solve needs a few more.
Both say: the step at which the Krylov space is complete in exact arithmetic is not a step
at which anything may be asserted.

## 4. The near-breakdown window crosses two thresholds, and neither is where it looks

`b = v_7 + delta * v_23` sits inside a one-dimensional invariant subspace to relative
accuracy `delta`, so `beta_2` is `O(delta)`. Sweeping `delta` down through the guard at
`sqrt(eps) = 1.5e-08`:

| delta | beta2 / \|\|A v\|\| | iterations | forward error | residual |
|---|---|---|---|---|
| 1e-02 | 8.2e-03 | 2 | 2.5e-14 | 1.3e-14 |
| 1e-08 | 8.2e-09 | 2 | 2.6e-14 | 1.3e-14 |
| 1e-10 | 8.2e-11 | 2 | 2.5e-14 | 1.3e-14 |
| 1e-12 | 8.2e-13 | 1 | 4.4e-12 | 8.2e-13 |
| 1e-14 | 8.4e-15 | 1 | 6.1e-14 | 8.4e-15 |
| 0     | 1.6e-15 | 1 | 2.5e-14 | 1.6e-15 |

Crossing `sqrt(eps)` at `delta = 1e-08` changes nothing, because the guard it belongs to
does not terminate anything. What does change the iteration count is the tolerance: at
`delta = 1e-12` the residual after one step is 8.2e-13, which meets `tol = 1e-12`, so the
solve stops there. That is the only row in the sweep where the forward error is not at
rounding level, and 4.4e-12 is exactly what a residual of 8.2e-13 entitles it to on this
operator. Every row converges and certifies `pass`.

The definition agrees with the iterates through the whole window to 4.4e-12 at three steps,
at the same `delta` where a single step is just good enough to stop on. The test asserts
both references over the window at 1e-10.

## 5. What this changes in the package

Nothing in `R/`. The implementation was already right about all of it, and the value of the
exercise was finding out which of its behaviours are load-bearing:

- the residual test rather than a `beta` threshold as the stopping rule, per section 3
- the outer loop measuring `b - A x` rather than reading `phibar`, which is what makes the
  near-breakdown rows report a residual they have actually reached

`tests/testthat/test-solve-minres.R` gains the two references and seven tests. The existing
`krylov_argmin()` was real-only, absolute in its exhaustion threshold and used `qr.solve`;
it is now complex-aware, relative and truncated-SVD, so one routine serves the real and
complex cases and the breakdown fixtures at once. The pre-existing test that used it, the
minimal-residual identity at up to 12 steps on a kappa-8 spectrum, passes unchanged, and
section 1 says why 12 steps was safe there and would not be at kappa 1e4.

## 6. Gate 2

The line is met. MINRES agrees with a trusted reference on ill-conditioned problems, at five
conditionings from kappa 1e2 to 1e10 over five seeds, against two references at once, and
the breakdown and near-breakdown cases are both constructed rather than hoped for. The
external validation of the transcription against SciPy 1.17.1 lives in the spike rather than
in the suite, so the suite has no Python dependency, which is how the LSMR line was cleared
as well.

Reproducing it:

```
Rscript dev_notes/spikes/minres_reference.R                   # 1.7 s, both references
MINRES_XDIR=<dir> Rscript dev_notes/spikes/minres_export.R write
MINRES_XDIR=<dir> python  dev_notes/spikes/minres_scipy.py
MINRES_XDIR=<dir> Rscript dev_notes/spikes/minres_export.R compare
```
