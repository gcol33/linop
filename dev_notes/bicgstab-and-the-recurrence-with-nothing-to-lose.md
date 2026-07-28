# BiCGSTAB, and the only recurrence here that agrees with its reference bitwise

Written 2026-07-28. The seventh and last of the Krylov methods of section 7.1, and the
second that requires nothing of the operator. Three findings, one of which is about the
other six methods rather than about this one.

## 1. It agrees with the published recurrence bitwise, at every step count

Every other method in the package has a ceiling on how far a reference test can be carried.
CG's is the highest, MINRES and GMRES hold for a while, and LSQR and LSMR lose agreement
with a second implementation of their own recurrence by four steps. BiCGSTAB has no ceiling:

| fixture | k=1 | k=2 | k=4 | k=8 | k=12 | k=16 | k=24 |
|---|---|---|---|---|---|---|---|
| convdiff mu=.3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| convdiff mu=.7 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| laplacian | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| kms rho=.7 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| complex dense | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Exactly zero, not "below 1e-15": `expect_identical` on the iterates.

The reason is worth stating because it explains the other five rows of that comparison.
**What drifts in a Krylov method is orthogonality, and orthogonality is a property of a
stored basis.** GMRES stores one and orthogonalises against all of it. LSQR and LSMR store
two vectors that are supposed to stay orthogonal to everything that came before and do not.
BiCGSTAB stores `r`, `p`, `v`, `s`, `t` and a shadow vector, orthogonalises against nothing,
and asks no vector to remain orthogonal to any other. There is no invariant to lose, so
there is nothing for two implementations to disagree about: the same arithmetic operations
happen in the same order and produce the same bits.

That is not a claim that BiCGSTAB is more accurate. It is the opposite kind of fact: the
method has no orthogonality to lose because it never had any, and what it gives up in
exchange is section 2.

## 2. What it gives up is the minimisation

GMRES minimises the residual over the whole Krylov space, so its true residual cannot rise,
and its own suite asserts that. BiCGSTAB minimises over the single direction of the
stabilising half, which is one step of GMRES(1) on the intermediate residual `s`. Over 8
seeds and 20 steps:

| fixture | steps where the true residual rose | worst single-step rise |
|---|---|---|
| convdiff mu=.3 | 4 to 7 | +1171% |
| convdiff mu=.7 | 5 to 10 | +4079% |
| laplacian | 0 to 2 | +4% |
| kms rho=.7 | 0 to 1 | +590% |
| complex dense | 0 to 4 | +480% |

A factor of 41 in one step on the nonsymmetric fixture. This is the strongest case in the
package for the outer loop measuring rather than trusting: a stopping test reading only the
recurrence could stop on a step whose residual had just risen by that factor.

Against the cost, on the same fixtures and the same tolerance:

| fixture | BiCGSTAB | GMRES(30) |
|---|---|---|
| convdiff mu=.3 | 40 iterations, 80 applies | 207 iterations |
| laplacian | 40 iterations, 80 applies | 299 iterations |
| kms rho=.5 | 19 iterations, 38 applies | 28 iterations |
| shifted laplacian | 76 iterations, 152 applies | 1743 iterations |

Two applies per BiCGSTAB step against one per GMRES step, and no basis, no
orthogonalisation and no restart parameter. Both are in this package; the table is here so
the trade is on the record, not to rank them.

## 3. Breakdown is not an error, and one system makes that concrete

BiCGSTAB has three divisions by an inner product, and each is a breakdown when the product
goes to zero: `rho = <rhat0, r>` is the BiCG breakdown, `<rhat0, v>` is the alpha breakdown,
and `<t, s>` going to zero leaves `omega` at zero and the `p` recurrence dividing by it.

The threshold was measured on operators where no breakdown is present, per the MINRES rule
about thresholding a quantity rounding moves. Smallest relative value each reached over 12
seeds:

| fixture | rho | `<rhat0, v>` | `<t, s>` |
|---|---|---|---|
| convdiff mu=.3 | 2.8e-05 | 3.4e-04 | 9.5e-02 |
| convdiff mu=.7 | 1.2e-04 | 1.7e-05 | 9.1e-02 |
| laplacian | 1.0e-08 | 1.5e-08 | 8.8e-02 |
| kms rho=.7 | 1.1e-05 | 1.3e-04 | 3.7e-01 |
| shifted laplacian | 8.8e-12 | 2.3e-11 | 1.6e-04 |
| complex dense | 7.0e-03 | 1.5e-02 | 2.8e-01 |

`BICGSTAB_BREAKDOWN_TOL = 1e-14` is about 900 times below the closest approach a healthy
solve made, on a nearly singular indefinite operator, and six orders below it everywhere
else. It is also about 45 eps, which is where a relative inner product stops carrying
information at all, so the threshold sits at the arithmetic floor rather than at a level
tuned to these fixtures. A test asserts it is bitwise inert: with `breakdown_tol = 0` the
iterates and the iteration counts are identical on every healthy fixture.

**The cure is the restart, and the restart is the outer loop.** A column that breaks down
keeps the iterate it had and leaves. The loop then re-measures `b - A x` and starts a round
with the shadow vector re-seeded to that residual, which makes `rho_1 = ||r||^2` and cannot
break down at its first step. That is the standard cure, and it costs no new machinery: the
round structure was already there for GMRES's restarts.

A cure that recovers nothing is reported. For a real skew-symmetric `A`, `<A z, z> = 0` for
every real `z`, so `<rhat0, v>` is exactly zero at the first step of every round however the
shadow vector is re-seeded. BiCGSTAB cannot solve that system. What comes back is a `fail`
certificate whose convergence row says the recurrence broke down and a fresh shadow residual
did not recover it, a finite iterate, and no error thrown. GMRES solves the same system to
3.8e-16. That pair is the test: the two methods that require nothing of the operator are not
interchangeable, and the one that fails says so rather than returning a number.

## 4. The three sides, and the trap in testing them

Section 4.3 leaves BiCGSTAB unrestricted, and as with GMRES that is not a statement that the
sides agree. Right is Saad's Algorithm 7.7, on `A M^-1`, and its `r` is the true `b - A x`.
Left runs on `M^-1 A` and measures `M^-1 r`. Split runs on `L^-1 A L^-H` and measures
`sqrt(<r, M^-1 r>)`; both of the latter inherit MINRES's currency conversion.

Split would need three further preconditioner applies per step if `<u, M^-1 v>` were formed
on demand. It needs none. Under the uniform change of variable `w~ = L^-1 w`, every vector
of the split iteration maps to the corresponding vector of the *right* iteration and the
euclidean inner product maps to `<u, M^-1 v>`, so split is the right form with one
substitution. And `M^-1 r`, `M^-1 s` and `M^-1 p` satisfy the same linear recurrences their
unpreconditioned counterparts do, so carrying them alongside leaves only `M^-1 v` and
`M^-1 t` to be applied. Two applies per step on every side. This is the trick GMRES uses to
carry `u_j = M v_j`, one level up.

**The first draft of the test that the three differ proved nothing, and the reason is in
CLAUDE.md already.** The Jacobi preconditioner of `convdiff_1d` is built from its diagonal,
which is constant 2, so `M` is a scalar multiple of the identity, and a scalar commutes with
everything. All three sides produced *bitwise identical* iterates and the test passed while
asserting nothing. With a non-constant `M` the three differ by 8.1e-1, 1.7e-1 and 7.6e-1 at
three steps. The test now asserts both halves: that the three differ under a varying `M` and
that they coincide bitwise under a scalar one. The second assertion is what stops the first
from decaying back into the version that proved nothing.

## 5. Shared

`split_precond_not_hpd()` and `left_precond_singular()` moved from `solvers-gmres.R` to
`solvers-common.R` and take the method name. They were written for GMRES and are needed
verbatim by the second method that implements all three sides; the split message also takes
the quantity that gave the preconditioner away, which is `<w, M^-1 w>` for GMRES and
`<r, M^-1 r>` here.

## Files

- `R/solvers-bicgstab.R`, `R/solvers-common.R`, `R/solvers-gmres.R`
- `tests/testthat/test-solve-bicgstab.R`, 28 tests, 292 assertions
- `dev_notes/spikes/bicgstab_probe.R`, `bicgstab_reference.R`
