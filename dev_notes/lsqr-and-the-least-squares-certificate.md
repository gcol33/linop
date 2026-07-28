# LSQR, and the certificate for a problem that has no exact answer

Written 2026-07-28. The fifth of the seven Krylov methods, and the first whose problem is a
minimisation rather than an equation. Four things follow from that, three of which the
CG-MINRES-GMRES template had no position on, and one of which corrects the plan.

## 1. The certificate had a square-system shape, and it would have failed a correct solve

`solve_certificate()` as CG left it tests `||b - A x|| / ||b||` against `tol` and reports
Rigal-Gaches as the backward error. For a rectangular `A` neither line survives: `b` need
not lie in the range of `A`, so `b - A x` does not go to zero, and its size is a fact about
the problem rather than a measure of the solve. Run unchanged, a converged least-squares
solution certifies as `fail` on both lines.

This is the S0.6 shape again, structural this time rather than arithmetic. S0.6 found a
converged solve certifying as `fail` because the bound kept decaying past machine epsilon;
here it would certify as `fail` because the quantity being tested was never going to reach
the tolerance, however exact the answer.

**What goes to zero is `A^H r`, and Stewart's perturbation makes that a backward error
rather than a heuristic.** With

```
dA = -(r r^H A) / ||r||^2,      A + dA = (I - P_r) A,
```

the perturbed residual is `r + P_r A x`, and

```
(A + dA)^H (r + P_r A x) = A^H (I - P_r) r + A^H (I - P_r) P_r A x = 0
```

identically, because `(I - P_r) r = 0` and `(I - P_r) P_r = 0`. So `x` is the *exact*
least-squares solution of `min ||b - (A + dA) x||`, and the size of the perturbation is

```
||dA||_2 = ||r r^H A||_2 / ||r||^2 = ||A^H r|| / ||r||
```

exactly. Relative to `||A||` that is `||A^H r|| / (||A|| ||r||)`, which is Paige and
Saunders' second stopping rule. The rule and the backward error are the same number: S2 is
not a convergence heuristic that happens to work, it is a backward-error criterion.

Two consequences for the certificate.

**The rows do not change.** Both readings are backward errors, so one row carries both and
the column decides which it is entitled to. A `linop` certificate has the same shape whatever
method produced it, which matters once `solve()` dispatches.

**The lower bound on `||A||` stays the right direction.** It sits in the denominator here as
it does in Rigal-Gaches, so a lower bound enlarges the reported error. The existing rule
carries over without re-reading.

The residual row becomes `not_checked` where the problem is incompatible, with the measured
distance and the count of columns it applies to in the detail. That is decided *by
measurement inside the certificate*, not from a flag the solver passes: the solver's opinion
about which regime it is in is exactly the kind of thing the outer loop exists not to trust.
A column that misses both readings still reports `fail`.

`test-solve-lsqr.R` checks Stewart's perturbation directly rather than trusting the
derivation: it forms `dA`, confirms `x` satisfies the normal equations of `A + dA` to 1e-8
relative, and confirms `sigma_max(dA) / sigma_max(A)` is at most the number the certificate
reported.

## 2. The floor works out to `c eps` again, and it is load-bearing again

The computed `A^H r` inherits the residual's own error through `||A|| floor_abs` and adds its
own `c eps ||A|| ||r||`, so against the denominator `||A|| ||r||` the floor is

```
c eps + floor_abs / ||r||.
```

The second term matters only where `||r||` has itself fallen to the residual floor, which is
the regime where the compatible reading is in use, so in practice the least-squares line
carries the same plain `c eps` the backward-error line already carried.

Measured on a 40 x 12 problem with singular values 1 to 4, at a tolerance set to half what
the solve actually delivered:

```
achieved backward error 1.5954e-16, c eps = 8.8818e-16
identical iterates: TRUE
backward error   with floor qualified  without fail
convergence      with floor pass       without fail
```

Same iterates, byte for byte, certified twice. Without the floor a fully converged
least-squares solve reports `fail`.

The test sets its tolerance from what the run achieved rather than from a constant. The
window the floor covers is one `c eps` wide, and whether a fixed tolerance lands inside it
depends on the fixture; a fixed constant would have made the test pass or fail for reasons
that have nothing to do with the floor.

## 3. Two implementations of the recurrence do not stay close, and it is not a bug

The reference test was first written the way CG's is, asserting bitwise agreement with a
textbook single-column implementation over 12 steps. It failed by 1.6e-2. The recurrence was
right: agreement is at rounding level for the first few steps and then leaves.

Relative gap between this implementation and the published recurrence, same fixture, 45 x 15
with singular values 0.5 to 20, so `kappa(A) = 40` and nothing here is about ill conditioning:

|  | seed 1 | seed 2 | seed 3 | seed 4 | seed 5 |
|---|---|---|---|---|---|
| step 2 | 4.6e-16 | 2.8e-16 | 5.4e-16 | 1.7e-16 | 4.2e-16 |
| step 4 | 2.1e-15 | 5.7e-16 | 1.0e-15 | 5.5e-16 | 6.5e-16 |
| step 6 | 5.5e-14 | 1.4e-14 | 5.3e-14 | 2.0e-14 | 2.2e-15 |
| step 8 | 1.8e-11 | 1.0e-11 | 6.0e-11 | 5.3e-12 | 7.4e-13 |
| step 10 | 7.0e-08 | 6.4e-08 | 1.1e-07 | 2.7e-08 | 1.9e-09 |
| step 12 | **1.8e-02** | 9.8e-03 | 6.1e-03 | 5.7e-03 | 8.2e-05 |
| step 14 | 2.7e-05 | 1.8e-04 | 4.2e-04 | 1.1e-05 | 2.6e-05 |
| step 16 | 4.9e-04 | 1.7e-06 | 3.2e-04 | 3.7e-04 | 2.7e-05 |

Three orders of magnitude every two steps, a peak of O(1) relative around step 12, and then
it closes as both converge on the same answer. The mechanism is Paige's: the bidiagonal
vectors lose orthogonality, and what fills the lost direction is whichever rounding noise
arrived there first. The only difference between the two implementations is that one calls
`M %*% X` on a block and the other on a column, which is a different BLAS path and a
different last bit.

The two run bitwise identically against *themselves* at any block width, which is what the
lockstep test asserts and what it is for. The reference test now asserts agreement at four
steps, where the gap is 1e-15, and says in the file why it cannot be carried further.

**This is the sharpest argument yet for the outer loop measuring.** CG's recurrence residual
drifts; MINRES's is in a different norm; LSQR's mid-flight state is not reproducible to more
than a few digits *across implementations of the same algorithm*. Nothing the recurrence
believes about itself can be reported.

## 4. Finite termination does not survive floating point

In exact arithmetic the Krylov space of a 15-column operator is complete at step 15 and the
iterate there is the answer. Measured against the closed-form least-squares solution:

|  | seed 1 | seed 2 | seed 3 | seed 4 | seed 5 |
|---|---|---|---|---|---|
| step 14 | 5.4e-02 | 2.2e-02 | 9.5e-02 | 1.5e-01 | 1.4e-02 |
| step 15 | 3.8e-02 | 2.2e-02 | 9.5e-02 | 1.5e-01 | 1.3e-02 |
| step 16 | 1.7e-03 | 1.5e-03 | 9.8e-03 | 3.1e-03 | 1.7e-04 |
| step 18 | 1.9e-05 | 1.4e-04 | 7.2e-06 | 8.3e-06 | 7.0e-05 |
| step 20 | 2.4e-09 | 3.4e-09 | 6.6e-10 | 5.2e-09 | 4.5e-10 |
| step 25 | 5.3e-15 | 2.3e-15 | 5.5e-15 | 4.9e-15 | 1.6e-15 |

At the step where the method is supposed to be finished it is wrong in the second digit. It
needs about 25, and this implementation reaches 1e-12 in 22 or 23. So the iteration budget
defaults to a multiple of `n` rather than to `n`, and a test that asserted termination at `n`
would be asserting an exact-arithmetic theorem against floating-point code.

## 5. The least-squares preconditioner rows are narrowed, and the contract had a gap

Plan section 4.3 marks the LSQR and LSMR rows unrestricted, with the instruction that a row
stays open until a check narrows it. The check narrows both to `right`.

The iteration runs on `A M^-1` and recovers `x = M^-1 y`, so `M^-1` is applied to vectors of
the *domain* and `M` is `n x n`. A left preconditioner acts on the residual, which lives in
the codomain, and for a rectangular `A` that is not the same space: an `m x m` map is not
conformable with an `n x n` one. Where the two spaces do coincide it is still the wrong
object, because `min ||M^-1 (b - A x)||` is a weighted least-squares problem with a different
minimiser, so the side would be changing the answer rather than the path to it. Split has the
same defect, since `L^-1 A L^-H` puts `L^-1` on the codomain again.

Two rows now carry a `reason` string, and `check_preconditioner()` reads it from the table
rather than testing for FGMRES by name.

**The gap.** The adjoint of `A M^-1` is `M^-H A^H`, and the `preconditioner` object carried
`M^-1` only. Every other method in v0.1 applies `M^-1` alone, so nothing had needed the other
direction. `preconditioner()` takes an optional `apply_inverse_adjoint`; a preconditioner
declaring `hermitian = TRUE` supplies it for nothing, since `M^-H = M^-1` there, and one
declaring neither is refused by name rather than silently assumed. Both construction routes
carry it through: `as_preconditioner.linop()` builds it from `adjoint(x)`, and
`as_preconditioner.linsolve()` reads the `adjoint` field the `linsolve` contract already
had, which is the second time that field has turned out to be there for a reason.

The declaration is what makes the hermitian shortcut legitimate. Inferring `M^-H = M^-1` from
a preconditioner that happens to be symmetric in the cases tested is the laundering pattern
section 5.3 exists to prevent.

## 6. What carried over unchanged

| Decision | Held? |
|---|---|
| Several right-hand sides in lockstep, bitwise equal to per-column runs | yes, including a block whose columns stop on different tests at different steps |
| The outer loop measures, the inner loop trusts | yes, and see section 3 |
| Non-convergence returns a certificate, contradiction throws | yes, and there is nothing to contradict: like GMRES, LSQR requires no capability |
| `precond_applier()` as the single guard | extended rather than copied: `guarded()` now serves both directions |
| The condition limit at `1/eps` | yes. `GMRES_CONDITION_LIMIT` is now `KRYLOV_CONDITION_LIMIT` in `solvers-common.R`, reached by a second route: the singular values of the Golub-Kahan bidiagonal lie inside those of `A`, so its condition is a lower bound on `cond(A)` and the test stops later, never earlier |

Three copies of the solver preamble had accumulated, differing only in the method name in
their messages. `solver_setup()` in `solvers-common.R` replaces them and takes `square` as a
parameter, which is the only thing a rectangular method needed of it. The `x0` message now
names the shape the iterate has to have rather than the shape of the right-hand side, because
for a rectangular operator those are different.

## 7. Arithmetic: the third case

| Method | Krylov scalars | why |
|---|---|---|
| CG, MINRES | real | `A` is hermitian, so `<x, A y>` is real; `col_dot()` is `Re(<x, y>)` |
| GMRES, FGMRES | complex | Arnoldi's `<v_i, A v_j>` is genuinely complex; `col_cdot()` |
| LSQR | real | `alpha` and `beta` are norms *by construction*, so the bidiagonal is real and non-negative for a complex `A` as much as for a real one |

LSQR needs no inner product at all. Every scalar in the recurrence, the rotations included,
is a norm or is built from norms, so `col_norms()` is the whole requirement and neither
`col_dot()` nor `col_cdot()` appears in the file. That is worth stating because it looks like
an oversight otherwise.

## 8. What is not here

`damp`, the Tikhonov term of the published algorithm, solving
`min ||b - A x||^2 + damp^2 ||x||^2`. It is a small addition to the recurrence and a large one
to the certificate: the quantity that goes to zero becomes `A^H r - damp^2 x`, so Stewart's
perturbation has to be rederived for the damped normal equations before any line above means
what it says. Deferred deliberately, not overlooked.

## 9. Against the plan

| Plan says | Measured |
|---|---|
| 4.3, LSQR and LSMR rows unrestricted by side | Narrowed to `right`. Left and split act on the codomain, which for a rectangular operator is the wrong space, and where the spaces coincide they change the minimiser |
| 4.3, the preconditioner contract is `apply_inverse` | Insufficient. A method running on `A M^-1` needs `M^-H`, and `preconditioner()` had no way to express it |
| 6, the certificate table | No least-squares reading anywhere. The residual and backward-error lines as written report a converged least-squares solve as a failure |
| 7.1, "MINRES and LSMR are the two with real numerical risk" | LSQR carries the same risk. Its recurrence is not reproducible across implementations past about ten steps, and its finite termination is gone by step `n` |
| 10, closed-form fixtures | Two added, both rectangular: `lsq_prescribed()` with a dialled-in singular value spectrum, and `diff_1d()`, matrix free and rank deficient, whose minimum-norm solution is closed form |
