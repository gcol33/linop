# MINRES, and the norm the recurrence is not measured in

Written 2026-07-27. The second of the seven Krylov methods. Three things the CG template
did not anticipate, one of which was a bug in the first draft and was caught by measurement.

## 1. What carried over from CG unchanged

| Decision | Held? |
|---|---|
| Several right-hand sides in lockstep, one block apply per step | yes, bitwise per column |
| Converged columns leave the active set | yes |
| Non-convergence returns a `fail` certificate, it is not thrown | yes |
| A contradicted declaration throws, naming the capability | yes, but on a different kind of test; see 3 |
| The outer loop measures, the inner loop trusts | yes, but the two now measure different quantities; see 2 |

The recurrence itself needed no correction. The iterate at step `m` matched a dense
least-squares solve over the same Krylov space to 1e-16 relative on the first run, at every
`m` from 1 to 12, so the rotation algebra and the direction-vector update were right as
derived. `test-solve-minres.R` keeps that check, because it is an identity about every
intermediate iterate rather than a tolerance statement about the answer.

The cost that is genuinely different is storage. CG carries three persistent `n x k` blocks;
MINRES carries eight. That is inherent to the method, not to this implementation.

## 2. The preconditioned residual is in a different norm, and the stopping test has to know

For a hermitian positive definite `M = L L^H`, MINRES runs on `L^-1 A L^-H` and

```
|| L^-1 (b - A x) ||_2^2  =  (b - A x)^H M^-1 (b - A x),
```

so the quantity it minimises is `||r||_{M^-1}`, not `||r||_2`. The scalar the rotations
carry, `phibar`, is that quantity. Unpreconditioned the two coincide, and a probe confirms
`phibar` equals the recomputed `||b - A x||_2` to 1.4e-16 relative at every step.

This is where MINRES departs from CG. `dev_notes/cg-and-the-arithmetic-floor.md` records
that CG's three preconditioner sides generate the same iterates and that CG reports the true
`b - A x` regardless, so nothing downstream depends on the preconditioner. Neither half of
that is true here. A solver that reported `phibar` as its residual would claim a tolerance
it had not reached, by a factor that grows with the conditioning of `M`.

The fix keeps the CG architecture and converts the currency at the boundary. Each outer
round measures `||r||_2`, converts the caller's target into the `M^-1` norm by the ratio the
two have on that round's own residual, runs the recurrence against the converted target, and
then re-measures in the norm the caller asked about. The conversion is exact for the residual
it is computed from and approximate thereafter, so the outer loop is what makes it safe
rather than what makes it right.

Measured against a deliberately bad diagonal `M`:

| spread of `M` | iterations | restarts | true relative residual, asked 1e-12 |
|---|---|---|---|
| 1e0 | 65 | 0 | 4.6e-13 |
| 1e3 | 232 | 1 | 7.2e-13 |
| 1e6 | 1491 | 3 | 9.4e-13 |

The restarts are the conversion being wrong by enough to matter and the outer loop catching
it. Every case still lands under the tolerance that was requested.

Restarting costs more here than in CG, because MINRES minimises over the Krylov space and a
restart discards it. That is the reason the inner loop returns only when it believes the
target is met, rather than on a fixed inner budget.

## 3. The obvious symmetry test is a test of orthogonality, and it fails on hermitian operators

CG contradicts a false `positive_definite` on a sign: `p^H A p <= 0` cannot happen for a
definite operator. MINRES needs the analogue for `hermitian`, and non-hermitian is not a
one-sided condition, so it has to be a threshold. Which quantity to threshold turned out to
matter more than where to put the threshold.

The first draft used the identity the recurrence already carries,
`v_{j-1}^H A v_j = beta_j`, which follows from taking the `M` inner product of the recurrence
at step `j-1` with `v_j`. It is free. It is also wrong for this purpose: it holds only while
the `v_j` stay orthogonal, and classical Lanczos loses orthogonality as Ritz values converge.
It tests orthogonality, not symmetry.

Worst relative violation, over hermitian operators only:

| fixture | recurrence identity | adjoint identity |
|---|---|---|
| indefinite spread, kappa 1e5 | 8.8e-09 | 1.2e-14 |
| clustered at +-1 | 4.5e-09 | 1.1e-16 |
| tight cluster with a gap | **9.99e-01** | 5.3e-10 |
| near singular | 3.6e-12 | 4.1e-16 |
| complex hermitian, 5 seeds | 5e-02 to 9.9e-01 | 7.1e-11 |

Against a signal of 1.3e-06 at a skew part of relative size 1e-6, the recurrence identity has
no usable threshold at all: its noise on correct operators exceeds its signal on wrong ones
by six orders of magnitude. In the first draft this refused 8 of 40 hermitian fixtures.

The identity that works is the definition, `<x, A y> = <A x, y>`, on the two vectors the
iteration already holds: `A v_j` is this step's apply and `A v_{j-1}` was the previous one,
so it costs no apply either. It does not involve orthogonality, and it stays at 5.3e-10 in the
worst case across the same fixtures. `MINRES_SYMMETRY_TOL` is 1e-6, three orders above that
and at the signal for a 1e-6 skew part. Zero false alarms over 40 runs after the change.

The imaginary part of `v^H A v` is a third signal, independent and equally free, and vacuous
in real arithmetic.

**Below a skew part of 1e-6 the contradiction test is the wrong instrument and does not
fire.** What happens instead is the solve fails to converge and the certificate reports the
true residual it reached. That division is deliberate: the throw is for an operator that is
grossly not what it declared, and the certificate is for everything else.

`test-solve-minres.R` asserts both halves of this, so restoring the cheaper identity fails a
test rather than quietly refusing correct operators.

## 4. A scalar preconditioner cannot demonstrate preconditioning

Recorded because the first draft of the test suite contained it and it passed.

`M = c I` leaves the Krylov space unchanged, so preconditioned and plain MINRES produce
**bitwise identical** iterates. A test comparing them shows nothing. The first draft appeared
to pass only because the two sides were given different iteration budgets, which made it a
comparison between budgets.

The replacement is a badly scaled indefinite operator at kappa 3.2e7 with both sides on the
same budget: plain exhausts 600 iterations at a relative residual of 2.7e-02, Jacobi on
`|diag(M)|` converges in 111 to 4.5e-13. The absolute value is not decoration; a
preconditioner for MINRES has to stay positive definite where the operator does not, or the
`M^-1` norm it minimises is not a norm.

## 5. The MINRES side row, and why it was not narrowed

`dev_notes/fgmres-and-preconditioner-sides.md` section 3 leaves six rows unrestricted pending
a sourced pass over Saad chapter 9, and reserves narrowing for that pass. Implementing the
method produced a derivation that bears on the MINRES row, so it is recorded here rather than
acted on.

The derivation in `R/solvers-minres.R` is the split form, `M = L L^H`. Left preconditioning
in the `M` inner product gives the same iterates, by the same argument that covers CG. Right
preconditioning does not: MINRES needs a hermitian operator to run a three-term recurrence at
all, and `A M^-1` is not hermitian for a general hermitian positive definite `M`. So the
MINRES row looks narrower than the table currently states.

That is a derivation, not the citation the open item asks for, and the governing note is
explicit that narrowing a row is a deliberate edit against the book. **The row is unchanged
and the test asserting all three sides still passes.** This is the strongest candidate for the
sourced pass when it happens, and it is the one place where the conservative reading may be
doing more than failing to refuse.

## Corrections to the plan

| Plan location | Correction |
|---|---|
| 7.1, "MINRES and LSMR are the two with real numerical risk ... correct recurrence, breakdown detection, condition estimation and stopping tests" | The recurrence and the breakdown paths were right as derived and needed no correction. The risk landed on the two the plan does not name: the stopping test under a preconditioner, because the recurrence measures a different norm than the certificate reports, and the contradiction test, which is a linop-specific consequence of declaring `hermitian` and has no counterpart in a solver that does not carry evidence. |
| 6, certificate table | Unchanged. Every line behaves as it does for CG, including the arithmetic floor: at `tol = 0` a converged solve certifies `qualified` with the floor and `fail` without it, on byte-identical iterates. |
