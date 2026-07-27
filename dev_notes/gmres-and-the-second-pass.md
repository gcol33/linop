# GMRES, FGMRES, and the second pass that has to be earned

Written 2026-07-28. The third and fourth of the seven Krylov methods, in one implementation.
Four things the CG and MINRES template did not anticipate, two of which were wrong in the
first draft and were caught by measurement rather than by reasoning.

## 1. What carried over unchanged

| Decision | Held? |
|---|---|
| Several right-hand sides in lockstep, one block apply per step | yes, bitwise per column |
| Converged columns leave the active set | yes, and now mid-round rather than only at a restart |
| Non-convergence returns a `fail` certificate, it is not thrown | yes |
| The outer loop measures, the inner loop trusts | yes, and the currency conversion MINRES introduced is needed on two of the three sides |
| A contradicted declaration throws, naming the capability | only for the preconditioner; GMRES requires nothing of the operator, so it has nothing to contradict |

The last row is the structural difference. CG needs positive definiteness, MINRES needs
hermitian, and GMRES needs a square operator and an apply. It is therefore the method that
answers when the capability set says nothing at all, and the only one in the package that
works on an operator supplying no adjoint: every apply is mode `N`. `norm2()` falls back to
random probes there rather than failing, so the certificate still has a `||A||` to rest on.

The recurrence needed no correction. The iterate at step `m` matched a dense least-squares
solve over the same Krylov space to 1e-16 relative on the first run, at every `m` from 1 to
12, so the Arnoldi loop, the rotation algebra and the back substitution were right as
derived. `test-solve-gmres.R` keeps that check.

## 2. Complex arithmetic stops being optional

`aaa-utils.R` computes the Krylov scalars in real arithmetic, and the comment there says why:
`col_dot()` is `Re(<x, y>)`, which is the whole inner product for the quantities a hermitian
operator's recurrences form. Both hermitian methods rely on that.

Arnoldi has no such guarantee. `<v_i, A v_j>` for a general complex operator is complex, and
so is the sine of the rotation that eliminates `h_{j+1,j}`. Keeping only the real part
orthogonalises against the wrong direction, silently, and the solve simply fails to converge.

`col_cdot()` is the addition: the whole inner product, of which `col_dot()` is the real part
and `col_dot_im()` the imaginary one. The three are kept separate rather than layered,
because the two hermitian methods are written in real arithmetic deliberately and routing
them through a complex product to discard the imaginary half would cost them that.

The rotation is the LAPACK `zlartg` convention, `cs` real and `sn` complex. `g` is a norm and
so never negative, which is what lets `cs` stay real, and the real case falls out of the same
expressions with the phase of `f` reducing to its sign. One rotation serves both storage
modes.

## 3. `side` selects an algorithm, and this is the first method for which it does

CG's row in the section 4.3 table is unrestricted because the three sides *agree*: left,
split and right preconditioned CG generate the same iterates. Reading the GMRES row the same
way would be wrong. Its three sides are three different iterations minimising three different
quantities:

| side | Arnoldi on | minimises | conversion needed |
|---|---|---|---|
| right | `A M^-1` | `\|\|b - A x\|\|_2` | none, this is what the certificate reports |
| left | `M^-1 A` | `\|\|M^-1 (b - A x)\|\|_2` | yes |
| split | `L^-1 A L^-H` | `\|\|b - A x\|\|_{M^-1}` | yes |

The conversion is the one MINRES introduced, unchanged: the outer round measures `||r||_2`,
converts the caller's target by the ratio the two norms have on that round's own residual,
and re-measures in the euclidean norm afterwards.

Split preconditioning looked at first like something the contract could not express. The
`preconditioner` object carries `M^-1` as a single map and never `L`, and `L^-1 A L^-H` needs
both halves. It is reachable by the change of variable MINRES uses, `v_j = L^-H q_j`, which
turns euclidean orthogonality of the `q_j` into `M`-orthogonality of the `v_j`:

```
h_ij = <v_i, A v_j>,   w = A v_j - sum_i h_ij u_i,
h_{j+1,j} = sqrt(<w, M^-1 w>),
u_{j+1} = w / h_{j+1,j},   v_{j+1} = M^-1 w / h_{j+1,j},
```

carrying `u_j = M v_j` alongside `v_j` so that `M` itself is never applied. Two bases instead
of one, and only `M^-1` is ever called. **The three sides were verified to produce genuinely
different iterates**, so the test that asserts all three are accepted is not asserting that
one implementation was relabelled three times.

The split form exists only for a hermitian positive definite `M`, and the GMRES row asks only
for `fixed`. So the check is at run time on `<w, M^-1 w>`, and it produces two different
messages: a preconditioner that declared `positive_definite = TRUE` has been contradicted,
one that declared nothing has not, and reporting them the same way would accuse a caller who
never made the claim.

FGMRES is a flag, as plan 7.1 says. Right preconditioning with a fixed `M` recovers its
update by applying `M^-1` once to the assembled step; with a flexible one it has to store the
preconditioned basis. That is the whole difference, and the two agree to 9.5e-15 on a fixed
preconditioner at the same iteration count. They are not bitwise equal, because
`M^-1 (V_m y)` and `sum_i (M^-1 v_i) y_i` are the same number by different routes.

## 4. The second Gram-Schmidt pass, taken unconditionally, is worse than not taking it

The first draft reorthogonalised on every step. The reasoning was that the usual criterion,
a second pass only where the first one cancelled badly, is data dependent per column and
would cost the lockstep identity that CG and MINRES both hold.

That reasoning is wrong twice over.

It is wrong about the identity: the Daniel-Gragg-Kaufman-Stewart test is computed from a
column's own data, so applying it per column is exactly what a per-column solve would do. The
implementation zeroes the coefficient where the first pass sufficed, and subtracting an exact
zero multiple leaves those columns bit for bit as they were. The lockstep test passes on both
settings of `reorth`.

It is wrong about the arithmetic, which is the part that mattered. Once the first pass has
cancelled a direction to rounding level, the coefficients the second pass computes are noise,
and they are *added to the Hessenberg column* rather than discarded. On a Vandermonde matrix
with `kappa = 3.6e25`, unconditional reorthogonalisation reached a relative residual of
**1.7e3** where plain modified Gram-Schmidt reached **5.8e-1**.

| fixture | m | unconditional second pass | MGS only |
|---|---|---|---|
| Vandermonde n = 40 | 30 | 1.71e+03 | 5.77e-01 |
| Hilbert n = 60 | 40 | 2.28e+00 | 7.75e-01 |

## 5. Which turned out to be a condition estimate missing, not an orthogonalisation bug

Instrumenting the Vandermonde case step by step showed the two settings agreeing exactly
through step 18 and diverging after it. What diverges is not the orthogonality:

- with MGS only, the basis is broken, the iterate stops moving, and `||x||` pins at 1.9e5;
- with reorthogonalisation the basis is *correct*, the space keeps being extended, and `||x||`
  runs to 1.7e10 while the recurrence residual falls monotonically from 0.60 to 0.365 and the
  true residual climbs to 1.7e3.

Better orthogonality, worse iterate. MGS winning there is luck: its stagnation is a symptom
of a broken basis, not a safeguard. The real fault is that a Krylov space goes on admitting
directions long after those directions have stopped meaning anything, and the projected
residual — which is a true statement about the projected problem — cannot tell.

Plan 7.1 names condition estimation among the things a Krylov method has to get right, and
this is where GMRES needs it. `GMRES_CONDITION_LIMIT` stops extending the space once
`max|R_ii| / min|R_ii|` exceeds `1/eps`, retiring the column at the last step still worth
solving. That ratio is a lower bound on `cond(R)`, so it stops later than a sharp estimate
would and never earlier, which is the right direction for a rule that ends an iteration.

With it in place, measured over 6 seeds each:

| fixture | default | `conlim = Inf` | ratio |
|---|---|---|---|
| Vandermonde n = 40 | 5.8e-01 to 6.6e+00 | 1.5e+02 to 9.8e+03 | 22 to 12248 |
| Hilbert n = 60 | 7.9e-01 to 9.3e-01 | 1.4e+00 to 2.7e+00 | 1.7 to 3.1 |
| clustered SPD n = 60 | unchanged | unchanged | 1.0 |
| convection-diffusion, scaled, KMS | **bitwise identical** | | 1.000 |

It never fires on a well-posed problem, which is what a rule that ends an iteration has to be
able to say. `test-solve-gmres.R` asserts both halves.

The certificate reported the truth throughout: every divergent run above certifies as `fail`.
What the condition limit buys is the difference between a bad answer and a much worse one,
not the difference between a wrong certificate and a right one.

## 6. And the accuracy claim about reorthogonalisation is not supported

Having fixed the condition estimate, the obvious next test is that `reorth = TRUE` beats
`reorth = FALSE`. It does not, reliably. Over 12 seeds, drift of the recurrence scalar from
the true residual:

| fixture | median ratio F/T | worst case | seeds where TRUE is worse |
|---|---|---|---|
| convection-diffusion n = 50, scaled 1e6, m = 40 | 1.6 | 0.2 | 3 of 12 |
| convection-diffusion n = 80, scaled 1e8, m = 70 | 1.1 | 0.1 | 5 of 12 |

The first draft of this suite asserted `drift(TRUE) < drift(FALSE)` and passed, on seed 7,
where the ratio happens to be 13.5. That is the MINRES scalar-preconditioner mistake in a new
costume: a knob that plainly does something, asserted to do the thing that reads well, on the
one configuration where it does.

So the test asserts only that the setting changes the iterates, and says in the test body why
it asserts nothing more. `reorth = TRUE` stays the default because it is what makes the basis
orthogonal by construction rather than by luck, and because after the condition fix it is
never *worse* on achieved residual on any fixture tried — but that is a weaker statement than
the one a reader would expect, so it is written down rather than tested.

## 7. A bug the lockstep test found and the recovery tests did not

`shrink()` narrows every per-column vector when a column retires mid-round. `dmax` and `dmin`,
added for the condition estimate, were narrowed in one of the two call sites and not in the
other. The result was a fractional-recycling warning and a condition estimate comparing the
current columns against the widths the round started with.

Nothing in the recovery suite noticed: every solve still converged, because the corrupted
comparison only mis-times a stopping rule. What surfaced it was R's recycling warning under
the multi-column lockstep test, which is the only test that runs the mid-round retirement path
with more than one column in the block. Per-column state now narrows in exactly one place, and
the comment there says that a vector left out of the list is silently recycled.

## Corrections to the plan

| Plan location | Correction |
|---|---|
| 4.3, side table | The GMRES row being unrestricted is a statement that three different algorithms are all available, not, as with CG, that the three sides coincide. Implemented as three code paths, verified to produce different iterates. |
| 7.1, "restarted, MGS with reorthogonalisation" | Unconditional reorthogonalisation is worse than none on an ill-conditioned basis, by three orders of magnitude on the worst fixture measured. The criterion has to be the conditional one, and it is compatible with the lockstep identity, which is why the first draft rejected it. |
| 7.1, "MINRES and LSMR are the two with real numerical risk" | GMRES has a third: a Krylov space that has gone numerically singular, where the projected residual keeps falling while the true residual diverges. The condition estimate the plan lists for LSMR is needed here too, and `GMRES_CONDITION_LIMIT` is the same `1/eps` LSQR and LSMR will want. |
| 6, certificate table | Unchanged. Every line behaves as it does for CG and MINRES, including the arithmetic floor. Worth noting for the next method: on a badly enough scaled operator the floor swallows the tolerance and every run certifies as met, so a test comparing two configurations has to require `pass` rather than accept `qualified`. |

## Open, unchanged

`dev_notes/fgmres-and-preconditioner-sides.md` section 3 leaves six rows unrestricted pending
a sourced pass over Saad chapter 9, and `dev_notes/minres-and-the-preconditioned-norm.md`
section 5 parks a derivation that bears on the MINRES row. Implementing GMRES produced no
reason to narrow the GMRES row: all three sides are implemented and tested. FGMRES remains the
one sourced restriction, against PETSc `KSPFGMRES`. **The rows are unchanged and the tests
asserting them still pass.**
