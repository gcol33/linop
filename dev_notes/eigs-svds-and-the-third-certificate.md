# eigs(), svds(), and the certificate that is not a solve certificate

Written 2026-07-28. Section 7.2, the reference eigensolver and the reference singular value
solver. Six findings, three of which correct or extend the plan.

## 1. The eigenpair certificate is a third shape, and its forward line is the first bound

`solve_certificate()` has two readings, and both are backward errors on a linear system:
`b` is data, so a perturbation may touch `A` and `b`, and Rigal-Gaches divides by
`||A|| ||x|| + ||b||`. An eigenpair is not that. `theta` and `x` are both outputs of the
solve, so `A` is the only datum and the denominator has no second term. Trying to reach that
by passing a flag into `solve_certificate()` would have meant a denominator that switches
meaning on the flag, which is the shape section 6 is written against.

What makes it a certificate rather than a report is that the perturbation is exhibited, in
the way Stewart's is for LSQR. A Ritz value is the Rayleigh quotient of its own Ritz vector,
so for a unit `x` the residual `r = A x - theta x` satisfies `x^H r = 0`, and

```
E = -r x^H - x r^H
```

is hermitian with `E x = -r`, hence `(A + E) x = theta x` exactly. In the orthonormal basis
`{x, r/||r||}` it is `[[0, -||r||], [-||r||, 0]]`, so `||E||_2 = ||r||`. The pair is exact for
an operator `||r||` away from `A` **and in the same class as `A`**, which is what the next
step needs.

That step is new to the package. Weyl bounds the movement of a hermitian spectrum by the
norm of a hermitian perturbation, and `theta` is an exact eigenvalue of `A + E`, so

```
min_j |theta - lambda_j(A)| <= ||r|| / ||x||.
```

Every other `forward error` row in the package is `not_checked`, because the bound for a
linear solve needs `||A^-1||` and no residual implies one. Here there is a real bound, so
this is the first row in the package carrying `guarantee = "deterministic_bound"`, and the
first with `source = "theorem"`.

Two consequences.

**The summary line had a latent bug that only a deterministic bound could expose.**
`build_certificate()` computed its `without_deterministic_bound` list as
`guarantee != "identity"`, which was indistinguishable from correct while `identity` and
`estimate` were the only two guarantees any row used. With this row present it would list a
theorem's bound among the checks that have no deterministic bound behind them. The test is
now membership in `identity` or `deterministic_bound`. It changes no existing certificate,
because no existing row used the other three values.

**The bound is only as good as the hermitian declaration, and says so.** `cert_rows()` now
accepts an `evidence()` object instead of the three flat fields, and the forward row passes
`evidence("theorem", "deterministic_bound", 1, depends_on = list(<the hermitian evidence>))`.
So `evidence_satisfies()` sees the declaration underneath the row exactly as it already sees
it underneath a propagated capability. A dense symmetric leaf establishes symmetry by an
exact check on data it already holds and clears a requirement that excludes declarations; the
same operator behind a callback with `properties = c(hermitian = TRUE)` does not. That is
section 5.3's laundering case reaching the certificate, and `test-eigs.R` asserts both halves.

## 2. `eigs()` applies no evidence minimum, deliberately

Every other dispatcher in the package that applies one has somewhere to fall back to: `auto`
refuses CG and takes MINRES, refuses MINRES and takes GMRES. `eigs()` has one method and no
non-hermitian fallback until Arnoldi arrives with the backend, so a minimum would not select
between methods, it would only refuse to run -- and it would refuse on exactly the operators
the package exists for, since a `fun` leaf can only ever declare its own symmetry.

So the value gates the run and the evidence is reported. `require_capability()` needs
`hermitian` to be `TRUE`; what established it goes into the forward-error row and is
visible to any requirement a caller cares to apply. This is the first place in the package
where "carry the evidence" is the answer rather than "apply the minimum", and it is worth
saying that the alternative was not available rather than not preferred.

## 3. The svds certificate is the eigs certificate, on an operator that is hermitian for free

The augmented operator

```
H = [ 0    A  ]
    [ A^H  0  ]
```

is hermitian for every `A`, and its eigenpairs are `(sigma_i, [u_i; v_i]/sqrt(2))` and
`(-sigma_i, [u_i; -v_i]/sqrt(2))`. So a singular triplet **is** a hermitian eigenpair,
`svds()` builds `H` and calls `eigen_certificate()` on it, and every line means there what it
means for `eigs()`. `[u; v]/sqrt(2)` is already unit, so nothing is rescaled, and the two
applies the residual costs are one in each mode, which is what a triplet residual needs
anyway.

The asymmetry that falls out is the interesting part. `H` is hermitian by *construction* with
an empty `depends_on`, so the forward bound on a singular value rests on nothing anyone
declared. The same bound out of `eigs()` rests on whatever established the operator's own
symmetry, which for a callback leaf is a bare declaration. A singular value therefore carries
a stronger certificate than an eigenvalue of the same operator, and it costs nothing to get.

`||H||_2 = ||A||_2` exactly, since the spectrum of `H` is the singular values with both
signs, so the norm estimate is inherited structurally rather than re-estimated by a power
iteration on `H`. That is a structural rule over a possibly-estimated child, so it records
`construction <- [whatever norm2() found]` and the estimate stays visible at the top.

One line is genuinely weaker in the augmented reading and it is worth writing down: the
`orthogonality` row measures `Z^H Z - I` for `Z = [U; V]/sqrt(2)`, where a deviation in `U`
could in principle be cancelled by one in `V`. `test-svds.R` asserts `U` and `V` orthonormal
separately instead.

## 4. Thick restarting is not a refinement, it is what makes the method converge

The first draft restarted from the sum of the unconverged Ritz vectors, which is the obvious
thing and throws away the subspace that produced them. Measured with only that block swapped
(`dev_notes/spikes/restart-comparison.R`), everything else being the shipped code:

| fixture | restart | rounds | iterations | worst backward error | converged |
|---|---|---|---|---|---|
| `laplacian_1d(60)`, 4 smallest algebraic, `ncv = 24` | simple | 9 | 240 | 5.68e-4 | 0 of 4 |
| | thick | 6 | 96 | 5.06e-13 | 4 of 4 |
| `laplacian_1d(60)`, 4 largest algebraic, `ncv = 24` | simple | 6 | 168 | 9.55e-3 | 0 of 4 |
| | thick | 6 | 94 | 4.20e-11 | 4 of 4 |
| `spd_prescribed(40)`, 3 largest, `ncv = 12` | simple | 0 | 12 | 1.745e-15 | 3 of 3 |
| | thick | 0 | 12 | 1.745e-15 | 3 of 3 |

The simple restart does not converge slowly, it stalls: nine rounds and four fifths of the
budget to reach a backward error of 5.7e-4 and then stop improving, because each round has to
rediscover from one direction what the last one spent 24 applies learning.

The third row is the control the package's own rules ask for. Where the request fits in one
round there is no restart to differ over, and the two agree to the last digit -- so the
difference is confined to the block that was swapped, and it is decisive exactly where a
restart happens.

Wu and Simon's thick restart keeps the wanted Ritz vectors themselves as the start of the
next basis. They are already invariant to the accuracy they have reached, so the projected
problem keeps them exactly and one row and column of coupling ties them to the new
directions:

```
[ theta_1                b_1 ]
[        ...             ... ]        b_i = beta_last * s_i[last]
[            theta_p     b_p ]
[ b_1    ...     b_p   alpha ]  ->  tridiagonal from here on
```

For `svds()` the same restart is the broken arrow of the augmented form, with
`rho_i = beta_last * s_u,i[last]` from `A^H U = V B^H + beta_last v_next e^H` and `v_next`
orthogonal to `V`.

## 5. The Rayleigh quotient is reported, not the Ritz value, and that is what makes shift-invert honest

`sigma + 1/theta` inherits whatever the inner solve got wrong. `x^H A x / (x^H x)` is the
minimiser of `||A x - mu x||` over `mu`, so it is the value for which the exhibited
perturbation is smallest, and it is measured on `A` rather than on the transformed operator.
For a plain run the two agree to rounding, since a Ritz value of a hermitian operator is the
Rayleigh quotient of its own Ritz vector; under a shift they do not.

This is the outer-loop-measures discipline reaching the eigenproblem, and it is what lets the
shift-invert route be built out of the package's own solvers -- `A - sigma I` through
`solve(..., method = "minres")` -- without anything downstream having to trust the inner
tolerance. `test-eigs.R` checks that the shifted run finds the three eigenvalues nearest
`sigma = 2` on the 60-point Laplacian to 1e-9 against the closed form, and that an unshifted
run at the same budget finds different ones, so the knob is not a relabelling.

## 6. Storage is the whole of what "reference rather than production" means

Neither method restarts implicitly. A round holds at most `ncv` vectors and orthogonalises
each new direction against all of them, which is O(ncv) storage and O(ncv^2 n) work in the
orthogonalisation. That is the gap `linop.primme` closes in Phase 3, and it is a statement
about cost rather than about accuracy: the certificate is the same object every solver in the
package returns, and rests on the same measurements.

Full reorthogonalisation here is block classical Gram-Schmidt with a second pass under the
Daniel-Gragg-Kaufman-Stewart criterion, which is the criterion GMRES already used; the
constant moved to `solvers-common.R` as `REORTH_ETA` rather than being written twice. The
loop itself is not shared with GMRES and should not be: GMRES needs the coefficients, because
they are the entries of the Hessenberg it minimises over, and here they are a correction and
are discarded. Keeping them would turn a Lanczos recurrence into an Arnoldi one.

The same distinction separates `svds()` from `solvers-bidiag.R`, which builds the same
Golub-Kahan recurrence. LSQR and LSMR store no basis, which is why both lose orthogonality
after a handful of steps and why `lsqr-and-the-least-squares-certificate.md` can only assert
four. A singular triplet is read off the basis, so here the basis has to be kept and kept
orthogonal, and the reorthogonalisation a least-squares solve has no use for is the
load-bearing part. Sharing `bidiag_step()` would mean putting storage those methods do not
want into their inner loop.

## What is not here

- **The generalized problem `A x = lambda B x`.** Lanczos in the `B` inner product is a
  different recurrence rather than a flag on this one: every inner product changes and a
  solve against `B` enters each step. `eigs()` takes `B` and refuses it by name.
- **RSpectra delegation**, which plan section 7.2 puts in v0.1. Its interface is
  `function(x, args)`, single vector and real, so it can carry neither a complex operator nor
  a block apply, and the reference implementations have to exist regardless. Wiring it in
  before the section 3 backend registry exists would mean an ad-hoc branch the registry then
  has to absorb. Deferred to Phase 3 with the registry, which is a departure from the plan
  and is recorded here as one.
- **Arnoldi**, so a non-hermitian operator has no eigensolver. `eigs()` refuses by name.
- **A result padded to `k`.** A run whose subspace never reached `k` wide returns what it
  found; the convergence line counts against what was requested, so a short result still
  reads as incomplete.
