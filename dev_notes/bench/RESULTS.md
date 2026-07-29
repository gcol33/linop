# Benchmark results

Produced by `dev_notes/bench/run-bench.R`, which regenerates this file and
everything under `results/`. Gate 2's last line asks for a harness that runs
end to end with committed results, and these are they.

Two units. Operator applies is the primary one: it is what a matrix-free
caller pays and it is the same number on every machine. Wall time is
secondary, because it is a fact about the machine below and about no other.
Applies also do not price everything, which the GMRES and eigensolver rows
make visible: both orthogonalise against a stored basis, and that arithmetic
is invisible to a count of applies.

```
date                 2026-07-29
R                    R version 4.6.0 (2026-04-24 ucrt)
platform             x86_64-w64-mingw32
os                   Windows 10 x64
cpu                  Intel64 Family 6 Model 183 Stepping 1, GenuineIntel
BLAS                 no path reported by extSoftVersion()
LAPACK               no path reported by La_library()
LAPACK version       3.12.1
timing repetitions   3
harness runtime      2.7 minutes
```

## Seven methods, eight operators

Same tolerance, same budget, same right-hand side and the same zero start for every method on an operator. The roster per fixture is what the operator supplies rather than a ranking: `convdiff_1d` is nonsymmetric and offers the forward action alone, which are two separate exclusions, CG and MINRES for the missing symmetry and LSQR and LSMR for the missing adjoint. Those two are the only methods that require an adjoint. CG is also absent wherever the operator is indefinite. The `norm2()` row under each fixture is the operator-norm estimate the certificate's arithmetic floor rests on, priced on its own and already included in every solve row above it.

Forward error is measured against the closed form in the `truth` column of `results/solvers.csv` and is not something the certificate claims: `forward error` reports `not_checked` in every one of these runs, for the reason in section 6.1.

| fixture | kappa | method | iterations | applies | seconds | backward_error | forward_error | overall |
|---|---|---|---|---|---|---|---|---|
| laplacian_1d(200) | 1.64e+04 | cg | 200 | 213 | 0.0113 | 3e-16 | 1.92e-15 | qualified |
| laplacian_1d(200) | 1.64e+04 | minres | 200 | 213 | 0.0209 | 8.33e-15 | 9.96e-14 | qualified |
| laplacian_1d(200) | 1.64e+04 | gmres | 2000 | 2079 | 0.569 | 4.63e-08 | 0.000723 | fail |
| laplacian_1d(200) | 1.64e+04 | bicgstab | 212 | 436 | 0.0289 | 2.48e-13 | 1.06e-09 | qualified |
| laplacian_1d(200) | 1.64e+04 | lsqr | 1195 | 2407 | 0.119 | 9.61e-12 | 5.74e-10 | qualified |
| laplacian_1d(200) | 1.64e+04 | lsmr | 1204 | 2425 | 0.161 | 1.29e-11 | 2.85e-09 | qualified |
| laplacian_1d(200) | 1.64e+04 | norm2() | - | 10 | 0.000556 | - | - | - |
| shifted_laplacian_1d(200, 0.9) |  495 | minres | 200 | 213 | 0.0205 | 5.32e-14 | 7.52e-14 | qualified |
| shifted_laplacian_1d(200, 0.9) |  495 | gmres | 2000 | 2079 | 0.71 | 0.000757 | 0.254 | fail |
| shifted_laplacian_1d(200, 0.9) |  495 | bicgstab | 386 | 784 | 0.05 | 1.36e-10 | 2.59e-08 | qualified |
| shifted_laplacian_1d(200, 0.9) |  495 | lsqr | 356 | 729 | 0.0449 | 2.43e-10 | 7.14e-10 | qualified |
| shifted_laplacian_1d(200, 0.9) |  495 | lsmr | 357 | 731 | 0.0425 | 1.53e-10 | 6.15e-10 | qualified |
| shifted_laplacian_1d(200, 0.9) |  495 | norm2() | - | 10 | 0.000516 | - | - | - |
| convdiff_1d(80, 0.1) |  900 | gmres | 355 | 370 | 0.109 | 1.93e-10 | 2.42e-08 | qualified |
| convdiff_1d(80, 0.1) |  900 | bicgstab | 80 | 163 | 0.00986 | 1.25e-14 | 6.97e-13 | qualified |
| convdiff_1d(80, 0.1) |  900 | norm2() | - | 1 | 8.8e-05 | - | - | - |
| kms(300, 0.7) | 32.1 | cg | 54 | 67 | 0.00798 | 3.36e-10 | 3.25e-09 | qualified |
| kms(300, 0.7) | 32.1 | minres | 53 | 66 | 0.0157 | 3.61e-10 | 7.31e-09 | qualified |
| kms(300, 0.7) | 32.1 | gmres | 54 | 68 | 0.0239 | 3.84e-10 | 7.87e-09 | qualified |
| kms(300, 0.7) | 32.1 | bicgstab | 41 | 95 | 0.0142 | 1.77e-10 | 4.78e-09 | qualified |
| kms(300, 0.7) | 32.1 | lsqr | 217 | 451 | 0.0958 | 4.04e-10 | 8.71e-09 | qualified |
| kms(300, 0.7) | 32.1 | lsmr | 218 | 453 | 0.085 | 4.85e-10 | 1.35e-08 | qualified |
| kms(300, 0.7) | 32.1 | norm2() | - | 10 | 0.00102 | - | - | - |
| spd_prescribed(300, kappa 1e4) | 10000 | cg | 685 | 700 | 0.0938 | 4.11e-12 | 4.9e-09 | qualified |
| spd_prescribed(300, kappa 1e4) | 10000 | minres | 654 | 669 | 0.117 | 4.36e-12 | 2.42e-08 | qualified |
| spd_prescribed(300, kappa 1e4) | 10000 | gmres | 2497 | 2595 | 1.04 | 4.38e-12 | 4.12e-08 | qualified |
| spd_prescribed(300, kappa 1e4) | 10000 | bicgstab | 573 | 1160 | 0.141 | 2.31e-12 | 1.62e-08 | qualified |
| spd_prescribed(300, kappa 1e4) | 10000 | lsqr | 3000 | 6019 | 1.04 | 0.0184 | 0.303 | fail |
| spd_prescribed(300, kappa 1e4) | 10000 | lsmr | 3000 | 6019 | 1.02 | 0.000159 | 0.423 | fail |
| spd_prescribed(300, kappa 1e4) | 10000 | norm2() | - | 12 | 0.00125 | - | - | - |
| hpd_prescribed(200, kappa 1e3) | 1000 | cg | 228 | 243 | 0.03 | 3.74e-11 | 6.93e-09 | qualified |
| hpd_prescribed(200, kappa 1e3) | 1000 | minres | 224 | 239 | 0.038 | 3.71e-11 | 2.11e-08 | qualified |
| hpd_prescribed(200, kappa 1e3) | 1000 | gmres | 378 | 405 | 0.15 | 3.92e-11 | 3.33e-08 | qualified |
| hpd_prescribed(200, kappa 1e3) | 1000 | bicgstab | 171 | 357 | 0.0496 | 3.44e-11 | 2.5e-08 | qualified |
| hpd_prescribed(200, kappa 1e3) | 1000 | lsqr | 2000 | 4019 | 0.56 | 0.0302 | 2.05e-05 | fail |
| hpd_prescribed(200, kappa 1e3) | 1000 | lsmr | 2000 | 4019 | 0.55 | 0.00134 | 3.28e-05 | fail |
| hpd_prescribed(200, kappa 1e3) | 1000 | norm2() | - | 12 | 0.00197 | - | - | - |
| lsq_prescribed(600 x 200, kappa 1e3) | 1000 | lsqr | 2000 | 4019 | 0.816 | 1.07e-07 | 1.44e-05 | fail |
| lsq_prescribed(600 x 200, kappa 1e3) | 1000 | lsmr | 1977 | 3973 | 0.815 | 9.95e-09 | 3.07e-05 | qualified |
| lsq_prescribed(600 x 200, kappa 1e3) | 1000 | norm2() | - | 12 | 0.00158 | - | - | - |
| diff_1d(400) |  255 | lsqr | 399 | 815 | 0.0418 | 7.35e-16 | 2.63e-15 | qualified |
| diff_1d(400) |  255 | lsmr | 399 | 815 | 0.0457 | 7.89e-16 | 4.56e-15 | qualified |
| diff_1d(400) |  255 | norm2() | - | 10 | 0.000296 | - | - | - |

## Right-hand sides in lockstep

k columns through one solve against k columns through k solves. The block takes one apply per step where the separate route takes k, and pays for it in column-applies, because every column keeps being carried until the last one converges. The final column is how far apart the two answers are.

| fixture | k | route | applies | column_applies | seconds | agreement_vs_lockstep |
|---|---|---|---|---|---|---|
| laplacian_1d(200), matrix free | 1 | lockstep | 213 | 213 | 0.0121 | - |
| laplacian_1d(200), matrix free | 1 | one column at a time | 213 | 213 | 0.00985 | 0 |
| laplacian_1d(200), matrix free | 2 | lockstep | 213 | 416 | 0.0152 | - |
| laplacian_1d(200), matrix free | 2 | one column at a time | 426 | 426 | 0.0251 | 0 |
| laplacian_1d(200), matrix free | 4 | lockstep | 213 | 822 | 0.0134 | - |
| laplacian_1d(200), matrix free | 4 | one column at a time | 852 | 852 | 0.0496 | 0 |
| laplacian_1d(200), matrix free | 8 | lockstep | 213 | 1634 | 0.0201 | - |
| laplacian_1d(200), matrix free | 8 | one column at a time | 1704 | 1704 | 0.0922 | 0 |
| laplacian_1d(200), matrix free | 16 | lockstep | 213 | 3258 | 0.0319 | - |
| laplacian_1d(200), matrix free | 16 | one column at a time | 3408 | 3408 | 0.216 | 0 |
| spd_prescribed(300, kappa 1e4), stored | 1 | lockstep | 698 | 698 | 0.0739 | - |
| spd_prescribed(300, kappa 1e4), stored | 1 | one column at a time | 698 | 698 | 0.0729 | 0 |
| spd_prescribed(300, kappa 1e4), stored | 2 | lockstep | 698 | 1376 | 0.0946 | - |
| spd_prescribed(300, kappa 1e4), stored | 2 | one column at a time | 1388 | 1388 | 0.143 | 0 |
| spd_prescribed(300, kappa 1e4), stored | 4 | lockstep | 698 | 2739 | 0.135 | - |
| spd_prescribed(300, kappa 1e4), stored | 4 | one column at a time | 2775 | 2775 | 0.293 | 0 |
| spd_prescribed(300, kappa 1e4), stored | 8 | lockstep | 708 | 5505 | 0.268 | - |
| spd_prescribed(300, kappa 1e4), stored | 8 | one column at a time | 5589 | 5589 | 0.576 | 0 |
| spd_prescribed(300, kappa 1e4), stored | 16 | lockstep | 708 | 10975 | 0.422 | - |
| spd_prescribed(300, kappa 1e4), stored | 16 | one column at a time | 11155 | 11155 | 1.14 | 0 |

## The size that stops being a problem

tridiag(-1, 3, -1) solved through its action, through a dense matrix and through a sparse one. The shift keeps the spectrum in (1, 5) at every n, so the iteration count is flat down the ladder and what moves is size. The sparse column is the fastest route at every size here, which is what a banded operator with an O(n) direct solve should look like.

`object_mb` is `object.size()` of what the route holds: the closure and its counter for the operator, the matrix for the other two.

| n | route | seconds | iterations | applies | object_mb | forward_error | note |
|---|---|---|---|---|---|---|---|
| 500 | matrix free, cg | 0.00561 | 20 | 33 | 0.0484 | 6.2e-09 |  |
| 500 | sparse, form and solve | 0.00112 | - | - | 0.0205 | 0 | the reference the forward errors use |
| 500 | dense, form the matrix | 0.000381 | - | - | 1.91 | - |  |
| 500 | dense, LU solve | 0.0112 | - | - | - | 1.28e-16 | the factorisation copies the matrix |
| 1000 | matrix free, cg | 0.00304 | 20 | 33 | 0.0484 | 5.92e-09 |  |
| 1000 | sparse, form and solve | 0.000568 | - | - | 0.0396 | 0 | the reference the forward errors use |
| 1000 | dense, form the matrix | 0.00131 | - | - | 7.63 | - |  |
| 1000 | dense, LU solve | 0.0878 | - | - | - | 1.22e-16 | the factorisation copies the matrix |
| 2000 | matrix free, cg | 0.00499 | 20 | 33 | 0.0484 | 5.94e-09 |  |
| 2000 | sparse, form and solve | 0.000833 | - | - | 0.0777 | 0 | the reference the forward errors use |
| 2000 | dense, form the matrix | 0.00743 | - | - | 30.5 | - |  |
| 2000 | dense, LU solve | 0.719 | - | - | - | 1.19e-16 | the factorisation copies the matrix |
| 4000 | matrix free, cg | 0.00813 | 20 | 33 | 0.0484 | 6.07e-09 |  |
| 4000 | sparse, form and solve | 0.00148 | - | - | 0.154 | 0 | the reference the forward errors use |
| 4000 | dense, form the matrix | 0.0819 | - | - |  122 | - |  |
| 4000 | dense, LU solve | 6.04 | - | - | - | 1.18e-16 | the factorisation copies the matrix |
| 8000 | matrix free, cg | 0.0142 | 20 | 33 | 0.0484 | 6.18e-09 |  |
| 8000 | sparse, form and solve | 0.00253 | - | - | 0.307 | 0 | the reference the forward errors use |
| 8000 | dense, form the matrix | 0.172 | - | - |  488 | - |  |
| 8000 | dense, LU solve | 45.2 | - | - | - | 1.22e-16 | the factorisation copies the matrix |
| 16000 | matrix free, cg | 0.025 | 20 | 33 | 0.0484 | 6.18e-09 |  |
| 16000 | sparse, form and solve | 0.00471 | - | - | 0.612 | 0 | the reference the forward errors use |
| 16000 | dense, form the matrix | - | - | - | 1.95e+03 | - | not run: the matrix is 1.9 GB before the factorisation copies it |
| 50000 | matrix free, cg | 0.0718 | 20 | 33 | 0.0484 | 6.15e-09 |  |
| 50000 | sparse, form and solve | 0.014 | - | - | 1.91 | 0 | the reference the forward errors use |
| 50000 | dense, form the matrix | - | - | - | 1.91e+04 | - | not run: the matrix is 18.6 GB before the factorisation copies it |
| 100000 | matrix free, cg | 0.123 | 20 | 33 | 0.0484 | 6.17e-09 |  |
| 100000 | sparse, form and solve | 0.0263 | - | - | 3.82 | 0 | the reference the forward errors use |
| 100000 | dense, form the matrix | - | - | - | 7.63e+04 | - | not run: the matrix is 74.5 GB before the factorisation copies it |

## eigs and svds

Both verbs against closed-form spectra, at three subspace sizes, with a dense factorisation of the same operator underneath. The dense routes return the whole spectrum where the two verbs return k pairs, so the seconds column compares different deliverables and the value error column compares the same k values.

Read `nconv` before anything else in this table. Both ends of `laplacian_1d(400)` have a relative gap of 4.6e-5, and at ncv = 20 and 40 the run stalls with nothing converged. Raising `maxit` from 1000 to 5000 or loosening `tol` from 1e-10 to 1e-8 changes neither the iteration count nor the value error there, so the subspace is the binding knob and the stall detector is reporting that rather than giving up early. A row with `nconv` below k is returning its best unconverged pairs, and its value error is the size of that.

The shift-invert row counts every inner MINRES step through the same counter, so its applies column is the cost of the transformation and not of the outer recurrence.

| problem | verb | which | route | nconv | applies | seconds | backward_error | value_error |
|---|---|---|---|---|---|---|---|---|
| laplacian_1d(400) | eigs | largest | lanczos, ncv = 20 | 0 | 296 | 0.0475 | 0.00336 | 0.000999 |
| laplacian_1d(400) | eigs | largest | lanczos, ncv = 40 | 0 | 262 | 0.0394 | 0.000908 | 0.000744 |
| laplacian_1d(400) | eigs | largest | lanczos, ncv = 80 | 6 | 499 | 0.11 | 9.21e-11 | 1.11e-16 |
| laplacian_1d(400) | eigen() | largest | dense eigen(), whole spectrum | - | - | 0.00768 | - | 4.44e-16 |
| laplacian_1d(400) | eigs | smallest | lanczos, ncv = 20 | 0 | 175 | 0.0209 | 0.00314 | 6.44 |
| laplacian_1d(400) | eigs | smallest | lanczos, ncv = 40 | 0 | 178 | 0.0253 | 0.00574 |  4.4 |
| laplacian_1d(400) | eigs | smallest | lanczos, ncv = 80 | 6 | 502 | 0.114 | 1.25e-15 | 5.3e-15 |
| laplacian_1d(400) | eigen() | smallest | dense eigen(), whole spectrum | - | - | 0.00767 | - | 1.93e-11 |
| laplacian_1d(100) | eigs | sigma = 1 | lanczos on (A - sigma I)^-1 | 4 | 3029 | 0.437 | 2.02e-12 | 2.3e-16 |
| diff_1d(400) | svds | largest | golub-kahan, ncv = 20 | 0 | 384 | 0.0386 | 0.0027 | 0.00127 |
| diff_1d(400) | svds | largest | golub-kahan, ncv = 40 | 0 | 514 | 0.0639 | 0.000629 | 0.000376 |
| diff_1d(400) | svds | largest | golub-kahan, ncv = 80 | 6 | 986 | 0.18 | 6.98e-11 | 1.22e-14 |
| diff_1d(400) | svd() | largest | dense svd(), whole spectrum | - | - | 0.0269 | - | 4.44e-16 |
| lsq_prescribed(600 x 200) | svds | largest | golub-kahan, ncv = 20 | 6 | 96 | 0.0286 | 8.5e-12 | 1.58e-15 |
| lsq_prescribed(600 x 200) | svds | largest | golub-kahan, ncv = 40 | 6 | 96 | 0.0255 | 1.97e-15 | 1.15e-15 |
| lsq_prescribed(600 x 200) | svds | largest | golub-kahan, ncv = 80 | 6 | 176 | 0.0562 | 1.61e-15 | 1.28e-15 |
| lsq_prescribed(600 x 200) | svd() | largest | dense svd(), whole spectrum | - | - | 0.0112 | - | 9.86e-16 |

