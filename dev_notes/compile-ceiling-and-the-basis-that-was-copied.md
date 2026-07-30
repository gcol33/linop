# The compile ceiling, and the basis that was copied

The question was whether to move `linop` off zero compiled code. The probes that
answer it are `dev_notes/spikes/eigs-compilable-ceiling.R`,
`dev_notes/spikes/eigs-orth-crossover.R`, `dev_notes/spikes/orth-allocation-share.R`,
`dev_notes/spikes/orth-selftime-attribution.R` and
`dev_notes/spikes/orth-basis-slice.R`, with results beside them.

The answer is no, and the interesting part is that the number which argued for C was
measuring something else. It is not an interpreter floor. It is one missing BLAS
argument, in one routine, reachable from R only by not being R.

## 1. The solvers are settled and the eigensolver is not

After the two allocation fixes (`7eeb4b4`, `bbdfb90`), linop's own R self time in a
solve is 1.1% to 11.0% of the run, so the ceiling on compiling every line of it is
1.01x to 1.12x. Six of the eight verbs have nothing to gain and the matter is closed
for them.

`eigs()` did not follow. `eigs-compilable-ceiling.csv`, at the sizes the package is
for and with two operator costs:

| kind | n | ncv | compiled | linop R | apply | ceiling |
|---|---|---|---|---|---|---|
| stencil | 1e5 | 20 | 44.0% | 27.2% | 28.8% | 1.37x |
| stencil | 1e5 | 80 | 52.7% | 41.2% | 6.2% | 1.70x |
| stencil | 1e6 | 20 | 48.5% | 23.2% | 28.4% | 1.30x |
| stencil | 1e6 | 80 | 53.3% | 36.2% | 10.5% | 1.57x |
| fft | 1e5 | 20 | 75.9% | 17.7% | 6.5% | 1.21x |
| fft | 1e5 | 80 | 65.9% | 32.4% | 1.8% | 1.48x |
| fft | 1e6 | 20 | 82.8% | 15.1% | 2.1% | 1.18x |
| fft | 1e6 | 80 | 72.5% | 26.1% | 1.4% | 1.35x |

Median 1.36x. Two readings of that table were proposed and both were wrong.

**It is not small n and it is not a toy operator.** The grid runs at n = 1e5 and
1e6, against a stencil and an FFT, not against the benchmark harness's
`laplacian_1d(400)`. There was a prediction that the ceiling collapses toward 1.00x
once n is large and the apply is real, on the grounds that the residual R cost is
per-call and therefore fixed per iteration. The table is the refutation: n = 1e5 to
1e6 at ncv = 80 moves the stencil ceiling 1.70x to 1.57x, and the FFT, whose apply
costs n log n, still reads 1.35x.

**Nor is the share diluted by growing n**, and `eigs-orth-crossover.csv` says why in
one column. The fraction of the run inside `orth_against()` at ncv = 80 is 0.673 at
n = 1e4, 0.691 at n = 1e5, 0.671 at n = 1e6. Flat. A fixed per-call cost cannot be
flat in n. Whatever is in there scales with the work.

## 2. The self time is not work

On the largest cell the profile puts `orth_against()`'s own self time at 1.132 s of
3.63 s, the single largest entry, above `%*%` at 1.004 s and `crossprod` at 0.670 s.
Its four lines dispatch every operation to a primitive that the same profile lists
separately, and `-` is 6 ms. Re-running with `gc.profiling = TRUE` breaks out `<GC>`
at 7.3%, so collection is not the bulk either.

Timed rather than sampled, the routine has nothing in it. `orth-allocation-share.R`
compares it against the two products it cannot avoid:

| n | ncv | two products | `orth_against` | removable |
|---|---|---|---|---|
| 1e5 | 20 | 3.17 ms | 3.88 ms | 1.23x |
| 1e5 | 80 | 16.25 ms | 16.29 ms | 1.00x |
| 1e6 | 20 | 43.99 ms | 47.29 ms | 1.07x |
| 1e6 | 80 | 165.20 ms | 171.09 ms | 1.04x |

At the cell where the profile charges it a second, it runs at 1.04x of its own BLAS
floor. And profiled on its own, `orth-selftime-attribution.R` charges it **0.0%**
self, with `%*%` at 48.1% and `crossprod` at 48.1%, matching the wall clock's 98.4%.

So the attribution is correct in isolation and wrong inside the solver, which means
the difference is in how the solver calls it, not in the routine.

## 3. What the difference is

`R/eigs.R:455`, and `R/svds.R:334`, `:356`, `:370`:

```r
W <- orth_against(V[, seq_len(i), drop = FALSE], W)
```

Base R has no matrix view, so `V[, seq_len(i), drop = FALSE]` allocates and copies
`i` columns of the live basis. At n = 1e6 and ncv = 80 that grows to 610 MB per call
and about 26 GB over one round. The copy is one more pass over exactly the memory
each of the two products already streams, so it is a third of the routine rather
than a rounding error.

It lands on `orth_against`'s frame because of lazy evaluation: the promise for `Q`
is forced at the first use, which is `is.null(Q)`, inside the callee. The caller
writes the subscript and the callee is billed for it. Pre-materialising `Q`, as
sections 2's probes do, moves the cost to where it is written and the entry vanishes.

This is the third instance of one family. `7eeb4b4` removed a `Conj()` that
duplicated the basis once per pass; `bbdfb90` removed a `rep()` that built an n x k
array to carry k numbers; this is the same shape and the largest of the three.

## 4. Pure R cannot remove it

There is one candidate that does not change the answer. The unused columns of the
preallocated basis are zero, so their coefficients are zero and orthogonalising
against the whole of `V` returns the same vector as orthogonalising against its
first `i` columns. No copy, at the price of both products running over `ncv`
columns instead of `i`. `orth-basis-slice.R` confirms it is bitwise identical at
every `i` tried, and prices it summed across a round:

| n | ncv | slice (ships) | zeros (pure R) | floor (compiled) |
|---|---|---|---|---|
| 1e5 | 80 | 64.03 ms, 1.00x | 64.03 ms, 1.00x | 34.97 ms, 1.83x |
| 1e6 | 80 | 672.76 ms, 1.00x | 687.36 ms, 0.98x | 376.75 ms, 1.79x |

A wash. The copy is `O(n i)` and the extra BLAS is `O(n (ncv - i))`, and averaged
over a round those are the same quantity, which is why the two columns agree to
within noise. There is no third route: what is wanted is a product against a
submatrix of a matrix that already exists, and expressing that needs a leading
dimension, which base R does not expose.

The `floor` column is what a compiled routine issuing one `dgemm` with an `lda`
argument would take. Carried through to a whole run at the largest cell, where
`orth_against` is 67.1% of it, that is `0.671 / 1.79 + 0.329`, or **1.42x** against
a measured ceiling of 1.57x for the same cell. The ceiling was the right magnitude
for the wrong reason.

## 5. What this decides

Compiled code buys 1.01x to 1.12x on the six solvers, and on the two spectral verbs
it buys about 1.4x at large `n * ncv`, all of it one BLAS call away, none of it
reachable from R.

It does not follow that `linop` should carry `src/`. The 1.4x sits entirely in the
reference eigensolver's full reorthogonalisation, and `R/eigen-common.R:7-10` already
names that loop as the one Phase 3 replaces: a round keeps its whole basis, which is
`O(ncv)` vectors of storage and `O(ncv^2 n)` of work, and is the price of not having
an implicit restart. PRIMME restarts implicitly. It never forms the 610 MB basis and
never runs the orthogonalisation, so it does not fix the copy, it deletes the loop
the copy is in. S0.3 has it building under Rtools45 and on macOS arm64.

So the finding is an argument for bringing `linop.primme` forward, not for compiling
`linop`. Spending the zero-compiled-code claim, which README, CONTRIBUTING, NEWS,
GATE1 and plan section 4 all rest on, to buy 1.4x on the one loop Phase 3 removes is
the worst available trade.

What would change it: a caller who needs large-`ncv` `eigs()` before
`linop.primme` exists. The minimal surface is a single GEMM wrapper taking a
leading dimension, not a rewrite, and the 610 MB basis is the worse problem at that
size anyway. Both are the backend's job.

## 6. Corrections

| Where | Was | Is |
|---|---|---|
| `CLAUDE.md`, settled-decisions table | S0.2's callback figure carries the whole zero-compiled-code decision | It covers the apply path. The orthogonalisation is a separate question with its own measurement |
| `README.md` | 0.13% callback overhead is "why the package carries no compiled code" | That measurement does not reach the eigensolver. The reason is scoped in section 5 above |
| `eigs-compilable-ceiling.csv` | linop R self time, so the ceiling on compiling linop | Inflated by the caller's subscript under lazy evaluation. Read section 2 before quoting a cell |
| Ceiling stated for the package | 1.25x to 1.68x, median 1.48x | 1.01x to 1.12x for the solvers; 1.18x to 1.70x, median 1.36x for the eigensolver alone, and that figure is the basis copy |
| Prediction that the ceiling collapses at n >= 1e5 | Structural, from per-call overhead | Refuted. The share is flat in n, because the cost scales with `n * ncv` |
