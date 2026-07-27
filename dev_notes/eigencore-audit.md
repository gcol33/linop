# S0.4 -- eigencore source audit

Audited 2026-07-27. **eigencore 1.0.2**, published to CRAN 2026-07-25, Bradley Buchsbaum.
Source: `eigencore_1.0.2.tar.gz` from cloud.r-project.org.
Size: 18,358 lines of R across 30 files, plus ~700 KB of C++ in `src/` (21 files).
`Depends: R (>= 4.1)`, `Imports: Matrix, methods`, `NeedsCompilation: yes`,
`SystemRequirements: C++17`.

**Read in full:** `R/operator.R` (453), `R/classes.R` (316). **Read targeted:**
`R/operator_algebra.R` propagation sites. **Not read:** the numerical kernels
(`reference_*.R`, `src/*.cpp`) beyond their interfaces. Every claim below carries a line
citation into a file that was read. The audit's question is whether eigencore already
provides linop's Layer 1 contract; that is answered by the operator model and the
capability/certificate design, not by Lanczos internals.

Plan section 2 said every claim about eigencore's internals was second-hand and that this
spike verifies or refutes them. Several were wrong, and in eigencore's favour.

---

## 1. eigencore already has more of Layer 1 than the plan assumes

The plan (section 2, "Where the gaps are") says what is missing is "a composable, typed,
capability-carrying operator shared across solver families", plus "operator-native Krylov"
and "block matvecs in any matrix-free interface". eigencore has all three.

`linear_operator()` (operator.R:29-49) constructs an S3 object of class
`eigencore_operator` (operator.R:47) with slots `dim`, `apply`, `apply_adjoint`, `dtype`,
`structure`, `name`, `metadata`.

### The apply signature is already linop's tier-2 GEMM ABI

```r
apply = function(X, alpha = 1, beta = 0, Y = NULL)          # operator.R:20, 78, 315
```

computing `alpha * A %*% X + beta * Y`. Plan section 5.2 proposes exactly this as the
optional tier-2 `apply_gemm(X, Y, alpha, beta)`. eigencore requires it as the *only* tier.
The docstring says "Create a block-native linear operator" (operator.R:1) and `X` is a
block throughout.

So the plan's "block matvecs in any matrix-free interface" gap is closed by eigencore,
and its two-tier design is a usability improvement over eigencore's one-tier, not a new
capability.

### It has an adapter generic for third parties

`as_operator()` is an S3 generic (operator.R:60-62) with methods for `matrix`
(operator.R:70), `Matrix` classes dispatched by `inherits` in `.default`
(operator.R:92-107), and identity on itself (operator.R:65). A third party writes
`as_operator.myclass()`. This is structurally the same adoption mechanism as linop's
`linop.myclass()` (plan section 1.1).

`adjoint()` is likewise an S3 generic (operator.R:136-138) with the operator method at
operator.R:141.

### Complex support is real but partial, same wall linop hit

Complex dense goes through a dedicated `zgemm` path
(`complex_dense_matrix_as_operator`, operator.R:176-197). Complex *sparse* is explicitly
refused with a "future scope" error (operator.R:110-129). This is the same Matrix 1.7.5
limitation recorded in S0.1 section 7, reached independently.

---

## 2. What linop genuinely adds, verified by full read

### 2a. Capabilities carry no evidence, and there are only two of them

```r
general   <- function() structure(list(kind = "general"),   class = "eigencore_structure")
hermitian <- function() structure(list(kind = "hermitian"), class = "eigencore_structure")
```
(classes.R:253-255, 262-264)

That is the entire capability system: **one field, two values.** There is no
`positive_definite`, no `orthogonal`, no `triangular`, no `full_rank`, no three-valued
`NA`, and above all no record of *why* a flag holds.

The consequence is visible at the point where a user asserts symmetry:

```r
if (isTRUE(validate)) check_adjoint(A, trials = 5L, tol = tol)
A$structure <- hermitian()
```
(operator_algebra.R:408-411)

After this line, an operator whose Hermitian flag came from a bare user assertion is
indistinguishable from one that came from five random probes, which is in turn
indistinguishable from `crossprod_operator()`'s unconditional construction
(operator_algebra.R:441). Plan section 5.3 exists precisely to keep those three apart, and
its argument that CG must not accept "ten positive Rayleigh quotients" as evidence of
positive definiteness has a concrete target here.

**This is linop's clearest differentiator and it survives contact with the source.**

### 2b. Propagation is right where it overlaps, and thinner elsewhere

eigencore's composition rules, and how they compare to plan section 5.6:

| eigencore | line | rule | linop 5.6 |
|---|---|---|---|
| `compose(A, B)` | operator_algebra.R:54 | always `general()` | agrees (`NA` unless provably adjoint) |
| `operator_sum` | operator_algebra.R:101 | hermitian iff all terms hermitian | agrees |
| scalar `scale` | operator_algebra.R:145 | **preserves `A$structure` unconditionally** | linop: hermitian iff `a` is *real* |
| `scale_rows`, `scale_cols` | operator_algebra.R:188, 232 | `general()` | n/a |
| `center` | operator_algebra.R:351 | `general()` | n/a |
| `crossprod_operator` | operator_algebra.R:441 | always `hermitian()` | agrees |

The comment at operator_algebra.R:48-53 is worth quoting, because it shows the same
reasoning the plan reaches independently:

> "AB is Hermitian only if A and B are Hermitian AND commute (e.g. when B == A and A is
> Hermitian). The previous heuristic compared A$name to B$name, but names are
> user-controlled metadata and must never drive certificate-relevant flags."

The scalar row is the one divergence. Preserving `hermitian` through multiplication by a
scalar is correct for real scalars and wrong for complex ones, since `(cA)^H = c-bar A^H`.
Whether eigencore restricts scalars to real was **not verified** -- I did not read the
argument validation in that function. Do not repeat this as a defect without reading it;
it is flagged here only as a place where linop's rule is stated more carefully.

### 2c. No linear solvers at all

Searched `R/` and `src/` for `minres|lsqr|lsmr|bicgstab|gmres`: **zero hits.** eigencore is
partial eigenvalue and singular value computation only; the exported verbs are
`eig_partial`, `eig_full`, `svd_partial`, `eigs`, `eigs_sym`, `generalized_svd`,
`generalized_schur`.

The plan's claim that MINRES, LSQR and LSMR are missing from CRAN is therefore not
threatened by eigencore, and the whole `solve()` half of linop is uncontested here.

### 2d. Preconditioners are bare functions

LOBPCG takes `preconditioner`, documented as "Optional function taking a residual block
and returning a preconditioned block with the same dimensions" (classes.R:203-205), with
validation only that it is a function (classes.R:219-221). There is no preconditioner
*object*, no `fixed` flag, no declared `hermitian`/`positive_definite`, and therefore no
way for a solver to refuse an inadmissible one. Plan section 4.3's enforcement table has
no counterpart.

### 2e. The API is the constructor-object DSL the plan rejects

eigencore's user-facing vocabulary is built from descriptor constructors: `largest()`,
`smallest()`, `largest_magnitude()`, `smallest_magnitude()`, `nearest(sigma)`,
`largest_real()`, `both_ends(k_low, k_high)` (classes.R:11-107); `auto()`, `lanczos()`,
`golub_kahan()`, `randomized()`, `lobpcg()`, `shift_invert(sigma)` (classes.R:119-247);
`general()`, `hermitian()`, `euclidean(dim)` (classes.R:253-275).

Plan section 1.1 says:

> "Not `target = largest()` or `method = shift_invert(...)`. Constructor objects are
> architecturally tidy and impose a miniature language on the user."

That sentence is now a direct, and evidently unwitting, contrast with the incumbent. It is
a legitimate and defensible difference -- linop takes strings and R's own generics -- but
the plan should state it as a deliberate choice against a named alternative rather than as
a self-evident principle. **And it must be stated collegially**: eigencore's descriptors
are typed, discoverable and self-documenting, which are real advantages. The linop tradeoff
is fewer names to learn against less structure to introspect.

### 2f. eigencore does not implement R's matrix generics

Nothing in eigencore makes an operator behave like a matrix. There is no `%*%`, `t`,
`crossprod`, `solve` or `dim` method on `eigencore_operator`; the class is a plain list
with a `print` method (operator.R:200-207). Users call `as_operator()` then `eig_partial()`.

linop's "a `linop` behaves like a matrix under R's normal arithmetic" is therefore a
genuine difference in kind, not degree, and after S0.1 it is known to be achievable with
one `matrixOps` group method.

---

## 3. Certificates: overlapping claim, needs a separate check

eigencore's DESCRIPTION claims:

> "Every result carries a numerical certificate with residuals, a backward-error bound,
> orthogonality loss, and a pass/fail flag, and bounds that can only be estimated are
> reported as such rather than passed."

That last clause is close to plan section 6's central discipline. `certification.R` (1,045
lines) and `validation.R` (1,087 lines) were **not read**. Exports include `certificate`,
`backward_error`, `diagnostics`, `check_adjoint`.

**Open question for Phase 2, not answerable from this audit:** whether eigencore's
certificate separates *status* from *evidence class* the way plan section 6 requires, and
whether it has an equivalent of `target_identity` defaulting to `not_checked` (section
6.1). If it does, linop's certificate design is a refinement rather than a new idea, and
the plan's framing must change accordingly. Read `certification.R` in full before writing
any certificate copy for the README or a paper.

---

## 4. Consequences for the plan

| Plan location | Says | Evidence |
|---|---|---|
| 2, "Where the gaps are" | missing: composable typed operator, operator-native Krylov, block matvecs | all three present in eigencore; operator.R:29-49, block-native by design |
| 2, "Not verified" | eigencore claims are second-hand | now verified; several were understatements |
| 5.2 | two-tier apply is the new idea | eigencore already requires the GEMM tier; linop's contribution is making tier 1 the *easy* path |
| 5.3 | evidence-bearing capabilities | **confirmed differentiator**; eigencore has one field, two values, no provenance (classes.R:253-264, operator_algebra.R:408-411) |
| 4.3 | preconditioner contract enforcement | **confirmed differentiator**; eigencore takes a bare function (classes.R:203-221) |
| 7.1 | MINRES/LSQR/LSMR absent | **confirmed**; zero hits in eigencore |
| 1.1 | rejects constructor-object DSL | now a contrast with a named incumbent; restate as a deliberate tradeoff, collegially |
| 14 | "real win condition: eigencore implementing the linop contract" | more plausible than the plan realised: `as_operator`/`adjoint` are already S3 generics with the same shape |
| 15 | reach out to Buchsbaum after v0.1 | the interface overlap is large enough that the conversation is more useful *before* the algebra freezes |

## 5. Presentation constraint, still binding

Sections 1 and 2 above are quotable because they cite lines from files read in full.
Section 3 is **not** quotable anywhere until `certification.R` and `validation.R` are read.
Nothing in this audit may be phrased as a deficiency of eigencore in a public artefact; the
plan's section 15 rule holds, and the honest framing throughout is that eigencore solves an
overlapping problem with a different set of tradeoffs.
