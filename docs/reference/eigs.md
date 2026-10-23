# Eigenvalues and eigenvectors of a hermitian operator

Lanczos with full reorthogonalisation, run against the operator's
action. The result records what converged, what residual each pair
reached, and what argument establishes the bound on the eigenvalues.

## Usage

``` r
eigs(
  A,
  k,
  which = "largest",
  sigma = NULL,
  B = NULL,
  tol = 1e-10,
  maxit = NULL,
  ncv = NULL,
  v0 = NULL,
  seed = 1L,
  method = "auto",
  inner = list(),
  section = list(),
  floor_const = SOLVE_FLOOR_CONST,
  norm_control = list()
)
```

## Arguments

- A:

  A square hermitian \`linop\`. No adjoint is required.

- k:

  How many eigenpairs.

- which:

  Which end of the spectrum: \`"largest"\` or \`"smallest"\` by
  magnitude, \`"largest_algebraic"\` or \`"smallest_algebraic"\` by
  value. RSpectra's \`"LM"\`, \`"SM"\`, \`"LA"\` and \`"SA"\` are
  accepted for the same four.

- sigma:

  A real shift. Given one, the eigenvalues nearest \`sigma\` are sought,
  by running the same recurrence on \`(A - sigma I)^-1\` through the
  package's own solvers. \`which\` does not apply and passing both is an
  error.

- B:

  Reserved for the generalized problem \`A x = lambda B x\`, which is
  not in this version.

- tol:

  Convergence tolerance. Every method holds each pair, measured on
  \`A\`, to the backward error \`\|\|A x - mu x\|\| / (\|\|A\|\|
  \|\|x\|\|)\` and to the relative residual \`\|\|A x - mu x\|\| /
  (\|mu\| \|\|x\|\|)\` wherever \`\|mu\|\` is above the arithmetic
  floor, each plus that floor: the certificate's \`backward error\` and
  \`residual\` rows.

- maxit:

  Budget of applies across all restarts, of the operator the iteration
  runs on: \`A\`, or under \`sigma\` the inner solve.

- ncv:

  Subspace size per restart. Storage is \`ncv\` vectors. \`NULL\` takes
  the engine's default: \`max(2k + 1, k + 20)\` for \`"lanczos"\`, and
  PRIMME's own basis size for \`"primme"\`.

- v0:

  Starting vector, or \`NULL\` for a seeded random one.

- seed:

  Seed for the starting vector. The caller's random stream is restored
  afterwards.

- method:

  \`"auto"\`, \`"lanczos"\` or \`"primme"\`. \`"lanczos"\` is the
  reference implementation in this package. \`"primme"\` is the PRIMME
  library, available when the \`linop.primme\` package is installed,
  which restarts implicitly and so keeps a basis of fixed small size.
  \`"auto"\` takes PRIMME when it is available and Lanczos otherwise;
  the certificate's \`dispatch\` records every engine considered and
  why. Whichever runs, the reported values are Rayleigh quotients
  measured on \`A\` and the certificate is the same.

- inner:

  Arguments for the inner solve when \`sigma\` is given, such as
  \`method\` or \`tol\`.

- section:

  Arguments for the truncation, when \`A\` is an operator on a sequence
  space and has no matrix: \`n_start\` for how much free tail the first
  section leaves beyond the window, and \`n_max\` for the widest section
  the search may reach. The width itself is not among them, because
  \`eigs(finite_section(A, n))\` is already how a caller fixes it.

- floor_const:

  \`c\` in the arithmetic floor of the certificate.

- norm_control:

  Arguments for the \`\|\|A\|\|\` estimate.

## Value

An object of class \`linop_eigen\`, with \`values\`, \`vectors\` and the
\`certificate\`. A run that exhausts its budget returns the best pairs
it reached and a \`fail\` certificate rather than an error, and it may
hold fewer than \`k\` of them: nothing is padded. On a sequence space it
also carries the width \`n\` that was chosen and the \`widths\` that
were tried.

## Details

The eigenvalue reported for each pair is the Rayleigh quotient \`x^H A x
/ (x^H x)\`, measured on \`A\` itself. It is the value that minimises
\`\|\|A x - mu x\|\|\`, so no other choice of \`mu\` certifies better.

\`target identity\` is \`not_checked\`, always. A small residual proves
the pair is an approximate eigenpair; that it is the largest, smallest
or nearest one requires an inertia count, an enclosure or a separation
bound, and matrix-free there is generally none. See
\`vignette("certificates")\`.

## Examples

``` r
n <- 40
L <- linop(function(X) rbind(X[-1, , drop = FALSE], 0) +
                       rbind(0, X[-n, , drop = FALSE]) - 2 * X,
           dim = c(n, n), properties = c(hermitian = TRUE))
fit <- eigs(L, k = 3, which = "largest_algebraic")
fit$values

## on an operator with no matrix, where the width is chosen rather than given
H <- linop_jacobi(diagonal = 1)
fit <- eigs(H, k = 1, which = "largest_algebraic")
fit$n
fit$values - sqrt(5)
```
