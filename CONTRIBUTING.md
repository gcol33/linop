# Contributing to linop

Thanks for taking the time. This document covers how to get the package running
locally, what the tests expect, and the two contributions that are most useful
right now.

## Code of Conduct

Participation is governed by the [Code of Conduct](CODE_OF_CONDUCT.md).

## Getting set up

```r
# clone, then from the package root
devtools::load_all(".")
devtools::test(".")
```

The package has no `Imports` and no compiled code, so there is nothing to build.
`Matrix` is suggested and is needed for the sparse tests. R 4.4.0 or newer is
required, because `crossprod()` and `tcrossprod()` became S3 generic in 4.4.0.

After changing anything in `R/`:

```r
devtools::load_all(".")
roxygen2::roxygenise(".", clean = TRUE)   # regenerates NAMESPACE and man/
devtools::test(".")
```

Before opening a pull request:

```
R CMD build .
R CMD check --as-cran --no-manual linop_0.0.0.9000.tar.gz
```

`R CMD check` should end at `Status: 1 NOTE`, which is the new-submission and
development-version note. Anything beyond that is a regression.

## What the tests expect

Tests are recovery and contract tests. A test that only asserts a shape, a class
or a dimension does not add coverage here, because the package's claims are about
numerical agreement and about soundness.

- A new **node type** needs the conformance suite run against it in both dtypes,
  and its four apply modes checked against a materialised reference.
- A new **propagation rule** needs a brute-force soundness check: build random
  operators, materialise them, and confirm the claimed capability holds of the
  materialised matrix over many draws. The existing propagation suite is 9,892
  such checks and is the pattern to follow.
- A new **solver** needs parameter recovery against closed-form truth, run to
  convergence rather than to a small iteration cap, plus certificate coverage
  over at least 20 seeds.

Two suites will fail deliberately if you are not expecting them:

- `test-api-budget.R` asserts the exported set exactly. Adding an export fails it
  until the budget is edited, which is the point. Never edit the budget as a side
  effect of adding a function.
- The same file asserts that the deferred node types are absent. `lowrank`,
  `kron`, stacks, `blockdiag`, `perm`, `power` and `inverse` are held back until
  an external adapter has exercised the registry.

The expression fuzzer defaults to 400 trees so a routine check stays inside
CRAN's time budget. To run it at full strength:

```r
Sys.setenv(LINOP_FUZZ_N = "10000")
testthat::test_file("tests/testthat/test-fuzzer.R")
```

## Two contributions that would help most

**An adapter for a storage format you already have.** Write
`linop.myclass()`, run `verify()`, and open an issue with the certificate. The
registry and the conformance suite exist to be exercised from outside the
package, and the deferred node types are waiting on exactly that. See
`vignette("adapters")`.

**A conformance check the suite is missing.** `verify()` tests the contract every
solver relies on. If you can write an operator that is wrong in a way the eleven
current checks do not notice, that is a bug report worth more than a fix. The
`complex linearity` check came from exactly that: a real operator that discarded
the imaginary part of its input passed the other ten.

## Documentation

Vignettes are executed at build time, so every printed result in them is real
output. Keep it that way: if a chunk cannot run in a few seconds, shrink the
example rather than setting `eval = FALSE`.

Roxygen and `.Rd` files are ASCII only. Superscripts, mathematical operators and
arrows outside Latin-1 break the LaTeX pass that CRAN and win-builder run, and
the failure shows up only there. Write `^2`, `<=`, `->` instead.

Comments in code describe what the code does and why the domain requires it.
Notes about history, status or process belong in the commit message or the pull
request.

## Reporting a bug

An operator that behaves incorrectly is best reported as a minimal `linop()` call
plus the `verify()` certificate it produces. For a numerical disagreement,
include the gap and the reference you compared against.
