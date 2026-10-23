# Collapse cheap subexpressions into dense leaves

Cost optimisation as an inspectable transformation rather than hidden
behaviour, so a wrong cost model costs a suboptimal plan and never a
surprise allocation.

## Usage

``` r
collapse(A, expected_applies = 100, max_entries = 1e+06)
```

## Arguments

- A:

  A \`linop\`.

- expected_applies:

  How many times the operator will be applied.

- max_entries:

  Largest leaf this is allowed to build.

## Value

A \`linop\` carrying a \`report\` attribute.
