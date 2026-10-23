# Estimated flops per column

An ordering heuristic only. S0.2 measured cases where a wider block is
\*slower\* per column for memory-bound operators, which a flop count
cannot express, so this must never be read as a predictor of wall time.

## Usage

``` r
linop_cost(A)
```

## Arguments

- A:

  A \`linop\`.

## Value

A numeric flop estimate.
