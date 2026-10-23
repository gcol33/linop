# State what evidence a solver will accept

State what evidence a solver will accept

## Usage

``` r
requirement(
  sources = EVIDENCE_SOURCES,
  guarantees = EVIDENCE_GUARANTEES,
  min_confidence = NA_real_
)
```

## Arguments

- sources:

  Acceptable evidence sources.

- guarantees:

  Acceptable guarantee kinds.

- min_confidence:

  Minimum confidence, or \`NA\` to not require one.

## Value

An object of class \`linop_requirement\`.
