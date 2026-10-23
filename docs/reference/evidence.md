# Describe the argument behind a capability value

Describe the argument behind a capability value

## Usage

``` r
evidence(source, guarantee, confidence = NA_real_, depends_on = list())
```

## Arguments

- source:

  One of \`construction\`, \`adapter_contract\`, \`user_declaration\`,
  \`probe\`, \`computation\`, \`theorem\`.

- guarantee:

  One of \`identity\`, \`deterministic_bound\`, \`probabilistic_bound\`,
  \`estimate\`, \`heuristic\`.

- confidence:

  \`1\`, a stated probability, or \`NA\`.

- depends_on:

  List of evidence objects this one rests on. Empty when the rule is
  unconditional.

## Value

An object of class \`linop_evidence\`.
