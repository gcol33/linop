# Does this evidence satisfy a requirement?

Recurses into \`depends_on\`, so a composite built by an unconditional
rule from user-declared inputs fails any requirement the inputs would
have failed directly. Without the recursion, \`crossprod(A) +
crossprod(B)\` would report bare \`construction\` and launder the
declarations underneath it.

## Usage

``` r
evidence_satisfies(actual, required)
```

## Arguments

- actual:

  Evidence from \[evidence()\], or \`NULL\`.

- required:

  A \[requirement()\].

## Value

\`TRUE\` or \`FALSE\`.
