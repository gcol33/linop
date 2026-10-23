# How fast an eigenvector decays

Outside the window the eigenvalue equation is the free recurrence \`b
(u\[j-1\] + u\[j+1\]) + a u\[j\] = lambda u\[j\]\`, whose solutions are
\`z^j\` for the two roots of \`z^2 - mu z + 1\` with \`mu = (lambda -
a_inf) / b_inf\`. The roots are reciprocal, so for \`\|mu\| \> 2\`
exactly one lies inside the unit disc, and membership in \`l^2\` forces
it on each side. An eigenvector is therefore \*exactly\* geometric
outside the window, at the rate this returns, rather than asymptotically
so.

## Usage

``` r
decay_rate(A, lambda)
```

## Arguments

- A:

  A \`linop\` from \[linop_jacobi()\], or a \[finite_section()\] of one.

- lambda:

  One or more real values.

## Value

A numeric vector in \`(0, 1)\`, with \`NA\` where \`lambda\` lies in the
essential spectrum.

## Details

\`NA\` for a value inside the essential spectrum: there both roots lie
on the unit circle, nothing decays, and there is no eigenvector to
decay.

## Examples

``` r
H <- linop_jacobi(diagonal = 1)
## the golden-ratio conjugate, for the eigenvalue sqrt(5)
decay_rate(H, sqrt(5))
decay_rate(H, 0)
```
