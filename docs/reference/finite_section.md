# Truncate an operator on a sequence space

`P_n A P_n` on the indices `-n` to `n`, as an operator of dimension
`2n + 1`. The result is an ordinary `linop`: it has a matrix,
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) runs on
it, and the arithmetic composes with anything else.

## Usage

``` r
finite_section(A, n)
```

## Arguments

- A:

  An operator on a sequence space, from
  [`linop_jacobi()`](https://gillescolling.com/linop/reference/linop_jacobi.md)
  or
  [`linop_space()`](https://gillescolling.com/linop/reference/linop_space.md).

- n:

  Half-width of the window. The result is `2n + 1` by `2n + 1`.

## Value

A `linop`.

## Details

What it returns holds the operator it truncates, so the printed tree
shows the infinite operator underneath the finite one. That is also
where the certificate comes from:
[`eigs()`](https://gillescolling.com/linop/reference/eigs.md) on a
finite section certifies a statement about `A`, not about the
truncation, because the truncation is not the thing anyone wanted to
know about.

`n` has to be large enough that the window sits strictly inside it.
Outside the window the coefficients are at their limits and the
eigenvalue equation is the free recurrence, and every closed form in the
certificate reads the tail that way; a truncation that cuts through the
perturbation has none of them.

## Examples

``` r
H <- linop_jacobi(diagonal = c(1.5, -0.5, 2))
F <- finite_section(H, n = 12)
dim(F)
#> [1] 25 25
print(F)
#> <linop> 25 x 25, double
#>   hermitian, symmetric, real
#> 
#> section |j| <= 12
#> `- jacobi on l^2(Z), sigma_ess [-2, 2]
```
