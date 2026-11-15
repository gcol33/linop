# Read the certificate off a result

Every verb in this package that computes something returns a certificate
with it, and each puts it where its own return value has room. A solve
carries it in an attribute, so the result stays the matrix a matrix
solve returns; the spectral verbs carry it in a field of their result;
\[verify()\] returns one directly. This reads it off any of them.

## Usage

``` r
certificate(x, ...)
```

## Arguments

- x:

  A result from \[solve()\], \[eigs()\], \[svds()\] or \[verify()\].

- ...:

  Unused.

## Value

An object of class \`linop_certificate\`.

## Details

A certificate belongs to a result and not to an operator. An operator on
its own has nothing to certify: the claims a certificate makes are about
a computation, and the one claim an operator can make alone is that it
satisfies the contract, which is \[verify()\] and returns a certificate
of its own.

A solve's attribute survives arithmetic on the solution and is dropped
by indexing, \`as.numeric()\`, \`c()\` and the reductions. So the cost
of \[solve()\] returning what a matrix solve returns is not that the
certificate is easily lost: it is that \`2 \* x\` still carries a
certificate describing \`x\`.

## Examples

``` r
A <- linop(diag(c(2, 3, 4)), properties = c(hermitian = TRUE,
                                            positive_definite = TRUE))
certificate(solve(A, c(2, 3, 4)))
#> solved by gmres (no capability is established, and gmres requires none)
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         computation      identity              1
#> residual                       pass         computation      identity              1
#> backward error                 pass         computation      identity              1
#> convergence                    pass         computation      identity              1
#> forward error                  not_checked  computation      identity              -
#> ------------------------------------------------------------------------------------------
#> overall                        qualified   no deterministic bound on: forward error
certificate(eigs(A, k = 1))
#> solved by primme (hermitian operator, PRIMME with implicit restarting (highest priority))
#> ------------------------------------------------------------------------------------------
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> arithmetic floor               pass         computation      identity              1
#> residual                       pass         computation      identity              1
#> orthogonality                  pass         computation      identity              1
#> backward error                 pass         computation      identity              1
#> target identity                not_checked  computation      identity              -
#> convergence                    pass         computation      identity              1
#> forward error                  pass         theorem          deterministic_bound   1
#> ------------------------------------------------------------------------------------------
#> overall                        qualified   no deterministic bound on: target identity
certificate(verify(A))
#> check                          status       source           guarantee             conf
#> ------------------------------------------------------------------------------------------
#> adjoint consistency            pass         computation      identity              1
#> view identities                pass         computation      identity              1
#> real transpose equals adjoint  pass         computation      identity              1
#> linearity                      pass         computation      identity              1
#> block consistency              pass         computation      identity              1
#> shape and dtype                pass         computation      identity              1
#> complex linearity              pass         computation      identity              1
#> declared capabilities          pass         probe            heuristic             -
#> purity                         pass         computation      identity              1
#> gemm agreement                 pass         computation      identity              1
#> materialisation agreement      pass         computation      identity              1
#> ------------------------------------------------------------------------------------------
#> overall                        pass   no deterministic bound on: declared capabilities
```
