# Jacobian of the exponential parameter transformation

Computes the Jacobian matrix of
[transform1](https://ameras.sanderroberti.com/reference/transform1.md).
Note that lower limits do not need to be specified as the Jacobian is
independent of those

## Usage

``` r
transform1.jacobian(params, index.t=1:length(params), ... )
```

## Arguments

- params:

  input parameter vector (before transformation) to evaluate the
  Jacobian at

- index.t:

  indices of parameters to be transformed (default all)

- ...:

  not used

## Value

Jacobian matrix.

## Examples

``` r
params <- c(.1, .5, 1)
transform1.jacobian(params)
#>          [,1]     [,2]     [,3]
#> [1,] 1.105171 0.000000 0.000000
#> [2,] 0.000000 1.648721 0.000000
#> [3,] 0.000000 0.000000 2.718282
```
