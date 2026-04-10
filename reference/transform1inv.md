# Inverse of exponential parameter transformation

Inverse of `transform1` for the purpose of deriving initial values.

## Usage

``` r
transform1.inv(params, index.t=1:length(params), lowlimit=rep(0,length(index.t)), ... )
```

## Arguments

- params:

  full input parameter vector

- index.t:

  indices of parameters to be transformed (default all)

- lowlimit:

  lower limits to be applied (default zero), where the k-th component of
  `lowlimit` is applied to the k-th index in `index.t`

- ...:

  not used

## Value

Transformed parameter vector.

## Examples

``` r
params <- c(.1, .5, 1) # Desired initial values on original scale
transform1.inv(params, lowlimit=c(0, -1, 1)) # Initial values to use on transformed scale
#> [1] -2.3025851  0.4054651       -Inf
```
