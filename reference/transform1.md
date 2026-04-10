# Exponential parameter transformation

Applies exponential transformation \\f(\theta_i)=\exp(\theta_i)+L_i\\ to
one or multiple components of parameter vector \\\bm \theta\\, where
\\L_i\\ are lower limits that can be different for each component

## Usage

``` r
transform1(params, index.t=1:length(params), lowlimit=rep(0,length(index.t)), 
  boundcheck=FALSE, boundtol=1e-3, ... )
```

## Arguments

- params:

  full input parameter vector

- index.t:

  indices of parameters to be transformed (default all)

- lowlimit:

  lower limits to be applied (default zero), where the k-th component of
  `lowlimit` is applied to the k-th index in `index.t`

- boundcheck:

  whether to produce a warning when any of the transformed parameters
  are within `boundtol` of `lowlimit`

- boundtol:

  tolerance for producing a warning for reaching the boundary

- ...:

  not used

## Value

Transformed parameter vector.

## Examples

``` r
params <- c(.1, .5, 1)
transform1(params, lowlimit=c(0, -1, 1))
#> [1] 1.1051709 0.6487213 3.7182818
```
