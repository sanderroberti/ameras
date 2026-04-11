# Parameter transformations

``` r
  library(ameras)
#> Loading required package: nimble
#> nimble version 1.4.2 is loaded.
#> For more information on NIMBLE and a User Manual,
#> please visit https://R-nimble.org.
#> 
#> Attaching package: 'nimble'
#> The following object is masked from 'package:stats':
#> 
#>     simulate
#> The following object is masked from 'package:base':
#> 
#>     declare
```

## Introduction

A transformation can be used to reparametrize parameters internally
(i.e., such that the likelihoods are evaluated at
`transform(parameters)`, where `parameters` are unconstrained), and
should be specified when fitting linear excess relative risk and
linear-exponential models to ensure nonnegative odds/risk/hazard. The
included function `transform1` applies an exponential transformation to
the desired parameters, see below. When supplying a function to
`transform`, this should be a function of the full parameter vector,
returning a full (transformed) parameter vector. In particular, the full
parameter vector contains parameters in the following order (see
[`?ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
the model definitions):
$\alpha_{0},{\mathbf{α}},\beta_{1},\beta_{2},{\mathbf{β}}_{m1},{\mathbf{β}}_{m2},\sigma$,
where $\mathbf{α}$, ${\mathbf{β}}_{m1}$ and ${\mathbf{β}}_{m2}$ can be
vectors, with lengths matching `X` and `M`, respectively. $\sigma$ is
only included for the linear model (Gaussian family), and no intercept
is included for the Cox and conditional logistic models. For the
multinomial model, the full parameter vector is the concatenation of
$Z - 1$ parameter vectors in the order as given above, where $Z$ is the
number of outcome categories. When no transformation is specified and
the linear ERR model is used, `transform1` is used for ERR parameters
$\beta_{1}$ and $\beta_{2}$ by default, with lower limits $- 1/max(D)$
for linear dose-response and
$\left( 0, - 1/max\left( D^{2} \right) \right)$ for linear-quadratic
dose-response, respectively (see below). For the linear-exponential
model, a lower limit of 0 is used for $\beta_{1}$, and no transformation
is used for $\beta_{2}$. If effect modifiers `M` are specified, no
transformation is used for those parameters. When negative RRs are
obtained during optimization, an error will be generated and a different
transformation or bounds should be used. All output is returned in the
original parametrization given in
[`?ameras`](https://ameras.sanderroberti.com/reference/ameras.md). The
Jacobian of the transformation (`transform.jacobian`) is required when
using a transformation with methods other than BMA. For `transform1`,
the Jacobian is given by `transform1.jacobian`.

## Exponential transformation using transform1

The included function `transform1` applies the exponential
transformation
$f\left( \theta_{i} \right) = \exp\left( \theta_{i} \right) + LB_{i}$ to
one or multiple components of parameter vector $\mathbf{θ}$, where
$LB_{i}$ are lower limits that can be different for each component. In
particular, a vector of indices of parameters to be transformed and a
vector of corresponding lower bounds LB can be supplied to arguments
`index.t` and `lowlimit`, respectively, resulting in transformed
parameters
$f\left( \theta_{i} \right) = \exp\left( \theta_{i} \right) + \text{LB}_{i}$.

In particular, `transform1` and `transform1.jacobian` are defined as
follows:

``` r
transform1 <- function(params, index.t=1:length(params), lowlimit=rep(0,length(index.t)), 
                       boundcheck=FALSE, boundtol=1e-3, ...){
  
  if(length(index.t)!=length(lowlimit)) 
    stop("Length mismatch between index.t and lowlimit")
  if(any(!(index.t %in% 1:length(params)))) 
    stop("Incorrect indices for transformation specified")
  params[index.t] <- exp(params[index.t]) + lowlimit
  if(boundcheck){
    if(any(params[index.t]-lowlimit < boundtol)) 
      warning(paste0("WARNING: one or multiple parameter estimates within ", boundtol, " of 
        lower bounds. Try different bounds or starting values."))
  }
  return(params)
}

transform1.jacobian <- function(params, index.t=1:length(params), ...){ 
  if(any(!(index.t %in% 1:length(params)))) 
    stop("Incorrect indices for transformation specified")
  grad <- rep(1, length(params))
  grad[index.t] <- exp(params[index.t])
  if(length(params)>1){
    return(diag(grad))
  } else{
    return(matrix(grad))
  }
}
```

## Defining a custom transformation

If you wish to supply your own transformation, it is helpful to start
from the definition of `transform1`. It is also important to keep the
following in mind:

1.  At a minimum, the function should take arguments `params` and `...`,
    where `params` is the full parameter vector. Parameters are in a
    specific order (see above), in case of doubt it is always possible
    to run `ameras` with the default settings to verify correct the
    order from the result.
2.  Extra arguments can be used and should be supplied when calling
    `ameras`
3.  Optional: if `boundcheck` is an argument to the function, it should
    be a logical. When transforming parameters after the optimization,
    the transformation is called with `boundcheck=TRUE`. If `boundcheck`
    is not an argument of the function, this is ignored.
4.  Also define the Jacobian of the transformation, keeping in mind
    point 1 above.

See the definition of `transform1` above for an example of how to apply
the transformation only to specific parameters, and how to use extra
arguments.

As an example, suppose instead of the exponential transformation from
`transform1`, for the parameters $\beta_{1}$ and $\beta_{2}$ we wish to
use the sigmoid transformation
$\left. f:{\mathbb{R}}\rightarrow(a,b) \right.$ given by
$$f_{a_{i},b_{i}}\left( \theta_{i} \right) = a_{i} + \left( b_{i} - a_{i} \right)\frac{1}{1 + \exp\left( - \theta_{i} \right)}.$$
Then, using `transform1` as a starting point, we can define the
transformation as follows (note that since $\mathbf{a}$ and $\mathbf{b}$
act as bounds, we use an updated bound check):

``` r
transform.sigmoid <- function(params, index.t=1:length(params), a=rep(0,length(index.t)), 
                              b=rep(1,length(index.t)), boundcheck=FALSE, boundtol=1e-3, ...){
  
  if(length(index.t)!=length(a) | length(index.t) != length(b)) 
    stop("Length mismatch between index.t, a, and b")
  if(any(!(index.t %in% 1:length(params)))) 
    stop("Incorrect indices for transformation specified")
  
  params[index.t] <- a + (b-a) * 1/(1+exp(-1*params[index.t]))
  if(boundcheck){
    if(any( (params[index.t]-a < boundtol) | (b-params[index.t] < boundtol))) 
      warning(paste0("WARNING: one or multiple parameter estimates within ", boundtol, 
      " of bounds. Try different bounds or starting values."))
  }
  return(params)
}
```

Next, noting that
$df_{a_{i},b_{i}}/d\theta_{i} = \left( b_{i} - a_{i} \right)\exp\left( - \theta_{i} \right)/\{ 1 + \exp\left( - \theta_{i} \right)\}^{2}$,
we can define the Jacobian as follows, using `transform1.jacobian` as a
starting point:

``` r
transform.sigmoid.jacobian <- function(params, index.t=1:length(params), 
                                       a=rep(0,length(index.t)), b=rep(1,length(index.t)), ...){ 
  if(length(index.t)!=length(a) | length(index.t) != length(b)) 
    stop("Length mismatch between index.t, a, and b")
  if(any(!(index.t %in% 1:length(params)))) 
    stop("Incorrect indices for transformation specified")
  grad <- rep(1, length(params))
  grad[index.t] <- (b-a)*exp(-1*params[index.t])/(1+exp(-1*params[index.t]))^2
  if(length(params)>1){
    return(diag(grad))
  } else{
    return(matrix(grad))
  }
}
```

Now let us try this transformation on the example data using regression
calibration for a linear-quadratic ERR model. Note that all parameters
are returned after transformation, and so there should be no difference
between using `transform1` and `transform.sigmoid`. First, we fit the
model using the sigmoid transformation:

``` r
data(data, package="ameras")
dosevars <- paste0("V", 1:10)
fit.ameras.sigmoid <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", deg=2, doseRRmod = "ERR", methods="RC",
                         transform=transform.sigmoid, transform.jacobian=transform.sigmoid.jacobian,
                         index.t=4:5)
#> Fitting RC
#> Obtaining profile likelihood CI for dose
#> Warning in ameras.rc(family = family, dosevars = dosevars, data = data, :
#> WARNING: Lower bound for dose is < 0 and may not exist if rescaling the
#> variable does not help
#> Obtaining profile likelihood CI for dose_squared
summary(fit.ameras.sigmoid)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     X = c("X1", "X2"), methods = "RC", deg = 2, doseRRmod = "ERR", 
#>     transform = transform.sigmoid, transform.jacobian = transform.sigmoid.jacobian, 
#>     index.t = 4:5)
#> 
#> Total run time: 5.2 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     5.2
#> 
#> Summary of coefficients by method:
#> 
#>  Method         Term Estimate      SE CI.lowerbound CI.upperbound
#>      RC  (Intercept) -0.87358 0.09758            NA            NA
#>      RC           X1  0.44587 0.07672            NA            NA
#>      RC           X2 -0.33552 0.09610            NA            NA
#>      RC         dose  0.04875 0.21281        0.0000        0.5130
#>      RC dose_squared  0.28764 0.08100        0.1328        0.4121
```

Next with default settings, using `transform1`:

``` r
fit.ameras.transform1 <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", deg=2, doseRRmod = "ERR", methods="RC")
#> Fitting RC
#> Obtaining profile likelihood CI for dose
#> Warning in ameras.rc(family = family, dosevars = dosevars, data = data, :
#> WARNING: Lower bound for dose is < 0 and may not exist if rescaling the
#> variable does not help
#> Obtaining profile likelihood CI for dose_squared
summary(fit.ameras.transform1)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     X = c("X1", "X2"), methods = "RC", deg = 2, doseRRmod = "ERR")
#> 
#> Total run time: 5.5 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     5.5
#> 
#> Summary of coefficients by method:
#> 
#>  Method         Term Estimate      SE CI.lowerbound CI.upperbound
#>      RC  (Intercept) -0.87359 0.09759            NA            NA
#>      RC           X1  0.44587 0.07672            NA            NA
#>      RC           X2 -0.33552 0.09610            NA            NA
#>      RC         dose  0.04878 0.21283        0.0000        0.5115
#>      RC dose_squared  0.28763 0.08100        0.1325        0.4108
```

As expected, there is no difference between the results.
