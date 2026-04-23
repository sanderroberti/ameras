# Analyze multiple exposure realizations

Fit regression models accounting for exposure uncertainty using multiple
Monte Carlo exposure realizations. Six outcome model families are
supported. The first is the Gaussian family for continuous outcomes,
\$\$Y_i \sim N(\mu_i, \sigma^2),\$\$ with \\\mu_i = \alpha_0 + \bm X_i^T
\bm \alpha +\beta_1 D_i+\beta_2 D_i^2 + \bm M_i^T \bm \beta\_{m1}D_i +
\bm M_i^T \bm \beta\_{m2}D_i^2\\. Here \\\bm X_i\\ are covariates,
\\D_i\\ is the exposure with measurement error, and \\\bm M_i\\ are
binary effect modifiers. The quadratic exposure terms and effect
modification are optional.

For non-Gaussian families, three relative risk models for the main
exposure are supported, the usual exponential \\RR_i=\exp(\beta_1
D_i+\beta_2 D_i^2+ \bm M_i^T \bm \beta\_{m1}D_i + \bm M_i^T \bm
\beta\_{m2} D_i^2)\\ and the linear excess relative risk (ERR) model
\\RR_i= 1+\beta_1 D_i+\beta_2 D_i^2 + \bm M_i^T \bm \beta\_{m1}D_i + \bm
M_i^T \bm \beta\_{m2}D_i^2\\, where the quadratic and effect
modification terms are optional. Finally, the linear-exponential
relative risk model \\RR_i= 1+(\beta_1 + \bm{M}\_i^T \bm{\beta}\_{m1})
D_i \exp\\(\beta_2+ \bm{M}\_i^T \bm{\beta}\_{m2})D_i\\\\ is supported.

The second supported family is logistic regression for binary outcomes,
with probabilities \$\$p_i/(1-p_i)=RR_i\exp(\alpha_0+\bm X_i^T \bm
\alpha).\$\$

Third is Poisson regression for counts, \$\$Y_i \sim
\text{Poisson}(\mu_i),\$\$ where \\\mu_i=RR_i \exp(\alpha_0 +\bm X_i^T
\bm \alpha)\times \text{offset}\_i\\ with optional offset.

Fourth is proportional hazards regression for time-to-event data, with
hazard function \$\$h(t) = h_0(t)RR_i\exp(\bm X_i^T \bm \alpha),\$\$
with \\h_0\\ the baseline hazard.

Fifth is multinomial logistic regression for a categorical outcome with
\\Z\>2\\ outcome categories, with the last category as the referent
category (i.e., \\\alpha\_{0,Z}=\bm
\alpha\_{Z}=\beta\_{1,Z}=\beta\_{2,Z}=\bm \beta\_{m1,Z} = \bm
\beta\_{m2,Z}=0\\): \$\$P(Y_i=z)=RR_i\exp(\alpha\_{0,z}+\bm X_i^T \bm
\alpha\_{z})/\\1+\sum\_{s=1}^{Z-1} RR_i\exp(\alpha\_{0,s}+\bm X_i^T \bm
\alpha\_{s})\\\$\$

Sixth is conditional logistic regression for matched case control data,
for which \$\$P\left(Y_i = 1, Y_k = 0 \forall k \neq i \bigg\| \sum\_{i
\in \mathcal{R}} Y_i = 1\right) = RR_i\exp(\bm X_i^T \bm
\alpha)/\\\sum\_{k \in \mathcal{R}}RR_k\exp(\bm X_k^T \bm \alpha)\\,\$\$
where \\\mathcal{R}\\ is the matched set corresponding to individual
\\i\\.

Methods include regression calibration (Carroll et al. 2006
[doi:10.1201/9781420010138](https://doi.org/10.1201/9781420010138) ),
extended regression calibration (Little et al. 2023
[doi:10.1038/s41598-023-42283-y](https://doi.org/10.1038/s41598-023-42283-y)
), Monte Carlo maximum likelihood (Stayner et al. 2007
[doi:10.1667/RR0677.1](https://doi.org/10.1667/RR0677.1) ), frequentist
model averaging (Kwon et al. 2023
[doi:10.1371/journal.pone.0290498](https://doi.org/10.1371/journal.pone.0290498)
), and Bayesian model averaging (Kwon et al. 2016
[doi:10.1002/sim.6635](https://doi.org/10.1002/sim.6635) ).

## Usage

``` r
ameras(formula=NULL, data, family="gaussian", methods="RC", 
  Y=NULL, dosevars=NULL, doseRRmod=NULL, deg=NULL,
  M=NULL, X=NULL, offset=NULL, entry=NULL, exit=NULL,
  setnr=NULL,
  CI=NULL, params.profCI=NULL,
  maxit.profCI=NULL, tol.profCI=NULL,
  transform=NULL,
  transform.jacobian=NULL, inpar=NULL, loglim=1e-30, MFMA=100000, 
  prophaz.numints.BMA=10, ERRprior.BMA="doubleexponential", nburnin.BMA=5000, 
  niter.BMA=20000, nchains.BMA=2, thin.BMA=10, included.replicates.BMA=NULL, 
  optim.method="Nelder-Mead", control=NULL, keep.data=TRUE, ... )
```

## Arguments

- formula:

  an object of class `"formula"` containing the model specification. See
  Details.

- data:

  input data frame.

- family:

  outcome model family: `"gaussian"`, `"binomial"`, `"poisson"`,
  `"prophaz"`, `"multinomial"` or `"clogit"` (default `"gaussian"`).

- methods:

  character vector of one or multiple methods to apply. Options: `"RC"`,
  `"ERC"`, `"MCML"`, `"FMA"`, `"BMA"` (default `"RC"`).

- Y:

  **Deprecated**. Use the formula interface instead. Name or column
  index of the outcome variable for linear, binomial, Poisson,
  multinomial and conditional logistic models, or event indicator
  variable for the proportional hazards model.

- dosevars:

  **Deprecated**. Use the formula interface instead. Names or column
  indices of exposure replicate vectors.

- doseRRmod:

  **Deprecated**. Use the formula interface instead. The functional form
  of the dose-response relationship; options are exponential RR
  (`"EXP"`), linear ERR (`"ERR"`), or linear-exponential RR (`"LINEXP"`)
  (default `"ERR"`).

- deg:

  **Deprecated**. Use the formula interface instead. For
  `doseRRmod="ERR"` and `doseRRmod="EXP"`, whether to fit a linear
  (`deg=1`) or linear-quadratic (`deg=2`) dose-response model (default
  linear).

- M:

  **Deprecated**. Use the formula interface instead. Names or column
  indices of binary effect modifying variables (optional).

- X:

  **Deprecated**. Use the formula interface instead. Names or column
  indices of other covariates (optional).

- offset:

  **Deprecated**. Use the formula interface instead. Name or column
  index of offset variable for Poisson regression (optional).

- entry:

  **Deprecated**. Use the formula interface instead. Name or column
  index of left truncation time variable for proportional hazards
  regression (optional).

- exit:

  **Deprecated**. Use the formula interface instead. Name or column
  index of exit time variable, required when `family=prophaz`.

- setnr:

  **Deprecated**. Use the formula interface instead. Name or column
  index of integer-valued matched set variable, required when
  `family="clogit"`.

- CI:

  **Deprecated**. Use confint() to compute confidence intervals instead.
  Method for calculation of 95% confidence or credible intervals (see
  Details). For RC, ERC, and MCML, options are `"wald.orig"`,
  `"wald.transformed"`, `"proflik"` (default `"proflik"`). For FMA and
  BMA, options are `"percentile"` and `"hpd"` (default `"percentile"`).
  If `methods` contains at least one of RC, ERC, and MCML and at least
  one of FMA and BMA, `CI` must be length 2 and specify one method for
  RC, ERC, and MCML, and one for FMA and BMA (see Details).

- params.profCI:

  **Deprecated**. Use confint() to compute confidence intervals instead.
  When `CI="proflik"`, whether to obtain profile-likelihood CIs for all
  parameters (`"all"`) or only dose-related parameters (`"dose"`,
  default).

- maxit.profCI:

  **Deprecated**. Use confint() to compute confidence intervals instead.
  Maximum iterations for determining profile-likelihood CIs; passed to
  `uniroot` (default 20).

- tol.profCI:

  **Deprecated**. Use confint() to compute confidence intervals instead.
  Tolerance for determining profile-likelihood CIs; passed to `uniroot`
  (default `1e-2`).

- transform:

  function for internal parameter transformation (see Details).

- transform.jacobian:

  Jacobian of the transformation function (see Details).

- inpar:

  vector of initial values for log-likelihood optimization (optional).

- loglim:

  parameter used in likelihood computations to avoid taking the log of
  very small or negative numbers via `log(max(x, loglim))` (default
  `1e-30`).

- MFMA:

  number of samples for `"FMA"` to compute estimates and CIs (default
  100,000).

- prophaz.numints.BMA:

  for `methods="BMA"` with `family="prophaz"`, the number of
  subintervals with constant baseline hazard (default 10). Cut points
  are determined based on quantiles of the event time distribution among
  cases.

- ERRprior.BMA:

  prior for dose-related parameters when `doseRRmod="ERR"` or `"LINEXP"`
  and `methods="BMA"`. Options: `"truncated_normal"`,
  `"truncated_horseshoe"`, `"truncated_doubleexponential"`, `"normal"`,
  `"horseshoe"`, `"doubleexponential"`, see Details (default
  `"doubleexponential"`).

- nburnin.BMA:

  number of MCMC burn-in iterations for BMA (default 1,000).

- niter.BMA:

  number of MCMC iterations per chain for BMA (default 5,000).

- nchains.BMA:

  number of MCMC chains for BMA (default 2).

- thin.BMA:

  thinning rate for BMA (default 10).

- included.replicates.BMA:

  indices of exposure replicates used in BMA (defaults to all
  replicates).

- optim.method:

  method used for optimization by `optim`. Options are `"Nelder-Mead"`
  and `"BFGS"`. When using Nelder-Mead, a second optimization with BFGS
  is run to ensure an optimal fit.

- control:

  control list passed to `optim` (default `list(reltol=1e-10)`).

- keep.data:

  whether to attach data to the output object (default `TRUE`). When the
  data object is large, `keep.data` can be set to `FALSE` to preserve
  memory. The attached data is used to compute profile likelihood
  confidence intervals, but can also be supplied seperately when
  `keep.data=FALSE`. See
  [`confint`](https://ameras.sanderroberti.com/reference/confint.md).

- ...:

  other arguments, passed to functions such as `transform`.

## Value

The output is an object of class `amerasfit`. General components are
`call` (the function call to `ameras`), `formula` (the formula object
specifying the model), `num.rows` (the number of rows in `data`),
`num.replicates` (the number of dose replicates provided), `transform`
(the used transformation function, if applicable), `transform.jacobian`
(the used Jacobian function for the transformation, if applicable),
`other.args` (any other arguments passed to ...), `model` (a list
containing the specified model components parsed from the formula),
`CI.computed` (logical, whether confidence intervals have been attached
by [`confint`](https://ameras.sanderroberti.com/reference/confint.md)),
and `data` (either the data frame used for model fitting when
`keep.data=TRUE` or `NULL` otherwise).

For each method supplied to `methods`, the output contains a list with
components:

- coefficients:

  named vector of model coefficients.

- sd:

  named vector of standard deviations.

- runtime:

  string with the runtime in seconds.

For RC, ERC, and MCML the following additional output is included:

- vcov:

  covariance matrix for the full parameter vector.

- optim:

  a list object with results returned by optim. Components are `par`
  (raw parameters before applying a transformation if applicable),
  `hessian` (Hessian matrix for `par`), `convergence` (convergence code
  with 0 indicating convergence and 1 indicating that the maximal number
  of iterations was reached), and `counts` (the number of likelihood
  function evaluations used during optimization).

- loglik:

  log-likelihood value at the optimum.

For RC and ERC, the output additionally contains:

- ERC:

  logical, whether the output is for ERC (`ERC=TRUE`) or RC
  (`ERC=FALSE`).

For BMA the output additionally contains:

- samples:

  MCMC posterior samples, as obtained from `nimble`. This is a list
  object with `nchains.BMA` components, each a named matrix with the
  samples from one chain in its rows, with columns corresponding to
  model parameters.

- Rhat:

  data frame with two columns, `Rhat` and `n.eff`. The first column
  contains the Gelman-Rubin statistics \\\hat R \geq 1\\ that can be
  used to assess convergence of MCMC chains. A value of 1 indicates good
  convergence and values \\\>1.05\\ indicate poor convergence. The
  effective sample size `n.eff` is a measure of how many independent
  samples the auto-correlated MCMC samples correspond to. A low
  effective sample size indicates high correlations and/or poor mixing.

- included.replicates:

  indices of replicate exposures that were included to obtain the
  results.

- prophaz.timepoints:

  for `family="prophaz"`, time points defining the intervals on which
  the estimated baseline hazards is constant; these are
  `prophaz.numints.BMA + 1` time points covering the interval
  `(min(entry), max(exit))`, based on quantiles among observed event
  times. See Details.

Finally, for FMA the output additionally contains:

- samples:

  the samples generated from the normal distributions associated with
  each dose replicate.

- included.samples:

  the total number of samples included.

- included.replicates:

  indices of replicate exposures that were included to obtain results.
  Fits without a valid variance estimate (i.e., non-invertible Hessian
  or inverse that is not positive definite) or that reach the maximal
  number of iterations without convergence are filtered out and not used
  to obtain results.

The class `amerasfit` supports the methods
[`print`](https://ameras.sanderroberti.com/reference/print.md),
[`coef`](https://ameras.sanderroberti.com/reference/coef.md),
[`confint`](https://ameras.sanderroberti.com/reference/confint.md),
[`summary`](https://ameras.sanderroberti.com/reference/summary.md), and
[`traceplot`](https://ameras.sanderroberti.com/reference/traceplot.md).

## Details

Models are specified through formulas of the form
`Y~dose(dose_expression, model="ERR", deg=1, modifier=M1+M2)+X1+X2`.
Here `dose_expression` specifies the dose realization columns and is
parsed by `eval_select` from the tidyselect package. Useful examples are
`D1:D1000` if the doses are in a sequence of columns with sequential
names such as `D1`-`D1000`, and `all_of(dosevars)` where `dosevars` is a
vector with the names of all dose columns. Further, `model` specifies,
for non-Gaussian families, whether to use the exponential dose-response
model (`model="EXP"`), the linear-exponential model (`model="LINEXP"`)
or the linear ERR model (`model="ERR"`). Next, `deg` is used to specify
whether a quadratic dose term should (`deg=2`) or should not (`deg=1`)
be estimated for the exponential or linear ERR dose-response model. The
`modifier` term is optional and used to specify binary effect
modification variables. Note that interactions in the modifier term are
not allowed, e.g. `M1*M2`. When `deg`, `modifier`, and `model` are not
supplied, the defaults are `deg=1`, no effect modifiers, and
`model="ERR"`. Finally, `X1` and `X2` above represent optional
additional covariates, which can include factor variables and
interactions such as `X1*X2`. The matched set variable `setnr` required
for conditional logistic regression is specified on the right-hand side
of the formula through a term `strata(setnr)`, and an optional offset
variable `offset` for Poisson regression similarly through a term
`offset(offset)`. For proportional hazards regression, the left-hand
side of the formula should have the form `Surv(exit, status)` or
`Surv(entry, exit, status)`.

A transformation can be used to reparametrize parameters internally
(i.e., such that the likelihoods are evaluated at
`transform(parameters)`, where `parameters` are unconstrained), and
should be specified when fitting linear excess relative risk and
linear-exponential models to ensure nonnegative odds/risk/hazard. The
included function
[`transform1`](https://ameras.sanderroberti.com/reference/transform1.md)
applies an exponential transformation to the desired parameters, see
[`?transform1`](https://ameras.sanderroberti.com/reference/transform1.md).
When supplying a function to `transform`, this should be a function of
the full parameter vector, returning a full (transformed) parameter
vector. In particular, the full parameter vector contains parameters in
the following order: \\\alpha_0, \bm \alpha, \beta_1, \beta_2, \bm
\beta\_{m1}, \bm \beta\_{m2}, \sigma\\, where \\\bm \alpha\\,
\\\bm\beta\_{m1}\\ and \\\bm \beta\_{m2}\\ can be vectors, with lengths
matching \\\bm X\\ and \\\bm M\\, respectively. \\\sigma\\ is only
included for the linear model (Gaussian family), and no intercept is
included for the proportional hazards and conditional logistic models.
For the multinomial model, the full parameter vector is the
concatenation of \\Z-1\\ parameter vectors in the order as given above,
where \\Z\\ is the number of outcome categories, with the last category
chosen as the referent category. See
[`vignette("transformations", package="ameras")`](https://ameras.sanderroberti.com/articles/transformations.md)
for an example of how to specify a custom transformation function.

When no transformation is specified and the linear ERR model is used,
`transform1` is used for ERR parameters \\\beta_1\\ and \\\beta_2\\ by
default, with lower limits \\-1/max(D)\\ for \\\beta_1\\ in the linear
dose-response and \\(0,-1/max(D^2))\\ for \\(\beta_1,\beta_2)\\ in the
linear-quadratic dose-response, respectively. For the linear-exponential
model, a lower limit of 0 is used for \\\beta_1\\, and no transformation
is used for \\\beta_2\\. If effect modifiers `M` are specified, no
transformation is used for those parameters. When negative RRs are
obtained during optimization, an error will be generated and a different
transformation or bounds should be used. All output is returned in the
original parametrization. The Jacobian of the transformation
(`transform.jacobian`) is required when using a transformation. For
[`transform1`](https://ameras.sanderroberti.com/reference/transform1.md),
the Jacobian is given by
[`transform1.jacobian`](https://ameras.sanderroberti.com/reference/transform1jacobian.md).
No transformations are used in BMA, and FMA is applied on the parameters
using the parametrization as given in above with variances obtained
using the delta method with the provided Jacobian function.

For BMA, a prior distribution for exposure-response parameters can be
chosen when using linear or linear-exponential exposure-response model.
The options are normal, horshoe, and double exponential priors, and the
same priors truncated at 0 to yield positive values. In particular:

- Normal: \\\beta_j \sim N(0,1000)\\ for all exposure-response
  parameters \\\beta_j\\

- Horseshoe (shrinkage prior): \\\tau \sim \text{Cauchy}(0,1)^+;
  \lambda_j \sim \text{Cauchy}(0,1)^+; \beta_j \sim N(0, \tau^2
  \lambda_j^2)\\. Here \\\tau\\ is shared across all parameters

- Double exponential (shrinkage prior): \\\lambda_j \sim
  \text{Cauchy}(0,1)^+; \beta_j \sim
  \text{DoubleExponential}(0,\lambda_j)\\

For all other parameters, and when using the exponential
exposure-response model or the Gaussian outcome family, the prior is
\\N(0, 1000)\\. For the parameter \\\sigma\\ in the Gaussian family,
this prior is truncated at 0.

Because the proportional hazards model is not available in `nimble`,
`ameras` uses a piecewise constant baseline hazard for Bayesian model
averaging. The interval `min(entry), max(exit))` is divided into
`prophaz.numints.BMA` subintervals with cutpoints obtained as quantiles
of the distribution of event times among cases, and a baseline hazard
parameter is estimated for each subinterval.

## See also

[`confint`](https://ameras.sanderroberti.com/reference/confint.md) for
computing confidence intervals,
[`summary`](https://ameras.sanderroberti.com/reference/summary.md) for a
summary of the fitted model including confidence intervals if computed,
[`coef`](https://ameras.sanderroberti.com/reference/coef.md) for
extracting coefficients.

## References

Roberti, S., Kwon D., Wheeler W., Pfeiffer R. (in preparation). ameras:
An R Package to Analyze Multiple Exposure Realizations in Association
Studies

## Examples

``` r
# \donttest{
  data(data, package="ameras")
  ameras(Y.gaussian~dose(V1:V10, modifier=M1+M2)+X1+X2, data=data, family="gaussian") 
#> Fitting RC
#> Call:
#> ameras(formula = Y.gaussian ~ dose(V1:V10, modifier = M1 + M2) + 
#>     X1 + X2, data = data, family = "gaussian")
#> 
#> Number of rows: 3000
#> Number of dose replicates: 10
#> 
#> Total run time: 0.3 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.3
#> 
#> Estimated model parameters:
#> 
#>                  RC
#> (Intercept) -1.3796
#> X1           0.4966
#> X2          -0.5151
#> dose         0.9818
#> dose:M1      0.1739
#> dose:M2      0.5054
#> sigma        1.0661
# }
```
