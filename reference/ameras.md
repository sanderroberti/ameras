# Analyze Multiple Exposure Realizations

Fit regression models accounting for exposure uncertainty using multiple
Monte Carlo exposure realizations. Six outcome model families are
supported. The first is the Gaussian family for continuous outcomes,
\$\$Y_i \sim N(\mu_i, \sigma^2),\$\$ with \\\mu_i = \alpha_0 + \bm X_i^T
\bm \alpha +\beta_1 D_i+\beta_2 D_i^2 + \bm M_i^T \bm \beta\_{m1}D_i +
\bm M_i^T \bm \beta\_{m2}D_i^2\\. Here \\\bm X_i\\ are covariates,
\\D_i\\ is the exposure with measurement error, and \\\bm M_i\\ are
effect modifier design columns. The quadratic exposure terms and effect
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
  future.chunk.size.FMA=NULL,
  prophaz.numints.BMA=10, ERRprior.BMA="doubleexponential", nburnin.BMA=5000, 
  niter.BMA=20000, nchains.BMA=2, thin.BMA=10, included.realizations.BMA=NULL, 
  included.replicates.BMA=NULL, optim.method="Nelder-Mead", control=NULL, 
  keep.data=TRUE, na.action=getOption("na.action"), ... )
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
  indices of exposure realization vectors.

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

  function for parameter transformation during optimization (see
  Details).

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

- future.chunk.size.FMA:

  optional positive number controlling the average number of FMA dose
  realizations per future when future.apply is available. Parallel
  execution is controlled by the user's current
  [`future::plan()`](https://future.futureverse.org/reference/plan.html);
  with the default sequential plan, FMA is fit sequentially.

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

  number of MCMC burn-in iterations for BMA (default 5,000).

- niter.BMA:

  number of MCMC iterations per chain for BMA (default 20,000).

- nchains.BMA:

  number of MCMC chains for BMA (default 2).

- thin.BMA:

  thinning rate for BMA (default 10).

- included.realizations.BMA:

  indices of exposure realizations used in BMA (defaults to all
  realizations).

- included.replicates.BMA:

  **Deprecated**. Use `included.realizations.BMA` instead. Indices of
  exposure realizations used in BMA (defaults to all realizations).

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
  confidence intervals, but can also be supplied separately when
  `keep.data=FALSE`. See
  [`confint`](https://ameras.sanderroberti.com/reference/confint.md).

- na.action:

  function or function name controlling missing-value handling for model
  inputs after formula terms have been expanded. By default,
  `getOption("na.action")` is used, typically `na.omit`. Use `na.fail`
  to reject missing model inputs explicitly, or `na.pass` to skip row
  removal; missing model inputs left by `na.pass` will generally cause
  an error during input checks or model fitting. `na.exclude` is
  accepted for fitting, but residuals and diagnostic plots are returned
  for the fitted rows rather than padded back to the originally supplied
  row count.

- ...:

  other arguments, passed to functions such as `transform`.

## Value

The output is an object of class `amerasfit`. General components are
`call` (the function call to `ameras`), `formula` (the formula object
specifying the model), `num.rows` (the number of rows used for fitting
after applying `na.action` and any family-specific filtering),
`num.rows.original` (the number of rows supplied before applying
`na.action`), `na.action` (the omitted-row object returned by the
missing-value action, when applicable), `num.realizations` (the number
of dose realizations provided), `transform` (the used transformation
function, if applicable), `transform.jacobian` (the used Jacobian
function for the transformation, if applicable), `other.args` (any other
arguments passed to ...), `model` (a list containing the specified model
components parsed from the formula), `CI.computed` (logical, whether
confidence intervals have been attached by
[`confint`](https://ameras.sanderroberti.com/reference/confint.md)), and
`data` (either the data frame used for model fitting when
`keep.data=TRUE` or `NULL` otherwise).

For each method supplied to `methods`, the output contains a list with
components:

- coefficients:

  named vector of model coefficients.

- sd:

  named vector of standard deviations.

- vcov:

  covariance matrix for the full parameter vector. Based on the observed
  information (negative second derivative of log-likelihood) for RC,
  ERC, and MCML, and on the obtained samples for FMA and BMA.

- runtime:

  string with the total CPU runtime currently available for the method.
  After fitting this is the fitting time; after
  [`confint`](https://ameras.sanderroberti.com/reference/confint.md)
  this includes confidence interval computation time.

- timing:

  list with CPU and elapsed timings for fitting, confidence interval
  computation, and their total. CPU time is used for printed runtime
  summaries.

For RC, ERC, and MCML the following additional output is included:

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

  MCMC posterior samples, as obtained from `nimble`. With
  `nchains.BMA > 1`, this is a list with `nchains.BMA` components, each
  a named matrix with the samples from one chain in its rows and columns
  corresponding to model parameters. With `nchains.BMA = 1`, this is a
  single named matrix.

- Rhat:

  data frame with two columns, `Rhat` and `n.eff`. The first column
  contains the Gelman-Rubin statistics \\\hat R \geq 1\\ that can be
  used to assess convergence of MCMC chains. A value of 1 indicates good
  convergence and values \\\>1.05\\ indicate poor convergence. The
  effective sample size `n.eff` is a measure of how many independent
  samples the auto-correlated MCMC samples correspond to. A low
  effective sample size indicates high correlations and/or poor mixing.
  For single-chain BMA fits, these diagnostics cannot be estimated and
  are returned as `NA`.

- included.realizations:

  indices of realization exposures that were included to obtain the
  results.

- prophaz.timepoints:

  for `family="prophaz"`, time points defining the intervals on which
  the estimated baseline hazard is constant; these are
  `prophaz.numints.BMA + 1` time points covering the interval
  `(min(entry), max(exit))`, based on quantiles among observed event
  times. See Details.

Finally, for FMA the output additionally contains:

- samples:

  data frame with samples generated from the normal distributions
  associated with parameters estimated for each dose realization.

- weights:

  vector of weights corresponding to the models fit to each included
  realization

- included.samples:

  the total number of samples included.

- included.realizations:

  indices of realization exposures that were included to obtain results.
  Fits without a valid variance estimate (i.e., non-invertible Hessian
  or inverse that is not positive definite) or that reach the maximal
  number of iterations without convergence are filtered out and not used
  to obtain results.

The class `amerasfit` supports the methods
[`print`](https://ameras.sanderroberti.com/reference/print.md),
[`coef`](https://ameras.sanderroberti.com/reference/coef.md),
[`confint`](https://ameras.sanderroberti.com/reference/confint.md),
[`summary`](https://ameras.sanderroberti.com/reference/summary.md),
[`summary_table`](https://ameras.sanderroberti.com/reference/summary.md),
[`vcov`](https://ameras.sanderroberti.com/reference/vcov.md),
[`residuals`](https://ameras.sanderroberti.com/reference/residuals.md),
[`plot`](https://ameras.sanderroberti.com/reference/plot.md),
[`traceplot`](https://ameras.sanderroberti.com/reference/traceplot.md),
[`convergence`](https://ameras.sanderroberti.com/reference/convergence.md),
[`Rhat`](https://ameras.sanderroberti.com/reference/Rhat.md), and
[`included_realizations`](https://ameras.sanderroberti.com/reference/included_realizations.md).

## Details

Missing values in model inputs are handled by `na.action` after formula
terms have been expanded. Thus missing values in raw variables and in
derived covariate columns, such as factor contrasts, interactions, and
spline bases, are handled consistently before model fitting. The
omitted-row action is stored on the fitted object. When
`keep.data=FALSE` and data are supplied later to methods such as
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), the same
missing-value policy is reapplied. Unlike
[`lm`](https://rdrr.io/r/stats/lm.html) and
[`glm`](https://rdrr.io/r/stats/glm.html), `na.exclude` does not pad
residuals or diagnostic plots back to the originally supplied row count;
those outputs are computed on the fitted rows.

Models are specified through formulas of the form
`Y~dose(dose_expression, model="ERR", deg=1, modifier=~M1+M2)+X1+X2`.
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
`modifier` term is optional and used to specify effect modification
variables through a formula. Use `modifier=~M` for
reference-plus-contrast coding, and `modifier=~0+M` or `modifier=~M-1`
for subgroup-specific dose coefficients. Formula modifiers may be binary
numeric/logical variables or factors; factor modifiers use the existing
factor level order, with the first level as the reference. Interactions,
transformed modifier terms, continuous numeric modifiers, character
modifiers, and no-intercept coding with multiple modifier variables are
not currently supported. The older syntax `modifier=M1+M2` is deprecated
and remains binary-only. When `deg`, `modifier`, and `model` are not
supplied, the defaults are `deg=1`, no effect modifiers, and
`model="ERR"`.

Additional right-hand-side covariates are expanded using
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html), so standard
formula features such as factor contrasts, interactions,
[`I()`](https://rdrr.io/r/base/AsIs.html) terms, and basis functions can
be used. For example, formulas can include terms such as `X1*X2`,
`I(age^2)`, `splines::ns(age, df=3)`, or `splines::bs(age, df=4)`. These
terms must preserve one value per row of the input data; terms that drop
rows, such as `stats::na.omit(X)`, are not supported. When
`keep.data=FALSE`, later calls that use supplied data apply the same
expanded covariate columns as the original fit. If covariates appear
poorly scaled or ill-conditioned, ameras emits a diagnostic warning;
centering or scaling continuous covariates such as calendar year can
improve numerical optimization.

The matched set variable `setnr` required for conditional logistic
regression is specified on the right-hand side of the formula through a
term `strata(setnr)`, and an optional offset variable `offset` for
Poisson regression similarly through a term `offset(offset)`. For
conditional logistic regression, matched sets of size 1 and matched sets
with no cases are excluded, and matched sets with more than one case are
not currently supported. For proportional hazards regression, the
left-hand side of the formula should have the form `Surv(exit, status)`
or `Surv(entry, exit, status)`. When entry times are supplied, rows with
`entry == exit` are excluded because they have zero follow-up time under
the `entry < t <= exit` risk-set convention.

A transformation can be used to reparametrize parameters for
optimization, such that the likelihoods are evaluated at
`transform(parameters)`, where `parameters` are unconstrained. This
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
is used for \\\beta_2\\. If reference-plus-contrast effect modifiers are
specified, no transformation is used for modifier contrast parameters.
With subgroup-specific modifier coding, the default transformation is
applied to the subgroup dose-effect parameters. For ERR and
linear-exponential models with effect modifiers, subgroup-specific
coding can therefore be more numerically robust than
reference-plus-contrast coding, because bounds are applied directly to
each subgroup dose-effect parameter. When negative RRs are obtained
during optimization, an error will be generated and a different
transformation or bounds should be used. All output is returned in the
original parametrization. The Jacobian of the transformation
(`transform.jacobian`) is required when using a transformation. For
[`transform1`](https://ameras.sanderroberti.com/reference/transform1.md),
the Jacobian is given by
[`transform1.jacobian`](https://ameras.sanderroberti.com/reference/transform1jacobian.md).
No transformations are used in BMA, and FMA is applied on the parameters
using the parametrization as given above with variances obtained using
the delta method with the provided Jacobian function.

For BMA, a prior distribution for exposure-response parameters can be
chosen when using linear or linear-exponential exposure-response model.
The options are normal, horseshoe, and double exponential priors, and
the same priors truncated at 0 to yield positive values. In particular:

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
this prior is truncated at 0. When subgroup-specific modifier coding is
used in BMA, priors and MCMC sampling remain on the internal
reference-plus-contrast scale; reported coefficients, variance matrices,
samples, diagnostics, and sample-based intervals are converted to the
subgroup-specific scale.

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
extracting coefficients,
[`vcov`](https://ameras.sanderroberti.com/reference/vcov.md) for
extracting variance-covariance matrices, and
[`included_realizations`](https://ameras.sanderroberti.com/reference/included_realizations.md)
for inspecting which realizations contributed to model averaging.

## References

Roberti, S., Kwon D., Wheeler W., Pfeiffer R. (in preparation). ameras:
An R Package to Analyze Multiple Exposure Realizations in Association
Studies

## Examples

``` r
# \donttest{
  data(data, package="ameras")
  ameras(Y.gaussian~dose(V1:V10, modifier=~M1+M2)+X1+X2, data=data, family="gaussian")
#> Fitting RC
#> Call:
#> ameras(formula = Y.gaussian ~ dose(V1:V10, modifier = ~M1 + M2) + 
#>     X1 + X2, data = data, family = "gaussian")
#> 
#> Number of rows: 3000
#> Number of dose realizations: 10
#> 
#> Total CPU runtime: 0.5 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method   Fit  CI Total
#>      RC 0.472 0.0 0.472
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
#> 
# }
```
