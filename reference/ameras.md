# Analyze Multiple Exposure Realizations

Fit regression models that account for exposure uncertainty using
multiple exposure realizations. Supported outcome models are Gaussian,
binomial, Poisson, proportional hazards, multinomial, and conditional
logistic regression. Estimation methods include regression calibration
(RC), extended regression calibration (ERC), Monte Carlo maximum
likelihood (MCML), frequentist model averaging (FMA), and Bayesian model
averaging (BMA).

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
  realizations per future. When `NULL` (the default), approximately one
  future per available worker is used. Parallel execution is controlled
  by the user's current
  [`future::plan()`](https://future.futureverse.org/reference/plan.html);
  the default plan is sequential.

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
  is normally run to refine the fit. If this automatic BFGS refinement
  fails because its finite-difference evaluations are non-finite, the
  Nelder-Mead solution is retained and checked using the usual Hessian
  and convergence diagnostics. A warning is emitted if the selected
  final optimizer result has a nonzero convergence code.

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

## Formula Syntax

Models are specified through formulas of the form
`Y ~ dose(dose_expression, model="ERR", deg=1, modifier=~M1+M2) + X1 + X2`.
Exactly one `dose()` term is required. Its unnamed `dose_expression`
argument selects the exposure realization columns using
[`eval_select`](https://tidyselect.r-lib.org/reference/eval_select.html)
syntax. For example, use `D1:D1000` for sequentially named columns or
`all_of(dosevars)` for a character vector of column names.

The recognized named arguments inside `dose()` are `model`, `deg`, and
`modifier`; their names must be written exactly. For non-Gaussian
families, `model` selects the exponential relative risk model (`"EXP"`),
linear-exponential model (`"LINEXP"`), or linear excess relative risk
model (`"ERR"`, the default). For exponential and linear ERR models,
`deg=1` fits a linear exposure effect and `deg=2` adds a quadratic term.

The optional `modifier` argument uses a one-sided formula. Use
`modifier=~M` for reference-plus-contrast coding and `modifier=~0+M` or
`modifier=~M-1` for subgroup-specific dose coefficients. Modifier
variables may be binary numeric, logical, or factor with the first level
used as the reference category. Interactions, transformed modifier
terms, continuous numeric modifiers, character modifiers, and
no-intercept coding with multiple modifier variables are not supported.

The special terms are parsed more narrowly than ordinary model-formula
terms. For non-survival models, the response must be a single,
untransformed data column. The argument to `strata()` for conditional
logistic regression and the argument to
[`offset()`](https://rdrr.io/r/stats/offset.html) for Poisson regression
must each be one untransformed data column. Proportional hazards models
require `Surv(exit, status)` or `Surv(entry, exit, status)`, with each
argument naming a data column. Expressions such as
`offset(log(person_time))` or `strata(interaction(center, set))` are
therefore not supported; create the derived column in `data` first.

Additional right-hand-side covariates are expanded with
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) and may use
standard formula features such as factor contrasts, interactions,
[`I()`](https://rdrr.io/r/base/AsIs.html) terms, and basis functions.
For example, `X1*X2`, `I(age^2)`, `splines::ns(age, df=3)`, and
`splines::bs(age, df=4)` are supported. These terms must preserve one
value per input row; terms that drop rows internally, such as
`stats::na.omit(X)`, are not supported. Removing the overall intercept
with `-1` or `+0` is also not supported.

## Supported Outcome Models

For the Gaussian family, \$\$Y_i \sim N(\mu_i, \sigma^2),\$\$ with
\$\$\mu_i = \alpha_0 + \bm X_i^T \bm \alpha + \beta_1 D_i + \beta_2
D_i^2 + \bm M_i^T \bm \beta\_{m1}D_i + \bm M_i^T \bm
\beta\_{m2}D_i^2.\$\$ Here \\\bm X_i\\ are covariates, \\D_i\\ is the
exposure with measurement error, and \\\bm M_i\\ are effect modifier
design columns. The quadratic exposure terms and effect modification are
optional.

For non-Gaussian families, the supported main-exposure relative risk
models are the exponential model \$\$RR_i=\exp(\beta_1 D_i+\beta_2
D_i^2+ \bm M_i^T \bm \beta\_{m1}D_i + \bm M_i^T \bm
\beta\_{m2}D_i^2),\$\$ the linear excess relative risk (ERR) model
\$\$RR_i=1+\beta_1 D_i+\beta_2 D_i^2+ \bm M_i^T \bm \beta\_{m1}D_i + \bm
M_i^T \bm \beta\_{m2}D_i^2,\$\$ and the linear-exponential model
\$\$RR_i=1+(\beta_1+\bm M_i^T \bm \beta\_{m1})D_i \exp\\(\beta_2+\bm
M_i^T \bm \beta\_{m2})D_i\\.\$\$ Quadratic and effect modification terms
are optional where applicable.

For binary outcomes, logistic regression uses
\$\$p_i/(1-p_i)=RR_i\exp(\alpha_0+\bm X_i^T \bm \alpha).\$\$

For counts, Poisson regression uses \$\$Y_i \sim \text{Poisson}(\mu_i),
\qquad \mu_i=RR_i\exp(\alpha_0+\bm X_i^T \bm \alpha)\times
\text{offset}\_i,\$\$ with an optional offset.

For time-to-event outcomes, proportional hazards regression uses
\$\$h(t)=h_0(t)RR_i\exp(\bm X_i^T \bm \alpha),\$\$ where \\h_0\\ is the
baseline hazard.

For a categorical outcome with \\Z\>2\\ categories, multinomial logistic
regression uses the last category as the referent (\\\alpha\_{0,Z}=\bm
\alpha_Z=\beta\_{1,Z}=\beta\_{2,Z}= \bm \beta\_{m1,Z}=\bm
\beta\_{m2,Z}=0\\): \$\$P(Y_i=z)=\frac{RR_i\exp(\alpha\_{0,z}+\bm X_i^T
\bm \alpha_z)} {1+\sum\_{s=1}^{Z-1}RR_i\exp(\alpha\_{0,s}+\bm X_i^T \bm
\alpha_s)}.\$\$

For matched case-control data, conditional logistic regression uses
\$\$P\left(Y_i=1, Y_k=0\\ \forall k\ne i\\ \bigg\|\\ \sum\_{i\in\mathcal
R}Y_i=1\right)= \frac{RR_i\exp(\bm X_i^T \bm \alpha)}
{\sum\_{k\in\mathcal R}RR_k\exp(\bm X_k^T \bm \alpha)},\$\$ where
\\\mathcal R\\ is the matched set containing individual \\i\\.

The estimation methods are described by Carroll et al. (2006)
[doi:10.1201/9781420010138](https://doi.org/10.1201/9781420010138) for
RC, Little et al. (2023)
[doi:10.1038/s41598-023-42283-y](https://doi.org/10.1038/s41598-023-42283-y)
for ERC, Stayner et al. (2007)
[doi:10.1667/RR0677.1](https://doi.org/10.1667/RR0677.1) for MCML, Kwon
et al. (2023)
[doi:10.1371/journal.pone.0290498](https://doi.org/10.1371/journal.pone.0290498)
for FMA, and Kwon et al. (2016)
[doi:10.1002/sim.6635](https://doi.org/10.1002/sim.6635) for BMA.

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
components parsed from the formula, including `model$data`, which is
either the data frame used for model fitting when `keep.data=TRUE` or
`NULL` otherwise), and `CI.computed` (logical, whether confidence
intervals have been attached by
[`confint`](https://ameras.sanderroberti.com/reference/confint.md)).

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

- timing:

  list with CPU and elapsed timings for fitting, confidence interval
  computation, and their total. CPU time is used for printed runtime
  summaries.

For RC, ERC, and MCML the following additional output is included:

- optim:

  a list object with results returned by optim. Components are `par`
  (raw parameters before applying a transformation if applicable),
  `hessian` (Hessian matrix for `par`), `convergence` (convergence code
  with 0 indicating convergence and a nonzero value indicating that the
  optimizer did not converge; code 1 indicates that the maximal number
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

When `keep.data=FALSE`, later calls that use supplied data apply the
same expanded covariate columns as the original fit. If covariates
appear poorly scaled or ill-conditioned, ameras emits a diagnostic
warning; centering or scaling continuous covariates such as calendar
year can improve numerical optimization.

For conditional logistic regression, matched sets of size 1 and matched
sets with no cases are excluded, and matched sets with more than one
case are not currently supported. For proportional hazards regression
with entry times, rows with `entry == exit` are excluded because they
have zero follow-up time under the `entry < t <= exit` risk-set
convention.

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

For extended examples, see
[`vignette("modelfitting", package="ameras")`](https://ameras.sanderroberti.com/articles/modelfitting.md)
for the general workflow,
[`vignette("relativeriskmodels", package="ameras")`](https://ameras.sanderroberti.com/articles/relativeriskmodels.md)
for relative risk models,
[`vignette("effectmodification", package="ameras")`](https://ameras.sanderroberti.com/articles/effectmodification.md)
for effect modification, and
[`vignette("confidenceintervals", package="ameras")`](https://ameras.sanderroberti.com/articles/confidenceintervals.md)
for confidence and credible intervals.

## References

Roberti, S., Kwon D., Wheeler W., Pfeiffer R. (in preparation). ameras:
An R Package to Analyze Multiple Exposure Realizations in Association
Studies

## Examples

``` r
## Fit a Gaussian regression-calibration model using ten exposure realizations
example_data <- ameras_simulated[1:500, ]
fit <- ameras(
  Y.gaussian ~ dose(V1:V10) + X1 + X2,
  data = example_data,
  family = "gaussian",
  methods = "RC"
)
#> Fitting RC

## Extract estimates
coef(fit)
#>                     RC
#> (Intercept) -1.4325161
#> X1           0.5776954
#> X2          -0.5216418
#> dose         1.1596141
#> sigma        1.1238400

## Add Wald confidence intervals and display a detailed summary
fit <- confint(fit, type = "wald.orig", print = FALSE)
summary(fit)
#> Call:
#> ameras(formula = Y.gaussian ~ dose(V1:V10) + X1 + X2, data = example_data, 
#>     family = "gaussian", methods = "RC")
#> 
#> Rows: 500
#> 
#> Total CPU runtime: 0.1 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method   Fit  CI Total
#>      RC 0.081 0.0 0.081
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper
#>      RC (Intercept)  -1.4325 0.09084   -1.611  -1.2545
#>      RC          X1   0.5777 0.10087    0.380   0.7754
#>      RC          X2  -0.5216 0.12672   -0.770  -0.2733
#>      RC        dose   1.1596 0.05252    1.057   1.2625
#>      RC       sigma   1.1238 0.03554    1.054   1.1935
#> 
```
