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
ameras(data, family="gaussian", Y, dosevars, M=NULL, X=NULL, offset=NULL, entry=NULL, 
  exit=NULL, setnr=NULL, methods="RC", deg=1, doseRRmod="ERR", transform=NULL,
  transform.jacobian=NULL, inpar=NULL, CI=c("proflik","percentile"),
  params.profCI="dose", maxit.profCI=20, tol.profCI=1e-2, loglim=1e-30, MFMA=100000, 
  prophaz.numints.BMA=10, ERRprior.BMA="doubleexponential", nburnin.BMA=5000, 
  niter.BMA=20000, nchains.BMA=2, thin.BMA=10, included.replicates.BMA=1:length(dosevars), 
  optim.method="Nelder-Mead", control=NULL, ... )
```

## Arguments

- data:

  input data frame.

- family:

  outcome model family: `"gaussian"`, `"binomial"`, `"poisson"`,
  `"prophaz"`, `"multinomial"` or `"clogit"` (default `"gaussian"`).

- Y:

  name or column index of the outcome variable for linear, binomial,
  Poisson, multinomial and conditional logistic models, or event
  indicator variable for the proportional hazards model.

- dosevars:

  names or column indices of exposure replicate vectors.

- M:

  names or column indices of binary effect modifying variables
  (optional).

- X:

  names or column indices of other covariates (optional).

- offset:

  name or column index of offset variable for Poisson regression
  (optional).

- entry:

  name or column index of left truncation time variable for proportional
  hazards regression (optional).

- exit:

  name or column index of exit time variable, required when
  `family=prophaz`.

- setnr:

  name or column index of integer-valued matched set variable, required
  when `family="clogit"`.

- methods:

  character vector of one or multiple methods to apply. Options: `"RC"`,
  `"ERC"`, `"MCML"`, `"FMA"`, `"BMA"` (default `"RC"`).

- deg:

  for `doseRRmod="ERR"` and `doseRRmod="EXP"`, whether to fit a linear
  (`deg=1`) or linear-quadratic (`deg=2`) dose-response model (default
  linear).

- doseRRmod:

  the functional form of the dose-response relationship; options are
  exponential RR (`"EXP"`), linear ERR (`"ERR"`), or linear-exponential
  RR (`"LINEXP"`) (default `"ERR"`).

- transform:

  function for internal parameter transformation (see Details).

- transform.jacobian:

  Jacobian of the transformation function (see Details).

- inpar:

  vector of initial values for log-likelihood optimization (optional).

- CI:

  method for calculation of 95% confidence or credible intervals (see
  Details). For RC, ERC, and MCML, options are `"wald.orig"`,
  `"wald.transformed"`, `"proflik"` (default `"proflik"`). For FMA and
  BMA, options are `"percentile"` and `"hpd"` (default `"percentile"`).
  If `methods` contains at least one of RC, ERC, and MCML and at least
  one of FMA and BMA, `CI` must be length 2 and specify one method for
  RC, ERC, and MCML, and one for FMA and BMA (see Details).

- params.profCI:

  when `CI="proflik"`, whether to obtain profile-likelihood CIs for all
  parameters (`"all"`) or only dose-related parameters (`"dose"`,
  default).

- maxit.profCI:

  maximum iterations for determining profile-likelihood CIs; passed to
  `uniroot` (default 20).

- tol.profCI:

  tolerance for determining profile-likelihood CIs; passed to `uniroot`
  (default `1e-2`).

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

  indices of exposure replicates used in BMA (default \\
  `1:length(dosevars)`).

- optim.method:

  method used for optimization by `optim`. Options are `"Nelder-Mead"`
  and `"BFGS"`. When using Nelder-Mead, a second optimization with BFGS
  is run to ensure an optimal fit.

- control:

  control list passed to `optim` (default `list(reltol=1e-10)`).

- ...:

  other arguments, passed to functions such as `transform`.

## Value

The output is an object of class `amerasfit` with a component `call` and
a component for every method supplied to `methods`. For each method, the
output is a list containing

- coefficients:

  named vector of model coefficients.

- sd:

  named vector of standard deviations.

- CI:

  data frame with columns `lower` and `upper` giving 95% confidence
  bounds or credible interval bounds. When the CI method is `"proflik"`,
  the data frame also has columns `pval.lower` and `pval.upper`
  (p-values to verify convergence of the root finder) and `iter.lower`
  and `iter.upper` (number of iterations used by `uniroot`).

- runtime:

  string with the runtime in seconds.

For RC, ERC, and MCML the following additional output is included:

- vcov:

  covariance matrix for the full parameter vector.

- convergence.optim:

  convergence code as returned by `optim`, with 0 indicating convergence
  and 1 indicating that the maximal number of iterations was reached.

- counts.optim:

  number of function evaluations used in the model fit returned by
  `optim`.

- loglik:

  log-likelihood value at the optimum.

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

- included.samples:

  the total number of samples included.

- included.replicates:

  indices of replicate exposures that were included to obtain results.
  Fits without a valid variance estimate (i.e., non-invertible Hessian
  or inverse that is not positive definite) or that reach the maximal
  number of iterations without convergence are filtered out and not used
  to obtain results.

The class `amerasfit` supports the methods `coef`, `summary`, and
[traceplot](https://ameras.sanderroberti.com/reference/traceplot.md).

## Details

A transformation can be used to reparametrize parameters internally
(i.e., such that the likelihoods are evaluated at
`transform(parameters)`, where `parameters` are unconstrained), and
should be specified when fitting linear excess relative risk and
linear-exponential models to ensure nonnegative odds/risk/hazard. The
included function
[transform1](https://ameras.sanderroberti.com/reference/transform1.md)
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
[transform1](https://ameras.sanderroberti.com/reference/transform1.md),
the Jacobian is given by
[transform1.jacobian](https://ameras.sanderroberti.com/reference/transform1jacobian.md).
No transformations are used in BMA, and FMA is applied on the parameters
using the parametrization as given in above with variances obtained
using the delta method with the provided Jacobian function.

Multiple options for confidence intervals are provided. For (extended)
regression calibration and Monte Carlo maximum likelihood, Wald and
profile likelihood intervals can be obtained. When a parameter
transformation \\\bm\theta = h(\bm\eta)\\ is used,
`CI="wald.transformed"` yields the CI \\h(\bm\eta \pm 1.96 \bm V)\\ with
\\\bm V\\ the vector of standard deviations estimated using the inverse
Hessian matrix, and `CI="wald.orig"` uses the delta method to obtain the
CI \\h(\bm\eta)\pm 1.96 \bm V\_\*\\ where \\\bm V\_\*\\ is the vector of
standard deviations estimated using \\J H^{-1} J^T\\ with \\J\\ the
Jacobian of the transformation and \\H\\ is the Hessian. When no
transformation is used, `CI="wald.orig"` should be used. The third
option is `proflik`, which uses the profile likelihood to compute
confidence bounds. For FMA and BMA, the options for confidence/credible
intervals are `CI="percentile"` which uses 2.5% and 97.5% percentiles,
and `CI="hpd"` which computes highest posterior density intervals using
`HPDinterval` from the `coda` package, both using the FMA samples or
Bayesian posterior samples.

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

## References

Roberti, S., Kwon D., Wheeler W., Pfeiffer R. (in preparation). ameras:
An R Package to Analyze Multiple Exposure Realizations in Association
Studies

## Examples

``` r
# \donttest{
  data(data, package="ameras")
  dosevars <- paste0("V", 1:10)
  ameras(data=data, family="gaussian", Y="Y.gaussian", dosevars=dosevars, 
  M=c("M1", "M2"), X=c("X1","X2")) 
#> Fitting RC
#> Obtaining profile likelihood CI for dose
#> Warning: P-value for dose upper bound more than 0.005 away from 0.05, reducing tol.profCI and/or increasing maxit.profCI is recommended
#> Warning: P-value for dose lower bound more than 0.005 away from 0.05, reducing tol.profCI and/or increasing maxit.profCI is recommended
#> Obtaining profile likelihood CI for dose:M1
#> Warning: P-value for dose:M1 upper bound more than 0.005 away from 0.05, reducing tol.profCI and/or increasing maxit.profCI is recommended
#> Warning: P-value for dose:M1 lower bound more than 0.005 away from 0.05, reducing tol.profCI and/or increasing maxit.profCI is recommended
#> Obtaining profile likelihood CI for dose:M2
#> Warning: P-value for dose:M2 upper bound more than 0.005 away from 0.05, reducing tol.profCI and/or increasing maxit.profCI is recommended
#> Warning: P-value for dose:M2 lower bound more than 0.005 away from 0.05, reducing tol.profCI and/or increasing maxit.profCI is recommended
#> Call:
#> ameras(data = data, family = "gaussian", Y = "Y.gaussian", dosevars = dosevars, 
#>     M = c("M1", "M2"), X = c("X1", "X2"))
#> 
#> Number of individuals: 3000
#> Number of dose replicates: 10
#> 
#> Total run time: 6.3 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     6.3
# }
```
