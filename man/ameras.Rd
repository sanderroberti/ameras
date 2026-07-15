
\name{ameras}
\alias{ameras}
\title{Analyze Multiple Exposure Realizations}

\description{
Fit regression models accounting for exposure uncertainty using multiple Monte Carlo exposure realizations. Six outcome model families are supported. The first is the Gaussian family for continuous outcomes, \deqn{Y_i \sim N(\mu_i, \sigma^2),} with \eqn{\mu_i = \alpha_0 + \bm X_i^T \bm \alpha +\beta_1 D_i+\beta_2 D_i^2 + \bm M_i^T \bm \beta_{m1}D_i + \bm M_i^T \bm \beta_{m2}D_i^2}. Here \eqn{\bm X_i} are covariates, \eqn{D_i} is the exposure with measurement error, and \eqn{\bm M_i} are effect modifier design columns. The quadratic exposure terms and effect modification are optional.

For non-Gaussian families, three relative risk models for the main exposure are supported, the usual exponential \eqn{RR_i=\exp(\beta_1 D_i+\beta_2 D_i^2+ \bm M_i^T \bm \beta_{m1}D_i + \bm M_i^T \bm \beta_{m2} D_i^2)} and the linear excess relative risk (ERR) model \eqn{RR_i= 1+\beta_1 D_i+\beta_2 D_i^2 + \bm M_i^T \bm \beta_{m1}D_i + \bm M_i^T \bm \beta_{m2}D_i^2}, where the quadratic and effect modification terms are optional. Finally, the linear-exponential relative risk model \eqn{RR_i= 1+(\beta_1 + \bm{M}_i^T \bm{\beta}_{m1}) D_i \exp\{(\beta_2+  \bm{M}_i^T \bm{\beta}_{m2})D_i\}} is supported.

The second supported family is logistic regression for binary outcomes, with probabilities \deqn{p_i/(1-p_i)=RR_i\exp(\alpha_0+\bm X_i^T \bm \alpha).}

Third is Poisson regression for counts, \deqn{Y_i \sim \text{Poisson}(\mu_i),} where \eqn{\mu_i=RR_i \exp(\alpha_0 +\bm X_i^T \bm  \alpha)\times \text{offset}_i} with optional offset.

Fourth is proportional hazards regression for time-to-event data, with hazard function \deqn{h(t) = h_0(t)RR_i\exp(\bm X_i^T \bm \alpha),} with \eqn{h_0} the baseline hazard.

Fifth is multinomial logistic regression for a categorical outcome with \eqn{Z>2} outcome categories, with the last category as the referent category (i.e., \eqn{\alpha_{0,Z}=\bm \alpha_{Z}=\beta_{1,Z}=\beta_{2,Z}=\bm \beta_{m1,Z} = \bm \beta_{m2,Z}=0}): \deqn{P(Y_i=z)=RR_i\exp(\alpha_{0,z}+\bm X_i^T \bm \alpha_{z})/\{1+\sum_{s=1}^{Z-1} RR_i\exp(\alpha_{0,s}+\bm X_i^T \bm \alpha_{s})\}}

Sixth is conditional logistic regression for matched case control data, for which \deqn{P\left(Y_i = 1, Y_k = 0 \forall k \neq i \bigg| \sum_{i \in \mathcal{R}} Y_i = 1\right) = RR_i\exp(\bm X_i^T \bm \alpha)/\{\sum_{k \in \mathcal{R}}RR_k\exp(\bm X_k^T \bm \alpha)\},} where \eqn{\mathcal{R}} is the matched set corresponding to individual \eqn{i}.

Methods include regression calibration (Carroll et al. 2006 \doi{10.1201/9781420010138}), extended regression calibration (Little et al. 2023 \doi{10.1038/s41598-023-42283-y}), Monte Carlo maximum likelihood (Stayner et al. 2007 \doi{10.1667/RR0677.1}), frequentist model averaging (Kwon et al. 2023 \doi{10.1371/journal.pone.0290498}), and Bayesian model averaging (Kwon et al. 2016 \doi{10.1002/sim.6635}).

}

\usage{
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
}
\arguments{
\item{formula}{an object of class \code{"formula"} containing the model specification. See Details.}
\item{data}{input data frame.}
\item{family}{outcome model family: \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, \code{"prophaz"}, \code{"multinomial"} or \code{"clogit"} (default \code{"gaussian"}).}
\item{methods}{character vector of one or multiple methods to apply. Options: \code{"RC"}, \code{"ERC"}, \code{"MCML"}, \code{"FMA"}, \code{"BMA"} (default \code{"RC"}).}
\item{Y}{\strong{Deprecated}. Use the formula interface instead. Name or column index of the outcome variable for linear, binomial, Poisson, multinomial and conditional logistic models, or event indicator variable for the proportional hazards model.}
\item{dosevars}{\strong{Deprecated}. Use the formula interface instead. Names or column indices of exposure realization vectors.}
\item{doseRRmod}{\strong{Deprecated}. Use the formula interface instead. The functional form of the dose-response relationship; options are exponential RR (\code{"EXP"}), linear ERR (\code{"ERR"}), or linear-exponential RR (\code{"LINEXP"}) (default \code{"ERR"}).}
\item{deg}{\strong{Deprecated}. Use the formula interface instead. For \code{doseRRmod="ERR"} and \code{doseRRmod="EXP"}, whether to fit a linear (\code{deg=1}) or linear-quadratic (\code{deg=2}) dose-response model (default linear).}
\item{M}{\strong{Deprecated}. Use the formula interface instead. Names or column indices of binary effect modifying variables (optional).}
\item{X}{\strong{Deprecated}. Use the formula interface instead. Names or column indices of other covariates (optional).}
\item{offset}{\strong{Deprecated}. Use the formula interface instead. Name or column index of offset variable for Poisson regression (optional).}
\item{entry}{\strong{Deprecated}. Use the formula interface instead. Name or column index of left truncation time variable for proportional hazards regression (optional).}
\item{exit}{\strong{Deprecated}. Use the formula interface instead. Name or column index of exit time variable, required when \code{family=prophaz}.}
\item{setnr}{\strong{Deprecated}. Use the formula interface instead. Name or column index of integer-valued matched set variable, required when \code{family="clogit"}.}
\item{CI}{\strong{Deprecated}. Use confint() to compute confidence intervals instead. Method for calculation of 95\% confidence or credible intervals (see Details).
For RC, ERC, and MCML, options are \code{"wald.orig"}, \code{"wald.transformed"}, \code{"proflik"} (default \code{"proflik"}).
For FMA and BMA, options are \code{"percentile"} and \code{"hpd"} (default \code{"percentile"}).
If \code{methods} contains at least one of RC, ERC, and MCML and at least one of FMA and BMA, \code{CI} must be length 2 and specify one method for RC, ERC, and MCML, and one for FMA and BMA (see Details).}
\item{params.profCI}{\strong{Deprecated}. Use confint() to compute confidence intervals instead. When \code{CI="proflik"}, whether to obtain profile-likelihood CIs for all parameters (\code{"all"}) or only dose-related parameters (\code{"dose"}, default).}
\item{maxit.profCI}{\strong{Deprecated}. Use confint() to compute confidence intervals instead. Maximum iterations for determining profile-likelihood CIs; passed to \code{uniroot} (default 20).}
\item{tol.profCI}{\strong{Deprecated}. Use confint() to compute confidence intervals instead. Tolerance for determining profile-likelihood CIs; passed to \code{uniroot} (default \code{1e-2}).}
\item{transform}{function for parameter transformation during optimization
(see Details).}
\item{transform.jacobian}{Jacobian of the transformation function (see Details).}
\item{inpar}{vector of initial values for log-likelihood optimization (optional).}
\item{loglim}{parameter used in likelihood computations to avoid taking the log of very small or negative numbers via \code{log(max(x, loglim))} (default \code{1e-30}).}
\item{MFMA}{number of samples for \code{"FMA"} to compute estimates and CIs (default 100,000).}
\item{future.chunk.size.FMA}{optional positive number controlling the average number of FMA dose realizations per future when \pkg{future.apply} is available. Parallel execution is controlled by the user's current \code{future::plan()}; with the default sequential plan, FMA is fit sequentially.}
\item{prophaz.numints.BMA}{for \code{methods="BMA"} with \code{family="prophaz"}, the number of subintervals with constant baseline hazard (default 10). Cut points are determined based on quantiles of the event time distribution among cases.}
\item{ERRprior.BMA}{prior for dose-related parameters when \code{doseRRmod="ERR"} or \code{"LINEXP"} and \code{methods="BMA"}.
Options: \code{"truncated_normal"}, \code{"truncated_horseshoe"}, \code{"truncated_doubleexponential"}, \code{"normal"}, \code{"horseshoe"}, \code{"doubleexponential"}, see Details (default \code{"doubleexponential"}).}
\item{nburnin.BMA}{number of MCMC burn-in iterations for BMA (default 5,000).}
\item{niter.BMA}{number of MCMC iterations per chain for BMA (default 20,000).}
\item{nchains.BMA}{number of MCMC chains for BMA (default 2).}
\item{thin.BMA}{thinning rate for BMA (default 10).}
\item{included.realizations.BMA}{indices of exposure realizations used in BMA (defaults to all realizations).}
\item{included.replicates.BMA}{\strong{Deprecated}. Use \code{included.realizations.BMA} instead. Indices of exposure realizations used in BMA (defaults to all realizations).}
\item{optim.method}{method used for optimization by \code{optim}. Options are \code{"Nelder-Mead"} and \code{"BFGS"}. When using Nelder-Mead, a second optimization with BFGS is run to ensure an optimal fit.}
\item{control}{control list passed to \code{optim} (default \code{list(reltol=1e-10)}).}
\item{keep.data}{whether to attach data to the output object (default \code{TRUE}). When the data object is large, \code{keep.data} can be set to \code{FALSE} to preserve memory. The attached data is used to compute profile likelihood confidence intervals, but can also be supplied separately when \code{keep.data=FALSE}. See \code{\link[=confint.amerasfit]{confint}}.}
\item{na.action}{function or function name controlling missing-value handling for model inputs after formula terms have been expanded. By default, \code{getOption("na.action")} is used, typically \code{na.omit}. Use \code{na.fail} to reject missing model inputs explicitly, or \code{na.pass} to skip row removal; missing model inputs left by \code{na.pass} will generally cause an error during input checks or model fitting. \code{na.exclude} is accepted for fitting, but residuals and diagnostic plots are returned for the fitted rows rather than padded back to the originally supplied row count.}
\item{\dots}{other arguments, passed to functions such as \code{transform}.}
}
\value{
The output is an object of class \code{amerasfit}. General components are \code{call} 
(the function call to \code{ameras}), \code{formula} (the formula object specifying 
the model), \code{num.rows} (the number of rows used for fitting after applying
\code{na.action} and any family-specific filtering), \code{num.rows.original} (the
number of rows supplied before applying \code{na.action}), \code{na.action} (the
omitted-row object returned by the missing-value action, when applicable),
\code{num.realizations} (the number of dose realizations provided),
\code{transform} (the used transformation 
function, if applicable), \code{transform.jacobian} (the used Jacobian function for 
the transformation, if applicable), \code{other.args} (any other arguments passed to 
\dots), \code{model} (a list containing the specified model components parsed from 
the formula), \code{CI.computed} (logical, whether confidence intervals have been 
attached by \code{\link[=confint.amerasfit]{confint}}), and \code{data} (either the 
data frame used for model fitting when \code{keep.data=TRUE} or \code{NULL} otherwise).

For each method supplied to \code{methods}, the output contains a list with components:
  \item{coefficients}{named vector of model coefficients.}
  \item{sd}{named vector of standard deviations.}
  \item{vcov}{covariance matrix for the full parameter vector. Based on the observed information (negative second derivative of log-likelihood) for RC, ERC, and MCML, and on the obtained samples for FMA and BMA.}
  \item{runtime}{string with the total CPU runtime currently available for the method. After fitting this is the fitting time; after \code{\link[=confint.amerasfit]{confint}} this includes confidence interval computation time.}
  \item{timing}{list with CPU and elapsed timings for fitting, confidence interval computation, and their total. CPU time is used for printed runtime summaries.}
  
For RC, ERC, and MCML  the following additional output is included:
  \item{optim}{a list object with results returned by optim. Components are \code{par} (raw parameters before applying a transformation if applicable), \code{hessian} (Hessian matrix for \code{par}), \code{convergence} (convergence code with 0 indicating convergence and 1 indicating that the maximal number of iterations was reached), and \code{counts} (the number of likelihood function evaluations used during optimization).}
  \item{loglik}{log-likelihood value at the optimum.}

For RC and ERC, the output additionally contains:
  \item{ERC}{logical, whether the output is for ERC (\code{ERC=TRUE}) or RC (\code{ERC=FALSE}).}
  
For BMA  the output additionally contains:
  \item{samples}{MCMC posterior samples, as obtained from \code{nimble}. With \code{nchains.BMA > 1}, this is a list with \code{nchains.BMA} components, each a named matrix with the samples from one chain in its rows and columns corresponding to model parameters. With \code{nchains.BMA = 1}, this is a single named matrix.}
  \item{Rhat}{data frame with two columns, \code{Rhat} and \code{n.eff}. The first  column contains the Gelman-Rubin statistics \eqn{\hat R \geq 1} that can be used to assess convergence of MCMC chains. A value of 1 indicates good convergence and values \eqn{>1.05} indicate poor convergence. The effective sample size \code{n.eff} is a measure of how many independent samples the auto-correlated MCMC samples correspond to. A low effective sample size indicates high correlations and/or poor mixing. For single-chain BMA fits, these diagnostics cannot be estimated and are returned as \code{NA}.}
  \item{included.realizations}{indices of realization exposures that were included to obtain the results.}
  \item{prophaz.timepoints}{for \code{family="prophaz"}, time points defining the intervals on which the estimated baseline hazard is constant; these are \code{prophaz.numints.BMA + 1} time points covering the interval \code{(min(entry), max(exit))}, based on quantiles among observed event times. See Details.}

Finally, for FMA the output additionally contains:
  \item{samples}{data frame with samples generated from the normal distributions associated with parameters estimated for each dose realization.}
  \item{weights}{vector of weights corresponding to the models fit to each included realization}
  \item{included.samples}{the total number of samples included.}
  \item{included.realizations}{indices of realization exposures that were included to obtain results. Fits without a valid variance estimate (i.e., non-invertible Hessian or inverse that is not positive definite) or that reach the maximal number of iterations without convergence are filtered out and not used to obtain results. }

The class \code{amerasfit} supports the methods \code{\link[=print.amerasfit]{print}}, \code{\link[=coef.amerasfit]{coef}}, \code{\link[=confint.amerasfit]{confint}}, \code{\link[=summary.amerasfit]{summary}}, \code{\link[=summary_table.amerasfit]{summary_table}}, \code{\link[=vcov.amerasfit]{vcov}}, \code{\link[=residuals.amerasfit]{residuals}}, \code{\link[=plot.amerasfit]{plot}}, \code{\link{traceplot}}, \code{\link[=convergence.amerasfit]{convergence}}, \code{\link[=Rhat.amerasfit]{Rhat}}, and \code{\link[=included_realizations.amerasfit]{included_realizations}}.
}

\details{
Missing values in model inputs are handled by \code{na.action} after formula
terms have been expanded. Thus missing values in raw variables and in derived
covariate columns, such as factor contrasts, interactions, and spline bases,
are handled consistently before model fitting. The omitted-row action is stored
on the fitted object. When \code{keep.data=FALSE} and data are supplied later to
methods such as \code{confint()}, \code{residuals()}, and \code{plot()}, the
same missing-value policy is reapplied. Unlike
\code{\link[stats]{lm}} and \code{\link[stats]{glm}}, \code{na.exclude} does not
pad residuals or diagnostic plots back to the originally supplied row count;
those outputs are computed on the fitted rows.

Models are specified through formulas of the form \code{Y~dose(dose_expression, model="ERR", deg=1, modifier=~M1+M2)+X1+X2}. Here \code{dose_expression} specifies the dose realization columns and is parsed by \code{eval_select} from the \pkg{tidyselect} package. Useful examples are \code{D1:D1000} if the doses are in a sequence of columns with sequential names such as \code{D1}-\code{D1000}, and \code{all_of(dosevars)} where \code{dosevars} is a vector with the names of all dose columns. Further, \code{model} specifies, for non-Gaussian families, whether to use the exponential dose-response model (\code{model="EXP"}), the linear-exponential model (\code{model="LINEXP"}) or the linear ERR model (\code{model="ERR"}). Next, \code{deg} is used to specify whether a quadratic dose term should (\code{deg=2}) or should not (\code{deg=1}) be estimated for the exponential or linear ERR dose-response model. The \code{modifier} term is optional and used to specify effect modification variables through a formula. Use \code{modifier=~M} for reference-plus-contrast coding, and \code{modifier=~0+M} or \code{modifier=~M-1} for subgroup-specific dose coefficients. Formula modifiers may be binary numeric/logical variables or factors; factor modifiers use the existing factor level order, with the first level as the reference. Interactions, transformed modifier terms, continuous numeric modifiers, character modifiers, and no-intercept coding with multiple modifier variables are not currently supported. The older syntax \code{modifier=M1+M2} is deprecated and remains binary-only. When \code{deg}, \code{modifier}, and \code{model} are not supplied, the defaults are \code{deg=1}, no effect modifiers, and \code{model="ERR"}.

Additional right-hand-side covariates are expanded using \code{\link[stats]{model.matrix}}, so standard formula features such as factor contrasts, interactions, \code{I()} terms, and basis functions can be used. For example, formulas can include terms such as \code{X1*X2}, \code{I(age^2)}, \code{splines::ns(age, df=3)}, or \code{splines::bs(age, df=4)}. These terms must preserve one value per row of the input data; terms that drop rows, such as \code{stats::na.omit(X)}, are not supported. When \code{keep.data=FALSE}, later calls that use supplied data apply the same expanded covariate columns as the original fit. If covariates appear poorly scaled or ill-conditioned, ameras emits a diagnostic warning; centering or scaling continuous covariates such as calendar year can improve numerical optimization.

The matched set variable \code{setnr} required for conditional logistic regression is specified on the right-hand side of the formula through a term \code{strata(setnr)}, and an optional offset variable \code{offset} for Poisson regression similarly through a term \code{offset(offset)}. For conditional logistic regression, matched sets of size 1 and matched sets with no cases are excluded, and matched sets with more than one case are not currently supported. For proportional hazards regression, the left-hand side of the formula should have the form \code{Surv(exit, status)} or \code{Surv(entry, exit, status)}. When entry times are supplied, rows with \code{entry == exit} are excluded because they have zero follow-up time under the \code{entry < t <= exit} risk-set convention.

A transformation can be used to reparametrize parameters for optimization, such that the likelihoods are evaluated at \code{transform(parameters)}, where \code{parameters} are unconstrained. This should be specified when fitting linear excess relative risk and linear-exponential models to ensure nonnegative odds/risk/hazard. The included function \code{\link{transform1}} applies an exponential transformation to the desired parameters, see \code{?transform1}. When supplying a function to \code{transform}, this should be a function of the full parameter vector, returning a full (transformed) parameter vector. In particular, the full parameter vector contains parameters in the following order: \eqn{\alpha_0, \bm \alpha, \beta_1, \beta_2, \bm \beta_{m1}, \bm \beta_{m2}, \sigma}, where \eqn{\bm \alpha}, \eqn{\bm\beta_{m1}} and \eqn{\bm \beta_{m2}} can be vectors, with lengths matching \eqn{\bm X} and \eqn{\bm M}, respectively. \eqn{\sigma} is only included for the linear model (Gaussian family), and no intercept is included for the proportional hazards and conditional logistic models. For the multinomial model, the full parameter vector is the concatenation of \eqn{Z-1} parameter vectors in the order as given above, where \eqn{Z} is the number of outcome categories, with the last category chosen as the referent category. See \code{vignette("transformations", package="ameras")} for an example of how to specify a custom transformation function. 

When no transformation is specified and the linear ERR model is used, \code{transform1} is used for ERR parameters \eqn{\beta_1} and \eqn{\beta_2} by default, with lower limits \eqn{-1/max(D)} for \eqn{\beta_1} in the linear dose-response and \eqn{(0,-1/max(D^2))} for \eqn{(\beta_1,\beta_2)} in the linear-quadratic dose-response, respectively. For the linear-exponential model, a lower limit of 0 is used for \eqn{\beta_1}, and no transformation is used for \eqn{\beta_2}. If reference-plus-contrast effect modifiers are specified, no transformation is used for modifier contrast parameters. With subgroup-specific modifier coding, the default transformation is applied to the subgroup dose-effect parameters. When negative RRs are obtained during optimization, an error will be generated and a different transformation or bounds should be used. All output is returned in the original parametrization. The Jacobian of the transformation (\code{transform.jacobian}) is required when using a transformation. For \code{\link{transform1}}, the Jacobian is given by \code{\link{transform1.jacobian}}. No transformations are used in BMA, and FMA is applied on the parameters using the parametrization as given above with variances obtained using the delta method with the provided Jacobian function. BMA does not currently support subgroup-specific modifier coding.

For BMA, a prior distribution for exposure-response parameters can be chosen when using linear or linear-exponential exposure-response model. The options are normal, horseshoe, and double exponential priors, and the same priors truncated at 0 to yield positive values. In particular:

\itemize{
  \item Normal: \eqn{\beta_j \sim N(0,1000)} for all exposure-response parameters \eqn{\beta_j}
  \item Horseshoe (shrinkage prior): \eqn{\tau \sim \text{Cauchy}(0,1)^+; \lambda_j \sim \text{Cauchy}(0,1)^+; \beta_j \sim N(0, \tau^2 \lambda_j^2)}. Here \eqn{\tau} is shared across all parameters
  \item Double exponential (shrinkage prior): \eqn{\lambda_j \sim \text{Cauchy}(0,1)^+; \beta_j \sim \text{DoubleExponential}(0,\lambda_j)}
}
For all other parameters, and when using the exponential exposure-response model or the Gaussian outcome family, the prior is \eqn{N(0, 1000)}. For the parameter \eqn{\sigma} in the Gaussian family, this prior is truncated at 0. 

Because the proportional hazards model is not available in \code{nimble}, \code{ameras} uses a piecewise constant baseline hazard for Bayesian model averaging. The interval \code{min(entry), max(exit))} is divided into \code{prophaz.numints.BMA} subintervals with cutpoints obtained as quantiles of the distribution of event times among cases, and a baseline hazard parameter is estimated for each subinterval.
}

\seealso{
  \code{\link[=confint.amerasfit]{confint}} for computing confidence intervals,
  \code{\link[=summary.amerasfit]{summary}} for a summary of the fitted model
  including confidence intervals if computed, \code{\link[=coef.amerasfit]{coef}} 
  for extracting coefficients, \code{\link[=vcov.amerasfit]{vcov}} for extracting
  variance-covariance matrices, and \code{\link[=included_realizations.amerasfit]{included_realizations}}
  for inspecting which realizations contributed to model averaging.
}

\examples{
\donttest{
  data(data, package="ameras")
  ameras(Y.gaussian~dose(V1:V10, modifier=~M1+M2)+X1+X2, data=data, family="gaussian")
}
}

\references{
  Roberti, S., Kwon D., Wheeler W., Pfeiffer R. (in preparation).
  ameras: An R Package to Analyze Multiple Exposure Realizations in Association Studies
}
