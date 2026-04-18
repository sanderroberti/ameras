
\name{ameras}
\alias{ameras}
\title{Analyze multiple exposure realizations}

\description{
Fit regression models accounting for exposure uncertainty using multiple Monte Carlo exposure realizations. Six outcome model families are supported. The first is the Gaussian family for continuous outcomes, \deqn{Y_i \sim N(\mu_i, \sigma^2),} with \eqn{\mu_i = \alpha_0 + \bm X_i^T \bm \alpha +\beta_1 D_i+\beta_2 D_i^2 + \bm M_i^T \bm \beta_{m1}D_i + \bm M_i^T \bm \beta_{m2}D_i^2}. Here \eqn{\bm X_i} are covariates, \eqn{D_i} is the exposure with measurement error, and \eqn{\bm M_i} are binary effect modifiers. The quadratic exposure terms and effect modification are optional.

For non-Gaussian families, three relative risk models for the main exposure are supported, the usual exponential \eqn{RR_i=\exp(\beta_1 D_i+\beta_2 D_i^2+ \bm M_i^T \bm \beta_{m1}D_i + \bm M_i^T \bm \beta_{m2} D_i^2)} and the linear excess relative risk (ERR) model \eqn{RR_i= 1+\beta_1 D_i+\beta_2 D_i^2 + \bm M_i^T \bm \beta_{m1}D_i + \bm M_i^T \bm \beta_{m2}D_i^2}, where the quadratic and effect modification terms are optional. Finally, the linear-exponential relative risk model \eqn{RR_i= 1+(\beta_1 + \bm{M}_i^T \bm{\beta}_{m1}) D_i \exp\{(\beta_2+  \bm{M}_i^T \bm{\beta}_{m2})D_i\}} is supported.

The second supported family is logistic regression for binary outcomes, with probabilities \deqn{p_i/(1-p_i)=RR_i\exp(\alpha_0+\bm X_i^T \bm \alpha).}

Third is Poisson regression for counts, \deqn{Y_i \sim \text{Poisson}(\mu_i),} where \eqn{\mu_i=RR_i \exp(\alpha_0 +\bm X_i^T \bm  \alpha)\times \text{offset}_i} with optional offset.

Fourth is proportional hazards regression for time-to-event data, with hazard function \deqn{h(t) = h_0(t)RR_i\exp(\bm X_i^T \bm \alpha),} with \eqn{h_0} the baseline hazard.

Fifth is multinomial logistic regression for a categorical outcome with \eqn{Z>2} outcome categories, with the last category as the referent category (i.e., \eqn{\alpha_{0,Z}=\bm \alpha_{Z}=\beta_{1,Z}=\beta_{2,Z}=\bm \beta_{m1,Z} = \bm \beta_{m2,Z}=0}): \deqn{P(Y_i=z)=RR_i\exp(\alpha_{0,z}+\bm X_i^T \bm \alpha_{z})/\{1+\sum_{s=1}^{Z-1} RR_i\exp(\alpha_{0,s}+\bm X_i^T \bm \alpha_{s})\}}

Sixth is conditional logistic regression for matched case control data, for which \deqn{P\left(Y_i = 1, Y_k = 0 \forall k \neq i \bigg| \sum_{i \in \mathcal{R}} Y_i = 1\right) = RR_i\exp(\bm X_i^T \bm \alpha)/\{\sum_{k \in \mathcal{R}}RR_k\exp(\bm X_k^T \bm \alpha)\},} where \eqn{\mathcal{R}} is the matched set corresponding to individual \eqn{i}.

Methods include regression calibration (Carroll et al. 2006 \doi{10.1201/9781420010138}), extended regression calibration (Little et al. 2023 \doi{10.1038/s41598-023-42283-y}), Monte Carlo maximum likelihood (Stayner et al. 2007 \doi{10.1667/RR0677.1}), frequentist model averaging (Kwon et al. 2023 \doi{10.1371/journal.pone.0290498}), and Bayesian model averaging (Kwon et al. 2016 \doi{10.1002/sim.6635}).

}

\usage{
ameras(formula, data, family="gaussian", methods="RC", transform=NULL,
  transform.jacobian=NULL, inpar=NULL, loglim=1e-30, MFMA=100000, 
  prophaz.numints.BMA=10, ERRprior.BMA="doubleexponential", nburnin.BMA=5000, 
  niter.BMA=20000, nchains.BMA=2, thin.BMA=10, included.replicates.BMA=NULL, 
  optim.method="Nelder-Mead", control=NULL, keep.data=TRUE, ... )
}
\arguments{
\item{formula}{an object of class \code{"formula"} containing the model specification. See Details.}
\item{data}{input data frame.}
\item{family}{outcome model family: \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, \code{"prophaz"}, \code{"multinomial"} or \code{"clogit"} (default \code{"gaussian"}).}
%\item{Y}{name or column index of the outcome variable for linear, binomial, Poisson, multinomial and conditional logistic models, or event indicator variable for the proportional hazards model.}
%\item{dosevars}{names or column indices of exposure replicate vectors.}
%\item{M}{names or column indices of binary effect modifying variables (optional).}
%\item{X}{names or column indices of other covariates (optional).}
%\item{offset}{name or column index of offset variable for Poisson regression (optional).}
%\item{entry}{name or column index of left truncation time variable for proportional hazards regression (optional).}
%\item{exit}{name or column index of exit time variable, required when \code{family=prophaz}.}
%\item{setnr}{name or column index of integer-valued matched set variable, required when \code{family="clogit"}.}
\item{methods}{character vector of one or multiple methods to apply. Options: \code{"RC"}, \code{"ERC"}, \code{"MCML"}, \code{"FMA"}, \code{"BMA"} (default \code{"RC"}).}
%\item{deg}{for \code{doseRRmod="ERR"} and \code{doseRRmod="EXP"}, whether to fit a linear (\code{deg=1}) or linear-quadratic (\code{deg=2}) dose-response model (default linear).}
%\item{doseRRmod}{the functional form of the dose-response relationship; options are exponential RR (\code{"EXP"}), linear ERR (\code{"ERR"}), or linear-exponential RR (\code{"LINEXP"}) (default \code{"ERR"}).}
\item{transform}{function for internal parameter transformation (see Details).}
\item{transform.jacobian}{Jacobian of the transformation function (see Details).}
\item{inpar}{vector of initial values for log-likelihood optimization (optional).}
\item{loglim}{parameter used in likelihood computations to avoid taking the log of very small or negative numbers via \code{log(max(x, loglim))} (default \code{1e-30}).}
\item{MFMA}{number of samples for \code{"FMA"} to compute estimates and CIs (default 100,000).}
\item{prophaz.numints.BMA}{for \code{methods="BMA"} with \code{family="prophaz"}, the number of subintervals with constant baseline hazard (default 10). Cut points are determined based on quantiles of the event time distribution among cases.}
\item{ERRprior.BMA}{prior for dose-related parameters when \code{doseRRmod="ERR"} or \code{"LINEXP"} and \code{methods="BMA"}.
Options: \code{"truncated_normal"}, \code{"truncated_horseshoe"}, \code{"truncated_doubleexponential"}, \code{"normal"}, \code{"horseshoe"}, \code{"doubleexponential"}, see Details (default \code{"doubleexponential"}).}
\item{nburnin.BMA}{number of MCMC burn-in iterations for BMA (default 1,000).}
\item{niter.BMA}{number of MCMC iterations per chain for BMA (default 5,000).}
\item{nchains.BMA}{number of MCMC chains for BMA (default 2).}
\item{thin.BMA}{thinning rate for BMA (default 10).}
\item{included.replicates.BMA}{indices of exposure replicates used in BMA (defaults to all replicates).}
\item{optim.method}{method used for optimization by \code{optim}. Options are \code{"Nelder-Mead"} and \code{"BFGS"}. When using Nelder-Mead, a second optimization with BFGS is run to ensure an optimal fit.}
\item{control}{control list passed to \code{optim} (default \code{list(reltol=1e-10)}).}
\item{keep.data}{whether to attach data to the output object (default \code{TRUE}). When the data object is large, \code{keep.data} can be set to \code{FALSE} to preserve memory. The attached data is used to compute profile likelihood confidence intervals, but can also be supplied seperately when \code{keep.data=FALSE}. See \code{\link[=confint.amerasfit]{confint}}.}
\item{\dots}{other arguments, passed to functions such as \code{transform}.}
}
\value{
The output is an object of class \code{amerasfit}. General components are \code{call} (the function call to \code{ameras}), \code{formula} (the formula object specifying the model), \code{num.rows} (the number of rows in \code{data}), \code{num.replicates} (the number of dose replicates provided), \code{transform} (the used transformation function, if applicable), \code{transform.jacobian} (the used Jacobian function for the transformation, if applicable), \code{other.args} (any other arguments passed to \dots), \code{model} (a list containing the specified model components parsed from the formula), and \code{CI.computed} (logical, whether confidence intervals have been attached by \code{\link[=confint.amerasfit]{confint}}).

For each method supplied to \code{methods}, the output contains a list with components:
  \item{coefficients}{named vector of model coefficients.}
  \item{sd}{named vector of standard deviations.}
  \item{runtime}{string with the runtime in seconds.}
  
For RC, ERC, and MCML  the following additional output is included:
  \item{vcov}{covariance matrix for the full parameter vector.}
  \item{optim}{a list object with results returned by optim. Components are \code{par} (raw parameters before applying a transformation if applicable), \code{hessian} (Hessian matrix for \code{par}), \code{convergence} (convergence code with 0 indicating convergence and 1 indicating that the maximal number of iterations was reached), and \code{counts} (the number of likelihood function evaluations used during optimization).}
  \item{loglik}{log-likelihood value at the optimum.}

For RC and ERC, the output additionally contains:
  \item{ERC}{logical, whether the output is for ERC (\code{ERC=TRUE}) or RC (\code{ERC=FALSE}).}
  
For BMA  the output additionally contains:
  \item{samples}{MCMC posterior samples, as obtained from \code{nimble}. This is a list object with \code{nchains.BMA} components, each a named matrix with the samples from one chain in its rows, with columns corresponding to model parameters.} 
  \item{Rhat}{data frame with two columns, \code{Rhat} and \code{n.eff}. The first  column contains the Gelman-Rubin statistics \eqn{\hat R \geq 1} that can be used to assess convergence of MCMC chains. A value of 1 indicates good convergence and values \eqn{>1.05} indicate poor convergence. The effective sample size \code{n.eff} is a measure of how many independent samples the auto-correlated MCMC samples correspond to. A low effective sample size indicates high correlations and/or poor mixing.}
  \item{included.replicates}{indices of replicate exposures that were included to obtain the results.}
  \item{prophaz.timepoints}{for \code{family="prophaz"}, time points defining the intervals on which the estimated baseline hazards is constant; these are \code{prophaz.numints.BMA + 1} time points covering the interval \code{(min(entry), max(exit))}, based on quantiles among observed event times. See Details.}

Finally, for FMA the output additionally contains:
  \item{samples}{the samples generated from the normal distributions associated with each dose replicate.}
  \item{included.samples}{the total number of samples included.}
  \item{included.replicates}{indices of replicate exposures that were included to obtain results. Fits without a valid variance estimate (i.e., non-invertible Hessian or inverse that is not positive definite) or that reach the maximal number of iterations without convergence are filtered out and not used to obtain results. }

The class \code{amerasfit} supports the methods \code{print}, \code{coef}, \code{\link[=confint.amerasfit]{confint}}, \code{summary}, and \code{\link{traceplot}}.
}

\details{
A transformation can be used to reparametrize parameters internally (i.e., such that the likelihoods are evaluated at \code{transform(parameters)}, where \code{parameters} are unconstrained), and should be specified when fitting linear excess relative risk and linear-exponential models to ensure nonnegative odds/risk/hazard. The included function \code{\link{transform1}} applies an exponential transformation to the desired parameters, see \code{?transform1}. When supplying a function to \code{transform}, this should be a function of the full parameter vector, returning a full (transformed) parameter vector. In particular, the full parameter vector contains parameters in the following order: \eqn{\alpha_0, \bm \alpha, \beta_1, \beta_2, \bm \beta_{m1}, \bm \beta_{m2}, \sigma}, where \eqn{\bm \alpha}, \eqn{\bm\beta_{m1}} and \eqn{\bm \beta_{m2}} can be vectors, with lengths matching \eqn{\bm X} and \eqn{\bm M}, respectively. \eqn{\sigma} is only included for the linear model (Gaussian family), and no intercept is included for the proportional hazards and conditional logistic models. For the multinomial model, the full parameter vector is the concatenation of \eqn{Z-1} parameter vectors in the order as given above, where \eqn{Z} is the number of outcome categories, with the last category chosen as the referent category. See \code{vignette("transformations", package="ameras")} for an example of how to specify a custom transformation function. 

When no transformation is specified and the linear ERR model is used, \code{transform1} is used for ERR parameters \eqn{\beta_1} and \eqn{\beta_2} by default, with lower limits \eqn{-1/max(D)} for \eqn{\beta_1} in the linear dose-response and \eqn{(0,-1/max(D^2))} for \eqn{(\beta_1,\beta_2)} in the linear-quadratic dose-response, respectively. For the linear-exponential model, a lower limit of 0 is used for \eqn{\beta_1}, and no transformation is used for \eqn{\beta_2}. If effect modifiers \code{M} are specified, no transformation is used for those parameters. When negative RRs are obtained during optimization, an error will be generated and a different transformation or bounds should be used. All output is returned in the original parametrization. The Jacobian of the transformation (\code{transform.jacobian}) is required when using a transformation. For \code{\link{transform1}}, the Jacobian is given by \code{\link{transform1.jacobian}}. No transformations are used in BMA, and FMA is applied on the parameters using the parametrization as given in above with variances obtained using the delta method with the provided Jacobian function.  

For BMA, a prior distribution for exposure-response parameters can be chosen when using linear or linear-exponential exposure-response model. The options are normal, horshoe, and double exponential priors, and the same priors truncated at 0 to yield positive values. In particular:

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
  for extracting coefficients.
}

\examples{
\donttest{
  data(data, package="ameras")
  ameras(Y.gaussian~dose(V1:V10, modifier=M1+M2)+X1+X2, data=data, family="gaussian") 
}
}

\references{
  Roberti, S., Kwon D., Wheeler W., Pfeiffer R. (in preparation).
  ameras: An R Package to Analyze Multiple Exposure Realizations in Association Studies
}