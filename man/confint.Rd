\name{confint.amerasfit}
\alias{confint.amerasfit}
\alias{add_confint}
\alias{add_confint.amerasfit}
\title{Compute or Attach Confidence Intervals for an amerasfit Object}

\description{
  Compute confidence or credible intervals for the parameters of a fitted
  \code{amerasfit} object. \code{confint()} returns standalone interval
  bounds. \code{add_confint()} attaches bounds and any available diagnostics
  to a copy of the fitted object for inclusion in \code{summary()} output,
  and records confidence-interval timing.
}

\usage{
\method{confint}{amerasfit}(object, parm="dose", level=0.95,
        type=c("proflik","percentile"), maxit.profCI=100,
        tol.profCI=1e-6, data=NULL, \dots)

add_confint(object, \dots)

\method{add_confint}{amerasfit}(object, parm="dose", level=0.95,
        type=c("proflik","percentile"), maxit.profCI=100,
        tol.profCI=1e-6, data=NULL, \dots)
}

\arguments{
  \item{object}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}.
  }
  \item{parm}{
    Either \code{"dose"} to compute profile likelihood intervals for
    dose-related parameters only, \code{"all"} for all parameters, or a
    character vector of parameter names. Only used when
    \code{type = "proflik"}. Defaults to \code{"dose"} because profile
    likelihood intervals can be computationally intensive.
  }
  \item{level}{
    Confidence level, as a single number between zero and one. Defaults to
    \code{0.95}.
  }
  \item{type}{
    Interval type. For RC, ERC, and MCML, use \code{"wald.orig"},
    \code{"wald.transformed"}, or \code{"proflik"}. For FMA and BMA, use
    \code{"percentile"} or \code{"hpd"}. When both groups of methods are
    present, supply one type for each group.
  }
  \item{maxit.profCI}{
    Maximum number of root-finding iterations per profile likelihood bound.
    Defaults to \code{100}. This does not control the
    separate bracketing search or the nuisance-parameter optimizer.
  }
  \item{tol.profCI}{
    Root-finding tolerance on the internal optimizer parameter scale. Defaults to \code{1e-6}. Larger values can reduce
    computation time, particularly for ERC and MCML, but may produce
    inaccurate bounds. A bound is returned as \code{NA} with a warning if
    root finding warns or its p-value differs from \code{1 - level} by more
    than \code{0.005}, even when a larger tolerance is explicitly requested.
  }
  \item{data}{
    Original fitting data. Required for profile likelihood intervals when the
    model was fitted with \code{keep.data = FALSE}.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  \describe{
    \item{\code{confint()}}{
      For an object containing one fitted method, a numeric matrix with
      parameters in rows and the lower and upper confidence limits in two
      percentage-labelled columns. For an object containing multiple fitted
      methods, a named list of such matrices, one per method. The result is
      returned visibly and contains only interval bounds.
    }
    \item{\code{add_confint()}}{
      A copy of the \code{amerasfit} object, returned invisibly, with a
      \code{CI} element added to each fitted method. Method-specific timing is updated to include
      interval computation. Each \code{CI} element contains \code{lower} and
      \code{upper}. Profile likelihood results additionally contain boundary
      p-values and iteration counts in \code{pval.lower}, \code{pval.upper},
      \code{iter.lower}, and \code{iter.upper}.
    }
  }
}

\details{
  Both functions use the same interval calculations and always compute the
  requested intervals, even if \code{object} already contains attached
  intervals. The input object is left unchanged.

  Use \code{confint()} when only standalone bounds are needed, for example
  for comparison or export. Use \code{add_confint()} when the fit will be
  retained for reporting or further inspection: \code{summary()} and
  \code{coef(summary(object))} include the attached bounds, profile-likelihood
  diagnostics are stored when available, and computation time is recorded.
  Assign the invisibly returned object to retain this information. Any
  previously attached intervals are replaced in that returned copy.
  Use \code{summary()} to display
  the stored intervals and \code{print(summary(fit), digits = 3)}, for example,
  to control the printed precision.
  Attaching intervals is particularly useful for expensive profile likelihood
  calculations because subsequent summaries use the stored results without
  repeating the calculation.
}

\section{Interval types}{
  For RC, ERC, and MCML, Wald intervals on the original parameter scale use
  the delta-method variance-covariance matrix. Wald intervals on the
  transformed scale are computed in the optimization parametrization and then
  back-transformed. More specifically, when
  \eqn{\boldsymbol\theta = h(\boldsymbol\eta)}, transformed-scale intervals
  apply \eqn{h} to the endpoints on the \eqn{\boldsymbol\eta} scale, whereas
  original-scale intervals use the covariance matrix obtained from
  \eqn{J H^{-1} J^T}, where \eqn{J} is the transformation Jacobian and
  \eqn{H} is the optimization Hessian. Transformed-scale intervals can
  therefore be asymmetric on the reported parameter scale.
  Profile likelihood intervals repeatedly re-evaluate the likelihood and can
  therefore be substantially more expensive.

  Transformations must be componentwise and strictly monotone, as required by
  \code{\link{ameras}}. Decreasing transformations reverse endpoint order;
  bounds and profile diagnostics are returned in reported lower/upper order.
  Transformed Wald and profile likelihood intervals check the transformation
  at the fitted parameters and determine its directions locally. Original-scale
  Wald intervals and sample-based intervals do not require this check.

  For FMA and BMA, percentile intervals use equal-tailed sample quantiles and
  HPD intervals use \code{\link[coda]{HPDinterval}}.
}

\section{Profile likelihood bracketing}{
  Before determining bounds, a separate adaptive search first finds a bracket for each lower and upper
  bound. A bracket is an interval within which the
  likelihood-ratio p-value crosses the target \code{1 - level}. The search
  proceeds outward from the fitted parameter value on the internal optimizer
  scale. Initial search distances are four approximate standard errors when
  the stored Hessian is usable; otherwise, fallback starting points are used.
  If necessary, the search distance is doubled up to 12 times. Each profile
  evaluation may require re-optimizing the nuisance parameters, so bracketing
  can itself be computationally expensive, particularly for ERC and MCML.

  Once a bracket is found, \code{\link[stats]{uniroot}} locates the bound
  within it using \code{tol.profCI} and \code{maxit.profCI}. These arguments
  control root finding only; they do not affect the bracketing search. If a
  bound cannot be bracketed or reliably solved, it is returned as \code{NA}
  with a warning reporting the last attempted parameter value and its
  p-value. A bracketing failure can reflect a limited search range, parameter
  scaling, numerical problems in the profile evaluations, or the absence of
  a finite bound.

  Nuisance optimizations use the effective fitting controls. Missing or
  \code{NULL} controls on older objects default to
  \code{list(reltol = 1e-10)}. When only
  one nuisance parameter remains, \code{optimize()} uses its own tolerance;
  \code{reltol} does not control that optimization.

  When a model was fitted with \code{keep.data = FALSE}, the original data must
  be supplied through \code{data} for profile likelihood intervals. The
  missing-value policy used during fitting is reapplied before likelihood
  evaluation. Wald and sample-based intervals do not require the original
  data. For models with formula effect modifiers, subgroup-specific profile
  intervals can be obtained by fitting subgroup-specific dose coefficients
  with \code{modifier = ~ 0 + M} or \code{modifier = ~ M - 1}.
}

\seealso{
  \code{\link{ameras}} for model fitting,
  \code{\link[=summary.amerasfit]{summary}} for fit summaries, and
  \code{\link[stats]{confint}} for the generic function.
  See also \code{vignette("confidenceintervals", package = "ameras")}.
}

\examples{
fit <- ameras(
  Y.binomial ~ dose(V1:V3, model = "ERR"),
  data = ameras_simulated,
  family = "binomial",
  methods = "RC"
)

## Compute standalone Wald bounds
confint(fit, type = "wald.orig")

## Attach Wald intervals for use in summaries
fit_with_ci <- add_confint(fit, type = "wald.orig")
summary(fit_with_ci)

## Attach expensive profile intervals for subsequent reporting
\donttest{% slow
fit_with_profile_ci <- add_confint(
  fit,
  type = "proflik",
  parm = "dose"
)
summary(fit_with_profile_ci)
}
}
