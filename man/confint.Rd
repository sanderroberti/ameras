\name{confint.amerasfit}
\alias{confint.amerasfit}
\title{Confidence Intervals for an amerasfit Object}

\description{
  Computes confidence or credible intervals for the parameters of a fitted
  \code{amerasfit} object and returns the interval bounds. The fitted object
  is not amended. Use \code{\link{add_confint}} to attach intervals and any
  available profile diagnostics to an \code{amerasfit} object for inclusion
  in \code{summary()} output and recording of confidence-interval timing.
}

\usage{
\method{confint}{amerasfit}(object, parm="dose", level=0.95,
        type=c("proflik","percentile"), maxit.profCI=20,
        tol.profCI=1e-2, data=NULL, \dots)
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
    Maximum number of root-finding iterations for profile likelihood
    intervals. Defaults to \code{20}.
  }
  \item{tol.profCI}{
    Root-finding tolerance for profile likelihood intervals. Defaults to
    \code{1e-2}.
  }
  \item{data}{
    Original fitting data. Required for profile likelihood intervals when the
    model was fitted with \code{keep.data = FALSE}.
  }
  \item{\dots}{
    Additional arguments, currently unused. To attach intervals to a fitted
    object or control printed output, use
    \code{\link{add_confint}}.
  }
}

\value{
  For an object containing one fitted method, a numeric matrix with parameters
  in rows and the lower and upper confidence limits in two percentage-labelled
  columns. For an object containing multiple fitted methods, a named list of
  such matrices, one per method.
}

\details{
  For RC, ERC, and MCML, Wald intervals on the original parameter scale use
  the delta-method variance-covariance matrix. Wald intervals on the
  transformed scale are computed in the optimization parametrization and then
  back-transformed. More specifically, when
  \eqn{\boldsymbol\theta = h(\boldsymbol\eta)}, transformed-scale intervals
  apply \eqn{h} to the endpoints on the \eqn{\boldsymbol\eta} scale, whereas
  original-scale intervals use the covariance matrix obtained from
  \eqn{J H^{-1} J^T}, where \eqn{J} is the transformation Jacobian and
  \eqn{H} is the optimization Hessian. Profile likelihood intervals repeatedly
  re-evaluate the likelihood and can therefore be substantially more
  expensive.

  For FMA and BMA, percentile intervals use equal-tailed sample quantiles and
  HPD intervals use \code{\link[coda]{HPDinterval}}.

  Profile likelihood diagnostics such as boundary p-values and iteration
  counts are not included in the returned matrices. Use
  \code{\link{add_confint}} to store these diagnostics in the method-specific
  \code{CI} components of a fitted object.

  The requested intervals are always computed, even if \code{object} already
  contains attached intervals. The input object is left unchanged. If a
  profile likelihood bound cannot be reliably bracketed or solved, that bound
  is returned as \code{NA} with a diagnostic warning.

  When a model was fitted with \code{keep.data = FALSE}, the original data must
  be supplied through \code{data} for profile likelihood intervals. The
  missing-value policy used during fitting is reapplied before likelihood
  evaluation. Wald and sample-based intervals do not require the original
  data. For models with formula effect modifiers, subgroup-specific profile
  intervals can be obtained by fitting subgroup-specific dose coefficients
  with \code{modifier = ~ 0 + M} or \code{modifier = ~ M - 1}.
}

\seealso{
  \code{\link{add_confint}} for attaching intervals and diagnostics to a fit,
  \code{\link{ameras}} for model fitting,
  \code{\link[=summary.amerasfit]{summary}} for fit summaries, and
  \code{\link[stats]{confint}} for the generic function.
}

\examples{
fit <- ameras(
  Y.binomial ~ dose(V1:V3, model = "ERR"),
  data = ameras_simulated,
  family = "binomial",
  methods = "RC"
)

## A two-column matrix of Wald interval bounds
confint(fit, type = "wald.orig")

## Profile likelihood intervals are slower
\donttest{% slow
confint(fit, type = "proflik", parm = "dose")
}
}
