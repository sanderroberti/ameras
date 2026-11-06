\name{add_confint}
\alias{add_confint}
\alias{add_confint.amerasfit}
\title{Attach Confidence Intervals to an amerasfit Object}

\description{
  Computes confidence or credible intervals, attaches them and any available
  diagnostics to a fitted \code{amerasfit} object, and optionally prints the
  bounds. Use \code{\link[=confint.amerasfit]{confint}} when only the interval
  bounds are required.
}

\usage{
add_confint(object, \dots)

\method{add_confint}{amerasfit}(object, parm="dose", level=0.95,
        type=c("proflik","percentile"), maxit.profCI=20,
        tol.profCI=1e-2, data=NULL, print=TRUE,
        digits=max(3, getOption("digits") - 3), \dots)
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
    \code{type = "proflik"}.
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
  \item{print}{
    Whether to print the attached interval bounds. Defaults to \code{TRUE}.
  }
  \item{digits}{
    Number of significant digits used when printing.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  A copy of the \code{amerasfit} object, invisibly, with a \code{CI} element
  added to each fitted method and \code{CI.computed} set to \code{TRUE}.
  Method-specific timing is updated to include interval computation. The input
  object is not modified.

  Each \code{CI} element contains \code{lower} and \code{upper}. Profile
  likelihood results additionally contain boundary p-values and iteration
  counts in \code{pval.lower}, \code{pval.upper}, \code{iter.lower}, and
  \code{iter.upper}.
}

\details{
  The requested intervals are always computed. Any intervals already attached
  to \code{object} are replaced in the returned object; the input object is
  left unchanged. This function is useful when the fitted object will be kept
  for reporting or further inspection: \code{summary()} and
  \code{coef(summary(object))} include the attached bounds, profile-likelihood
  diagnostics are stored when available, and confidence-interval computation
  time is recorded in the method-specific timing information. Use
  \code{\link[=confint.amerasfit]{confint}} when only standalone interval
  bounds are needed. Attaching intervals is particularly useful for
  computationally expensive profile likelihood intervals because subsequent
  summaries use the attached bounds and diagnostics without repeating the
  profile calculation. Computation and failure handling are otherwise
  identical.

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

  For FMA and BMA, percentile intervals use equal-tailed sample quantiles and
  HPD intervals use \code{\link[coda]{HPDinterval}}.
}

\seealso{
  \code{\link[=confint.amerasfit]{confint}} for returning interval bounds,
  \code{\link{ameras}} for model fitting, and
  \code{\link[=summary.amerasfit]{summary}} for summaries containing attached
  intervals.
}

\examples{
fit <- ameras(
  Y.binomial ~ dose(V1:V3, model = "ERR"),
  data = ameras_simulated,
  family = "binomial",
  methods = "RC"
)

fit_with_ci <- add_confint(fit, type = "wald.orig", print = FALSE)
summary(fit_with_ci)

## Profile likelihood intervals are slower
\donttest{% slow
fit_with_profile_ci <- add_confint(
  fit,
  type = "proflik",
  parm = "dose",
  print = FALSE
)
summary(fit_with_profile_ci)
}
}
