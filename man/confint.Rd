\name{confint.amerasfit}
\alias{confint.amerasfit}
\title{Confidence intervals for an amerasfit object}

\description{
  Computes confidence intervals for the parameters of a fitted
  \code{amerasfit} object. This is a separate step from model fitting,
  i.e., \code{\link{ameras}} fits the model and \code{confint} computes 
  intervals and attaches them to the fitted object.
}

\usage{
\method{confint}{amerasfit}(object, parm="dose", level=0.95, 
        type=c("proflik","percentile"), maxit.profCI=20, tol.profCI=1e-2, 
        data=NULL, \dots)
}

\arguments{
  \item{object}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}.
  }
  \item{parm}{
    Either \code{"dose"} to compute intervals for dose-related parameters
    only, \code{"all"} for all parameters, or a character vector of specific
    parameter names. Only used when \code{type = "proflik"} since Wald
    intervals are cheap to compute for all parameters simultaneously.
    Defaults to \code{"dose"} since profile likelihood computation can be
    extensive.
  }
  \item{level}{
    The confidence level. Defaults to \code{0.95}.
  }
  \item{type}{
    The type(s) of confidence intervals to determine. For
    RC, ERC, and MCML, this can be one of:
    \describe{
      \item{\code{"wald.orig"}}{
        Wald intervals on the original parameter scale using the delta
        method variance-covariance matrix.
      }
      \item{\code{"wald.transformed"}}{
        Wald intervals computed on the transformed (reparametrized) scale
        and then back-transformed. Only available when a transformation
        was used during fitting.
      }
      \item{\code{"proflik"}}{
        Profile likelihood intervals based on the chi-squared
        approximation. More accurate than Wald intervals but
        computationally intensive, especially for large datasets or
        complex models.
      }
    }
    For FMA and BMA, confidence intervals are based on the generated samples and
    possible confidence interval types are:
    \describe{
      \item{\code{"percentile"}}{
        Equal-tailed percentile intervals.
      }
      \item{\code{"hpd"}}{
        Highest posterior density intervals via
        \code{\link[coda]{HPDinterval}} from the \pkg{coda} package.
      }
    }
    If \code{object} contains results for at least one of RC, ERC, and MCML and at least one of FMA and BMA, \code{type} must be length 2 and specify one method for RC, ERC, and MCML, and one for FMA and BMA.
  }
  \item{data}{
    The original data frame used for fitting. Only required when
    \code{type = "proflik"} and the model was fitted with
    \code{keep.data = FALSE}
  }
  \item{maxit.profCI}{
    Maximum number of iterations for the root-finding algorithm used to
    locate profile likelihood interval bounds. Only used when
    \code{type = "proflik"}. Defaults to \code{20}.
  }
  \item{tol.profCI}{
    Tolerance for the root-finding algorithm. Only used when
    \code{type = "proflik"}. Defaults to \code{1e-2}. Reduce for more
    precise bounds at the cost of additional computation.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  The original \code{amerasfit} object with a \code{CI} element added to
  each fitted method result. For RC, ERC, and MCML the \code{CI} element
  is a data frame with columns:
  \describe{
    \item{\code{lower}}{Lower confidence bound.}
    \item{\code{upper}}{Upper confidence bound.}
  }
  When \code{type = "proflik"}, four additional columns are included:
  \describe{
    \item{\code{pval.lower}}{
      P-value at the lower bound, should be close to \eqn{1 -}
      \code{level}.
    }
    \item{\code{pval.upper}}{
      P-value at the upper bound, should be close to \eqn{1 -}
      \code{level}.
    }
    \item{\code{iter.lower}}{
      Number of iterations used by the root-finding algorithm for the
      lower bound.
    }
    \item{\code{iter.upper}}{
      Number of iterations used by the root-finding algorithm for the
      upper bound.
    }
  }
  For FMA and BMA the \code{CI} element is a data frame with columns
  \code{lower} and \code{upper}.
}

\details{

  For (extended) regression calibration and Monte Carlo maximum likelihood, 
  Wald and profile likelihood intervals can be obtained. When a parameter 
  transformation \eqn{\bm\theta = h(\bm\eta)} is used, \code{CI="wald.transformed"} 
  yields the CI \eqn{h(\bm\eta \pm 1.96 \bm V)} with \eqn{\bm V} the vector 
  of standard deviations estimated using the inverse Hessian matrix, 
  and \code{CI="wald.orig"} uses the delta method to obtain the CI 
  \eqn{h(\bm\eta)\pm 1.96 \bm V_*} where \eqn{\bm V_*} is the vector of 
  standard deviations estimated using \eqn{J H^{-1} J^T} with \eqn{J} the 
  Jacobian of the transformation and \eqn{H} is the Hessian. When no 
  transformation is used, \code{CI="wald.orig"} should be used. 
  The third option is \code{proflik}, which uses the profile likelihood to 
  compute confidence bounds. For FMA and BMA, the options for 
  confidence/credible intervals are \code{CI="percentile"} which uses 
  percentiles, and \code{CI="hpd"} which computes highest posterior density 
  intervals using \code{HPDinterval} from the \code{coda} package, both using 
  the FMA samples or Bayesian posterior samples.
  
  Profile likelihood intervals (\code{type = "proflik"}) require
  re-evaluating the likelihood repeatedly and can be time-consuming.
  The \code{parm} argument can be used to restrict computation to dose
  parameters only (the default) when intervals for the other parameters are
  not of interest.

  When the model was fitted with \code{keep.data = FALSE} and
  \code{type = "proflik"} is used for \code{confint}, the original data must be
  supplied via the \code{data} argument. Wald intervals do not require the
  data and can always be computed from the stored Hessian and parameter
  estimates alone.
}

\seealso{
  \code{\link{ameras}} for model fitting,
  \code{\link[=summary.amerasfit]{summary}} for a summary of the fitted model
  including confidence intervals if computed,
  \code{\link[stats]{confint}} for the generic function.
}

\examples{
data("data", package = "ameras")

## Fit the model
fit <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
               methods = "RC")

## Wald intervals (fast)
fit <- confint(fit, type = "wald.orig")
summary(fit)

## Profile likelihood intervals for dose parameters only (slower)
\donttest{
fit <- confint(fit, type = "proflik", parm = "dose")
summary(fit)
}

## With keep.data = FALSE, supply data explicitly for proflik
\donttest{
fit2 <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
               methods = "RC", keep.data = FALSE)
fit2 <- confint(fit2, type = "proflik", data = data)
}

## FMA and BMA with percentile intervals
\donttest{
fit3 <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial", 
               methods = c("FMA", "BMA"))
fit3 <- confint(fit3, type = "percentile")
summary(fit3)
}
}
