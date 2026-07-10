\name{dose_lrt}
\alias{dose_lrt}
\alias{dose_lrt.amerasfit}
\title{Likelihood-ratio tests for dose-related parameters}
\description{
  Compute global and individual likelihood-ratio tests for dose-related
  parameters in supported \code{amerasfit} methods.
}
\usage{
dose_lrt(object, ...)

\method{dose_lrt}{amerasfit}(
  object,
  methods = c("RC", "ERC", "MCML"),
  type = c("global", "individual"),
  data = NULL,
  ...
)
}
\arguments{
  \item{object}{An object of class \code{amerasfit}.}
  \item{methods}{
    Character vector of fitted methods to test. Supported methods are
    \code{"RC"}, \code{"ERC"}, and \code{"MCML"}.
  }
  \item{type}{
    Character vector selecting which tests to return. \code{"global"} tests all
    dose-related parameters jointly. \code{"individual"} tests each
    dose-related parameter separately. By default, both are returned.
  }
  \item{data}{
    Optional data frame required when \code{object} was fitted with
    \code{keep.data=FALSE}.
  }
  \item{\dots}{Currently unused.}
}
\value{
  A data frame with one row per likelihood-ratio test and columns:
  \describe{
    \item{\code{method}}{The ameras estimation method.}
    \item{\code{type}}{Either \code{"global"} or \code{"individual"}.}
    \item{\code{term}}{The tested dose term. For global tests this is
    \code{"all dose terms"}.}
    \item{\code{df}}{Degrees of freedom for the likelihood-ratio test.}
    \item{\code{logLik}}{The fitted model log likelihood.}
    \item{\code{logLik.null}}{The constrained null-model log likelihood.}
    \item{\code{statistic}}{The likelihood-ratio statistic.}
    \item{\code{p.value}}{The chi-square reference p-value.}
    \item{\code{null.optim.convergence}}{Optimizer convergence code for the
    constrained null fit. A warning is issued if this code is nonzero.}
  }
}
\details{
  \code{dose_lrt()} computes likelihood-ratio tests for dose-related parameters
  in fitted RC, ERC, and MCML models. The null hypothesis sets the tested
  dose-related parameter(s) equal to zero. FMA and BMA are not
  supported by this function.

  For transformed fits, nuisance parameters are evaluated through the fitted
  transformation. The tested dose parameters are then fixed to zero before the
  likelihood is evaluated. This is most straightforward for component-wise
  transformations, including the default transformations used by \pkg{ameras}.
  Results from coupled custom transformations should be interpreted with
  caution. For example, if a custom transformation defines one reported
  parameter as a function of multiple optimizer-scale parameters, or makes a
  non-dose reported parameter depend directly on a dose parameter, fixing only
  the reported dose parameter to zero may not correspond to the intended
  constrained model.

  For speed, constrained null fits report only the optimizer convergence code.
  They do not run the additional numerical-gradient diagnostics stored for the
  original RC, ERC, and MCML fits.
}
\examples{
set.seed(1)
dat <- data.frame(
  Y = 1 + 0.8 * seq(-1, 1, length.out = 30) + rnorm(30, sd = 0.2),
  D = seq(-1, 1, length.out = 30)
)

fit <- ameras(Y ~ dose(D), data = dat, family = "gaussian",
              methods = "RC", print = FALSE)

dose_lrt(fit)
dose_lrt(fit, type = "global")
}
