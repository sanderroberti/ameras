\name{plot.amerasfit}
\alias{plot.amerasfit}
\title{Diagnostic Plots for an amerasfit Object}

\description{
  Produces diagnostic plots for a fitted \code{amerasfit} object,
  including residuals versus fitted values and normal Q-Q plots.
  Plots are produced for each estimation method present in the fitted
  object.
}

\usage{
\method{plot}{amerasfit}(x,
     methods = c("RC", "ERC", "MCML", "FMA", "BMA"), 
     which = NULL, type = NULL, dose.col = NULL, 
     add.smooth = getOption("add.smooth", TRUE), 
     qqline = TRUE, id.n = 3, ask = NULL, data = NULL, \dots)
}

\arguments{
  \item{x}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}.
  }
  \item{methods}{
    Character vector specifying which estimation methods to produce
    plots for. One or more of \code{"RC"}, \code{"ERC"},
    \code{"MCML"}, \code{"FMA"}, and \code{"BMA"}. Defaults to all. 
    Methods not present in \code{x} are silently skipped.
  }
  \item{which}{
    Character vector specifying which plots to produce. For non-\code{"prophaz"}
    families, one or both of \code{"residuals-vs-fitted"} and \code{"qq"}.
    For \code{family="prophaz"}, \code{"schoenfeld"} produces a scaled
    Schoenfeld residual plot against event time. Defaults to both for
    non-\code{"prophaz"} families and to \code{"schoenfeld"} for
    \code{"prophaz"} families.
  }
  \item{type}{
    The type of residuals to use. One of \code{"pearson"},
    \code{"deviance"}, or \code{"response"} for non-\code{"prophaz"}
    models. For \code{family="prophaz"}, this must be
    \code{"schoenfeld"}. Defaults to \code{"pearson"} for non-\code{"prophaz"}
    families and to \code{"schoenfeld"} for \code{"prophaz"}. 
    See \code{\link{residuals.amerasfit}} for details.
  }
  \item{dose.col}{
    The dose realization to use for computation of residuals. If
    \code{NULL} (default), a best fitting realization is determined
    for each method (see Details).
  }
  \item{add.smooth}{
    Logical. If \code{TRUE}, a loess smooth line is added to the
    residuals versus fitted plot via \code{\link{panel.smooth}}.
    Defaults to \code{getOption("add.smooth", TRUE)}. Is not used 
    for the Schoenfeld residual plot for which a spline is always 
    shown, mirroring \link[survival]{plot.cox.zph}.
  }
  \item{qqline}{
    Logical. If \code{TRUE}, a reference line is added to normal
    Q-Q plots via \code{\link[stats]{qqline}}. Defaults to \code{TRUE}.
  }
  \item{id.n}{
    Integer. The number of extreme residuals to label in residuals vs fitted 
    and normal Q-Q plots. Labels show the row index of the observation. Defaults to
    \code{3}. Set to \code{0} to suppress labeling.
  }
  \item{ask}{
    Logical. If \code{TRUE}, the user is prompted before each new
    plot. Defaults to \code{NULL}, in which case prompting occurs
    automatically when the number of plots exceeds the number of
    available panels and the session is interactive.
  }
  \item{data}{
    The original data frame used for fitting. Only required when the
    model was fitted with \code{keep.data=FALSE}.
  }
  \item{\dots}{
    Additional graphical arguments passed to \code{\link{plot.default}} and
    \code{\link{qqnorm}}.
  }
}

\value{
  The \code{amerasfit} object \code{x} is returned invisibly.
}

\details{
  If not otherwise specified, the dose realization used to compute fitted
  values is selected for each estimation method as follows. For RC and 
  ERC the mean dose across realizations is used. For MCML and FMA the 
  realization yielding the largest likelihood at the final parameter 
  estimates is used (which for FMA corresponds to the highest model 
  averaging weight). For BMA the realization most frequently selected by 
  the MCMC sampler is used. The selected dose column is shown in the plot title.

  For the \code{"residuals-vs-fitted"} plot, the smooth line is added
  via \code{\link{panel.smooth}} using a loess smoother when \code{add.smooth=TRUE}. 
  For the \code{"qq"} plot, the residuals are plotted against the
  theoretical quantiles of the normal distribution when \code{qqline=TRUE}. For both plots, 
  the \code{id.n} most extreme residuals are labeled. For \code{family="multinomial"}, 
  one panel is produced per non-reference outcome category.

  For proportional hazards models, scaled Schoenfeld residuals are drawn
  against the observed event times to assess the proportional hazards assumption. 
  Under proportional hazards, the residuals should fluctuate randomly around
  zero with no systematic trend over time. Systematic patterns or smooth
  trends may indicate time-varying covariate effects and violation of the
  proportional hazards assumption. Note that the implementation here corresponds 
  to using \code{transform = "identity"} in \code{\link[survival]{cox.zph}}, 
  in contrast to the default setting for that function which applies a 
  Kaplan--Meier transformation of time.

  The data must either be stored on the object (\code{keep.data=TRUE} in 
  \code{ameras}, the default) or supplied via the \code{data} argument.

}

\seealso{
  \code{\link{residuals.amerasfit}} for computing residuals,
  \code{\link{ameras}} for model fitting,
  \code{\link[=confint.amerasfit]{confint}} for confidence intervals,
  \code{\link{plot.lm}} for the equivalent method for linear models.
}

\examples{
data("data", package="ameras")
dosevars <- paste0("V", 1:10)

## Binomial model
fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods="RC")

## Both diagnostic plots for RC
plot(fit)

## Residuals vs fitted only
plot(fit, which="residuals-vs-fitted")

## Deviance residuals
plot(fit, type="deviance")

## Multiple methods
\donttest{
fit2 <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="binomial", methods=c("RC", "ERC"))
plot(fit2)
}

## With keep.data=FALSE, supply data explicitly
\donttest{
fit3 <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="binomial", methods="RC",
               keep.data=FALSE)
plot(fit3, data=data)
}

## Schoenfeld residual plot
\donttest{
fit4 <- ameras(Surv(time, status) ~ dose(all_of(dosevars), model = "ERR"),
              data = data, family = "prophaz", methods = "RC")
plot(fit4)
}
}

