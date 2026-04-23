\name{summary.amerasfit}
\alias{summary.amerasfit}
\alias{print.summary.amerasfit}
\title{Summarize an amerasfit object}

\description{
  Produces a summary of a fitted \code{amerasfit} object, including
  parameter estimates, standard errors, and confidence intervals if
  computed via \code{\link[=confint.amerasfit]{confint}}.
}

\usage{
\method{summary}{amerasfit}(object, \dots)

\method{print}{summary.amerasfit}(x, digits = max(3, getOption("digits") - 3), \dots)
}

\arguments{
  \item{object}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}.
  }
  \item{x}{
    An object of class \code{summary.amerasfit}, as returned by
    \code{summary.amerasfit}.
  }
  \item{digits}{
    The number of significant digits to use. Defaults to
    \code{max(3, getOption("digits") - 3)}.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  \code{summary.amerasfit} returns an object of class
  \code{summary.amerasfit}, which is a list containing the following
  elements:
  \describe{
    \item{\code{call}}{
      The matched call from the original \code{\link{ameras}} invocation.
    }
    \item{\code{summary_table}}{
      A data frame with one row per parameter per method, containing
      columns:
      \describe{
        \item{\code{Method}}{The estimation method (RC, ERC, MCML, FMA, or BMA).}
        \item{\code{Term}}{The parameter name.}
        \item{\code{Estimate}}{The parameter estimate.}
        \item{\code{SE}}{The standard error.}
        \item{\code{CI.lower}}{
          The lower confidence bound, if confidence intervals have been
          computed via \code{\link[=confint.amerasfit]{confint}}.
        }
        \item{\code{CI.upper}}{
          The upper confidence bound, if confidence intervals have been
          computed via \code{\link[=confint.amerasfit]{confint}}.
        }
        \item{\code{pval.lower}}{
          p-value associated with the lower bound of the profile likelihood 
          CI, the assess validity of the obtained bound. Only included if
          profile likelihood confidence intervals were computed via
          \code{\link[=confint.amerasfit]{confint}}.
        }
        \item{\code{pval.upper}}{
          p-value associated with the upper bound of the profile likelihood 
          CI, the assess validity of the obtained bound. Only included if
          profile likelihood confidence intervals were computed via
          \code{\link[=confint.amerasfit]{confint}}.
        }
        \item{\code{Rhat}}{
          The Gelman-Rubin convergence diagnostic, included only when
          BMA results are present. Values above 1.05 indicate potential convergence problems.
        }
        \item{\code{n.eff}}{
          The effective sample size, included only when BMA results are present.
        }
      }
    }
    \item{\code{runtime_table}}{
      A data frame with columns \code{Method} and \code{Runtime},
      reporting the computation time in seconds for each method.
    }
    \item{\code{total_runtime_seconds}}{
      The total computation time in seconds across all methods.
    }
    \item{\code{CI.computed}}{
      Logical. \code{TRUE} if confidence intervals have been computed
      via \code{\link[=confint.amerasfit]{confint}}, \code{FALSE}
      otherwise.
    }
  }
}

\details{
  \code{summary.amerasfit} collects results from all estimation methods
  present in the fitted object into a single summary table. Columns for
  confidence intervals are only printed if they have been computed by
  \code{\link[=confint.amerasfit]{confint}}. When BMA results are present 
  in the fitted object, the summary table includes columns \code{Rhat} and 
  \code{n.eff}, with \code{NA} values for all other methods.
}

\seealso{
  \code{\link{ameras}} for model fitting,
  \code{\link[=confint.amerasfit]{confint}} for computing confidence
  intervals,
  \code{\link[=print.amerasfit]{print}} for a shorter printed summary,
  \code{\link[=coef.amerasfit]{coef}} for extracting coefficients.
}

\examples{
data("data", package = "ameras")

## Fit the model
fit <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
              methods = "RC")

## Summary without confidence intervals
summary(fit)

## Summary with confidence intervals
fit <- confint(fit, method = "wald.orig")
summary(fit)

## Access the summary table directly
s <- summary(fit)
s$summary_table

## Multiple methods
\dontrun{
fit2 <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
              methods = c("RC", "ERC", "MCML"))
fit2 <- confint(fit2, method = "wald.orig")
summary(fit2)
}
}
