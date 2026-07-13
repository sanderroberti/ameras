\name{summary.amerasfit}
\alias{summary.amerasfit}
\alias{print.summary.amerasfit}
\alias{summary_table}
\alias{summary_table.amerasfit}
\title{Summarize an amerasfit Object}

\description{
  Produces a summary of a fitted \code{amerasfit} object, including
  parameter estimates, standard errors, and confidence intervals if
  computed via \code{\link[=confint.amerasfit]{confint}}.
}

\usage{
\method{summary}{amerasfit}(object, \dots)

\method{summary_table}{amerasfit}(object, \dots)

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
      The matched call from the original call to \code{\link{ameras}}.
    }
    \item{\code{row_info}}{
      A list describing row handling during model fitting, including the number
      of rows supplied, omitted by \code{na.action}, and used for fitting.
      For conditional logistic regression, this also includes the number of
      rows excluded because they belong to uninformative matched sets. For
      proportional hazards regression with entry times, this also includes the
      number of rows excluded because they have zero follow-up time
      (\code{entry == exit}).
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
      A data frame reporting computation time in seconds for each method.
      When structured timing information is available this has columns
      \code{Method}, \code{Fit}, \code{CI}, and \code{Total}, with CPU time
      reported for each timing component. For objects that only contain the
      compatibility \code{runtime} field, this has columns \code{Method} and
      \code{Runtime}.
    }
    \item{\code{total_runtime_seconds}}{
      The total CPU computation time in seconds across all methods.
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
  \code{n.eff}, with \code{NA} values for all other methods. Profile
  likelihood diagnostic columns such as \code{pval.lower} and
  \code{pval.upper} remain available in the stored \code{CI} component
  but are not printed in the summary table.

  Printing the summary prints the original call to \code{ameras}, row counts,
  the CPU runtime (total and by method), and the table described above. Row
  counts are shown in the order applied during fitting. The supplied row count
  is always shown. Omitted or excluded row counts are shown only when nonzero,
  and the rows used for fitting are shown when this differs from the supplied
  row count. For conditional logistic regression, the summary reports rows
  excluded because they belong to uninformative matched sets, i.e., matched sets
  of size 1 or with no cases. For proportional hazards regression with entry
  times, the summary reports rows excluded because they have zero follow-up time
  (\code{entry == exit}). When confidence intervals have been computed, CI
  runtime is included in the total runtime.
  This table can also be accessed directly (i.e., to retrieve confidence
  intervals) using \code{summary_table}.
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
fit <- confint(fit, type = "wald.orig")
summary(fit)

## Access the summary table directly
s <- summary_table(fit)

## Multiple methods
\donttest{
fit2 <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
              methods = c("RC", "ERC", "MCML"))
fit2 <- confint(fit2, type = "wald.orig")
summary(fit2)
}
}
