\name{residuals.amerasfit}
\alias{residuals.amerasfit}
\title{Residuals for an amerasfit Object}

\description{
  Computes residuals for a fitted \code{amerasfit} object.
}

\usage{
\method{residuals}{amerasfit}(object,
          method = "RC", type = c("pearson", "deviance", "response", "schoenfeld"),
          data = NULL, dose_col = NULL, scaled_schoenfeld = TRUE, \dots)
}

\arguments{
  \item{object}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}.
  }
  \item{method}{
    Character string specifying which estimation method to compute
    residuals for. One of \code{"RC"}, \code{"ERC"}, \code{"MCML"},
    \code{"FMA"}, or \code{"BMA"}. Defaults to \code{"RC"}.
  }
  \item{type}{
    The type of residuals to compute. Let \eqn{\mu_i} denote the
    fitted value for individual \eqn{i}, defined as the predicted
    mean response given the covariates and specified \code{dose_col}.
    For conditional logistic regression, these are the predicted 
    conditional probabilities within matched sets.
    One of:
    \describe{
      \item{\code{"pearson"}}{
        Pearson residuals. For \code{"gaussian"}:
        \eqn{(Y_i - \mu_i) / \sigma} where \eqn{\sigma} is the
        estimated residual standard deviation (note that this differs from
        \code{\link{glm}} which returns the raw (response) residuals). For
        \code{"binomial"} or \code{"clogit"}: \eqn{(Y_i - \mu_i) / \sqrt{\mu_i(1 -
        \mu_i)}}. For \code{"poisson"}: \eqn{(Y_i - \mu_i) /
        \sqrt{\mu_i}}. For \code{"multinomial"}: per-category Pearson
        residuals \eqn{(Y_{iz} - \mu_{iz}) / \sqrt{\mu_{iz}(1 -
        \mu_{iz})}} returned as an \eqn{N \times (Z-1)} matrix, where
        \eqn{\mu_{iz}} is the fitted probability of category \eqn{z}
        for individual \eqn{i}.
      }
      \item{\code{"deviance"}}{
        Deviance residuals, defined as the signed square root of the
        individual contribution to the deviance. For
        \code{"gaussian"}: the raw residual \eqn{Y_i - \mu_i}. For
        \code{"binomial"} or \code{"clogit"}: \eqn{\pm\sqrt{-2\log(\hat{p}_i)}} where
        \eqn{\hat{p}_i = \mu_i} if \eqn{Y_i = 1} and
        \eqn{\hat{p}_i = 1 - \mu_i} if \eqn{Y_i = 0}, with sign
        equal to \eqn{\text{sign}(Y_i - \mu_i)}. For
        \code{"poisson"}: \eqn{\pm\sqrt{2(Y_i \log(Y_i/\mu_i) - (Y_i
        - \mu_i))}} for \eqn{Y_i > 0} and
        \eqn{\pm\sqrt{2\mu_i}} for \eqn{Y_i = 0}, with sign equal to
        \eqn{\text{sign}(Y_i - \mu_i)}. For \code{"multinomial"}:
        per-category deviance residuals returned as an \eqn{N \times
        (Z-1)} matrix.
      }
      \item{\code{"response"}}{
        Raw residuals \eqn{Y_i - \mu_i}. For \code{"multinomial"}
        the matrix \eqn{\bm{Y} - \hat{\bm{P}}} where \eqn{\bm{Y}} is
        the \eqn{N \times Z} indicator matrix of observed categories
        and \eqn{\hat{\bm{P}}} is the \eqn{N \times Z} matrix of
        fitted probabilities.
      }
      \item{\code{"schoenfeld"}}{
        Schoenfeld residuals for proportional hazards-type fits. This is
        a matrix with one column for each covariate in the model, and a row for 
        each event time. For each event time, the vector of raw Schoenfeld residuals 
        is the observed covariate vector for the individual experiencing the event 
        minus the risk-set weighted mean of the covariates at that event time. If
        \code{scaled_schoenfeld = TRUE}, the residuals are scaled using the
        full estimated coefficient covariance matrix, otherwise the raw Schoenfeld 
        residuals are returned. This residual type is only used for
        \code{family = "prophaz"}.
      }
    }
    Defaults to \code{"pearson"} for all families other than \code{"prophaz"}, for which the 
    default is \code{"schoenfeld"}.
  }
  \item{data}{
    The original data frame used for fitting. Only required when the
    model was fitted with \code{keep.data=FALSE}.
  }
  \item{dose_col}{
    Character string specifying the dose column to use when computing
    fitted values. If \code{NULL} (the default), the dose column is
    selected automatically based on the estimation method: the mean
    dose across realizations for RC and ERC, the realization with the
    highest likelihood for MCML and FMA, and the most frequently selected
    realization for BMA. Can be set to any dose column present in the
    data to override the automatic selection.
  }
  \item{scaled_schoenfeld}{
  Logical; whether to compute scaled (\code{scaled_schoenfeld=TRUE}) or raw 
  (\code{scaled_schoenfeld=FALSE}) Schoenfeld residuals. Not used for other 
  residual types.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  For \code{family} in \code{"gaussian"}, \code{"binomial"}, and
  \code{"poisson"}, a numeric vector containing the residuals.
  For \code{family="multinomial"}, a numeric matrix of dimension
  \eqn{N \times (Z-1)} where \eqn{Z} is the number of outcome
  categories and the reference category (the last factor level) is
  excluded. Column names correspond to the non-reference factor levels.
}

\details{
  Fitted values are computed using the mean dose across realizations
  (\code{rcdose_ameras}) for RC and ERC, or the single "best fitting" 
  realization for MCML, FMA, and BMA. For \code{family="multinomial"}, 
  residuals are computed per outcome category.
}

\seealso{
  \code{\link{plot.amerasfit}} for diagnostic plots using these
  residuals,
  \code{\link{ameras}} for model fitting,
  \code{\link[stats]{residuals}} for the generic function.
}

\examples{
data("data", package = "ameras")
dosevars <- paste0("V", 1:10)

fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model = "ERR"),
              data = data, family = "binomial", methods = "RC")

## Pearson residuals (default)
res <- residuals(fit)
summary(res)

## Deviance residuals
res_dev <- residuals(fit, type = "deviance")
summary(res_dev)

## Response residuals
res_raw <- residuals(fit, type = "response")
summary(res_raw)

## Specific method
\donttest{
fit2 <- ameras(Y.binomial ~ dose(all_of(dosevars), model = "ERR"),
               data = data, family = "binomial", methods = c("RC", "ERC"))
res_erc <- residuals(fit2, method = "ERC")
}

## Specific dose column
res_v1 <- residuals(fit, dose_col = "V1")

## With keep.data = FALSE
\donttest{
fit3 <- ameras(Y.binomial ~ dose(all_of(dosevars), model = "ERR"),
               data = data, family = "binomial", methods = "RC",
               keep.data = FALSE)
res <- residuals(fit3, data = data)
}
## Schoenfeld residuals for proportional hazards fits
\donttest{
fit4 <- ameras(Surv(time, status) ~ dose(all_of(dosevars), model = "ERR"),
              data = data, family = "prophaz", methods = "RC")
res_sch <- residuals(fit4, type = "schoenfeld", scaled_schoenfeld = TRUE)
}
}

