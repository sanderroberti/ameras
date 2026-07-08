\name{residuals.amerasfit}
\alias{residuals.amerasfit}
\title{Residuals for an amerasfit Object}

\description{
  Computes residuals for a fitted \code{amerasfit} object.
}

\usage{
\method{residuals}{amerasfit}(object,
          method = "RC",
          type = NULL,
          data = NULL,
          dose.col = NULL,
          scaled.schoenfeld = TRUE,
          \dots)
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
    The type of residuals to compute. Defaults to \code{"pearson"} for
    all families except \code{"prophaz"}, for which the default is
    \code{"schoenfeld"}. For types other than \code{"schoenfeld"}, let \eqn{Y_i} 
    and \eqn{\mu_i} denote the observed response and fitted mean response 
    for individual \eqn{i}, respectively. Then the residuals are defined as follows:
    \describe{
  \item{\code{"pearson"}}{
    Pearson residuals. For \code{"gaussian"}:
    \eqn{(Y_i - \mu_i) / \sigma} where \eqn{\sigma} is the
    estimated residual standard deviation (note: \code{\link{glm}}
    returns raw residuals for Gaussian). For \code{"binomial"} and \code{"clogit"}:
    \eqn{(Y_i - \mu_i) / \sqrt{\mu_i(1 - \mu_i)}}. For
    \code{"poisson"}: \eqn{(Y_i - \mu_i) / \sqrt{\mu_i}}. For
    \code{"multinomial"}: per-category Pearson
    residuals \eqn{(Y_{iz} - \mu_{iz}) / \sqrt{\mu_{iz}(1 -
    \mu_{iz})}} where
    \eqn{\mu_{iz}} is the fitted probability of category \eqn{z}
    for individual \eqn{i}.
    }
    \item{\code{"deviance"}}{
    Deviance residuals, defined as the signed square root of the
    individual contribution to the deviance. For
    \code{"gaussian"}: the raw residual \eqn{Y_i - \mu_i}. For
    \code{"binomial"} and \code{"clogit"}:
    \eqn{\pm\sqrt{-2\log(\hat{p}_i)}} where
    \eqn{\hat{p}_i = \mu_i} if \eqn{Y_i = 1} and
    \eqn{\hat{p}_i = 1 - \mu_i} if \eqn{Y_i = 0}, with sign
    equal to \eqn{\text{sign}(Y_i - \mu_i)}. For \code{"poisson"}:
    \eqn{\pm\sqrt{2(Y_i \log(Y_i/\mu_i) - (Y_i
    - \mu_i))}} for \eqn{Y_i > 0} and \eqn{\pm\sqrt{2\mu_i}} for
    \eqn{Y_i = 0}, with sign equal to
    \eqn{\text{sign}(Y_i - \mu_i)}. For \code{"multinomial"}:
    per-category deviance residuals as for \code{"binomial"}.
    }
      \item{\code{"response"}}{
        Raw residuals \eqn{Y_i - \mu_i}. For \code{"multinomial"}:
        the matrix \eqn{\bm{Y} - \hat{\bm{P}}} where \eqn{\bm{Y}} is
        the \eqn{N \times Z} indicator matrix of observed categories
        and \eqn{\hat{\bm{P}}} is the \eqn{N \times Z} matrix of
        fitted probabilities.
      }
      \item{\code{"schoenfeld"}}{
        Schoenfeld residuals for \code{family="prophaz"} only. For
        each event time, the unscaled residual for the individual 
        experiencing the event is the observed covariate vector 
        minus the risk-set weighted mean covariate vector at that 
        event time. See Details.
      }
    }
  }
  \item{data}{
    The original data frame used for fitting. Only required when the
    model was fitted with \code{keep.data=FALSE}.
  }
  \item{dose.col}{
    Character string specifying the dose column to use when computing
    fitted values. If \code{NULL} (the default), the dose column is
    selected automatically: the mean dose across realizations
    for RC and ERC, the realization with the
    highest likelihood for MCML and FMA, and the most frequently
    selected realization for BMA. Can be set to any dose column
    present in the data to override the automatic selection.
  }
  \item{scaled.schoenfeld}{
    Logical. If \code{TRUE} (the default), scaled Schoenfeld residuals
    are returned, obtained by multiplying the raw Schoenfeld residuals 
    by the estimated coefficient variance-covariance matrix following 
    Grambsch and Therneau (1994). If \code{FALSE}, raw Schoenfeld 
    residuals are returned. Only used when \code{type="schoenfeld"}.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  For families \code{"gaussian"}, \code{"binomial"},
  \code{"poisson"}, and \code{"clogit"}, a numeric vector of 
  length \eqn{N} containing the residuals, where \eqn{N} is the number
  of rows used for fitting after any \code{na.action} handling.

  For \code{family="multinomial"}, a numeric matrix of dimension
  \eqn{N \times Z} where \eqn{N} is the number of rows used for fitting
  after any \code{na.action} handling and \eqn{Z} is the number of outcome
  categories. Column names correspond to the factor levels.
  Note that when plotting via \code{\link[=plot.amerasfit]{plot}}, the 
  reference category is excluded since its residuals are a 
  linear combination of the other categories.

  For \code{family="prophaz"} with \code{type="schoenfeld"}, a data
  frame with columns \code{id} (row index of the event), \code{time}
  (event time), and one additional column per model parameter
  containing the Schoenfeld residuals. Only rows corresponding to
  events (\code{status=1}) are included.
}

\details{
  Fitted values are computed using the dose column specified by
  \code{dose.col}.

  If the model was fitted with \code{na.exclude}, residuals are returned for
  the fitted rows only; they are not padded with \code{NA} values back to the
  originally supplied row count. For proportional hazards models,
  Schoenfeld residuals are event-level quantities and are returned only for
  event rows.

  For \code{family="clogit"}, fitted values \eqn{\mu_i} are the
  conditional probabilities of being a case within each matched set,
  computed as \eqn{R_i / \sum_{j \in \text{set}(i)} R_j} where
  \eqn{R_i} is the relative risk for individual \eqn{i} and the sum
  is over all individuals in the same matched set.

  For \code{family="prophaz"} and \code{type="schoenfeld"}, the
  Schoenfeld residual for the individual failing at time \eqn{t} is:
  \deqn{\bm{r} = \bm{x} - \bar{\bm{x}}(t)}
  where \eqn{\bm{x}} is the observed covariate vector for the failing
  individual and
  \deqn{\bar{\bm{x}}(t) = \frac{\sum_{j \in \mathcal{R}(t)} \bm{x}_j R_j}
  {\sum_{j \in \mathcal{R}(t)} R_j}}
  is the weighted mean covariate vector in the risk set
  \eqn{\mathcal{R}(t)} at time \eqn{t}, with weights equal to the
  relative risks \eqn{R_j}. Ties in event times are handled using the
  Breslow approximation.
}

\references{
  Schoenfeld, D. (1982). Partial residuals for the proportional
  hazards regression model.
  \emph{Biometrika}, \bold{69}(1), 239--241.
  \doi{10.1093/biomet/69.1.239}

  Grambsch, P. M. and Therneau, T. M. (1994). Proportional hazards
  tests and diagnostics based on weighted residuals.
  \emph{Biometrika}, \bold{81}(3), 515--526.
  \doi{10.1093/biomet/81.3.515}
}

\seealso{
  \code{\link[=plot.amerasfit]{plot}} for diagnostic plots using these
  residuals,
  \code{\link{ameras}} for model fitting,
  \code{\link[stats]{residuals}} for the generic function.
}

\examples{
data("data", package="ameras")
dosevars <- paste0("V", 1:10)

fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods="RC")

## Pearson residuals (default)
res <- residuals(fit)
summary(res)

## Deviance residuals
res_dev <- residuals(fit, type="deviance")

## Response residuals
res_raw <- residuals(fit, type="response")

## Specific dose column
res_v1 <- residuals(fit, dose.col="V1")

## Multiple methods
\donttest{
fit2 <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="binomial", methods=c("RC", "ERC"))
res_erc <- residuals(fit2, method="ERC")
}

## With keep.data=FALSE
\donttest{
fit3 <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="binomial", methods="RC",
               keep.data=FALSE)
res <- residuals(fit3, data=data)
}

## Schoenfeld residuals for proportional hazards model
\donttest{
fit4 <- ameras(Surv(time, status) ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="prophaz", methods="RC")
res_sch <- residuals(fit4)
res_raw_sch <- residuals(fit4, scaled.schoenfeld=FALSE)
}
}
