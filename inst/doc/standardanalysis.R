## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE
)

## ----setup--------------------------------------------------------------------
library(ameras)

## ----data---------------------------------------------------------------------
data(data, package = "ameras")

## Use one exposure realization as the single observed dose.
data$D <- data$V1
head(data[c("D", "Y.gaussian", "Y.binomial", "Y.poisson", "X1", "X2")])

## ----helper-------------------------------------------------------------------
compare_coefficients <- function(ameras_fit, standard_fit) {
  ameras_terms <- c("(Intercept)", "X1", "X2", "dose")
  standard_terms <- c("(Intercept)", "X1", "X2", "D")

  out <- data.frame(
    term = ameras_terms,
    ameras = unname(ameras_fit$RC$coefficients[ameras_terms]),
    standard = unname(stats::coef(standard_fit)[standard_terms])
  )

  out$ameras <- round(out$ameras, 4)
  out$standard <- round(out$standard, 4)
  out
}

## ----gaussian-----------------------------------------------------------------
fit_ameras_gaussian <- ameras(Y.gaussian ~ dose(D) + X1 + X2,
                              data = data,
                              family = "gaussian")

fit_lm <- lm(Y.gaussian ~ D + X1 + X2, data = data)

compare_coefficients(fit_ameras_gaussian, fit_lm)

## ----gaussian-summary---------------------------------------------------------
fit_ameras_gaussian <- confint(fit_ameras_gaussian,
                               type = "wald.orig",
                               print = FALSE)
summary(fit_ameras_gaussian)

## ----binomial-----------------------------------------------------------------
fit_ameras_binomial <- ameras(Y.binomial ~ dose(D, model = "EXP") + X1 + X2,
                              data = data,
                              family = "binomial")

fit_glm_binomial <- glm(Y.binomial ~ D + X1 + X2,
                        data = data,
                        family = binomial())

compare_coefficients(fit_ameras_binomial, fit_glm_binomial)

## ----poisson------------------------------------------------------------------
fit_ameras_poisson <- ameras(Y.poisson ~ dose(D, model = "EXP") + X1 + X2,
                             data = data,
                             family = "poisson")

fit_glm_poisson <- glm(Y.poisson ~ D + X1 + X2,
                       data = data,
                       family = poisson())

compare_coefficients(fit_ameras_poisson, fit_glm_poisson)

## ----multiple-realizations----------------------------------------------------
fit_ameras_rc <- ameras(Y.binomial ~ dose(V1:V10, model = "EXP") + X1 + X2,
                        data = data,
                        family = "binomial")
summary(fit_ameras_rc)

## ----multiple-methods, eval = FALSE-------------------------------------------
# fit_ameras_all <- ameras(Y.binomial ~ dose(V1:V10, model = "EXP") + X1 + X2,
#                          data = data,
#                          family = "binomial",
#                          methods = c("RC", "ERC", "MCML", "FMA", "BMA"))

