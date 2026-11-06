## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(ameras)



## ----fits1, eval = identical(Sys.getenv("NOT_CRAN"), "true")------------------
dosevars <- paste0("V", 1:10)
fit.ameras <- ameras(Y.binomial~dose(V1:V10, model="ERR")+X1+X2, data = ameras_simulated,
                            family="binomial", methods=c("RC"))

ci.ameras.waldorig <- confint(fit.ameras, type="wald.orig")
ci.ameras.waldtransformed <- confint(fit.ameras, type="wald.transformed")
ci.ameras.proflik <- confint(fit.ameras, type="proflik", parm="all")

ci.ameras.waldorig
ci.ameras.waldtransformed
ci.ameras.proflik

# Attach intervals when they should be included in the fit summary
fit.ameras.waldorig <- add_confint(
  fit.ameras,
  type="wald.orig",
  print=FALSE
)
summary(fit.ameras.waldorig)


## ----fits2, eval = identical(Sys.getenv("NOT_CRAN"), "true")------------------
set.seed(123)
fit.ameras2 <- ameras(Y.binomial~dose(V1:V10, model="ERR")+X1+X2, data = ameras_simulated,
                            family="binomial", methods=c("FMA"))

ci.ameras.hpd <- confint(fit.ameras2, type="hpd")
ci.ameras.percentile <- confint(fit.ameras2, type="percentile")

ci.ameras.hpd
ci.ameras.percentile

