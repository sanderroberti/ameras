## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ameras)


## ----fits1, eval = identical(Sys.getenv("NOT_CRAN"), "true")------------------
data(data, package="ameras")
dosevars <- paste0("V", 1:10)
fit.ameras <- ameras(Y.binomial~dose(V1:V10, model="ERR")+X1+X2, data=data, 
                            family="binomial", methods=c("RC"))

fit.ameras.waldorig <- confint(fit.ameras, type="wald.orig")
fit.ameras.waldtransformed <- confint(fit.ameras, type="wald.transformed")
fit.ameras.proflik <- confint(fit.ameras, type="proflik", parm="all")

summary(fit.ameras.waldorig)
summary(fit.ameras.waldtransformed)
summary(fit.ameras.proflik)

## ----fits2, eval = identical(Sys.getenv("NOT_CRAN"), "true")------------------
fit.ameras2 <- ameras(Y.binomial~dose(V1:V10, model="ERR")+X1+X2, data=data, 
                            family="binomial", methods=c("FMA"))

fit.ameras.hpd <- confint(fit.ameras2, type="hpd")
fit.ameras.percentile <- confint(fit.ameras2, type="percentile")

summary(fit.ameras.hpd)
summary(fit.ameras.percentile)

