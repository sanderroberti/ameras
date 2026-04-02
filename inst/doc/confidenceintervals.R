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
fit.ameras.waldorig <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", methods=c("RC"), CI="wald.orig", doseRRmod="ERR")
fit.ameras.waldtransformed <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), 
                                     data=data, family="binomial", methods=c("RC"), 
                                     CI="wald.transformed", doseRRmod="ERR")
fit.ameras.proflik <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", methods=c("RC"), CI="proflik", doseRRmod="ERR", 
                            params.profCI = "all")
summary(fit.ameras.waldorig)
summary(fit.ameras.waldtransformed)
summary(fit.ameras.proflik)

## ----fits2, eval = identical(Sys.getenv("NOT_CRAN"), "true")------------------
fit.ameras.hpd <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", methods=c("FMA"), CI="hpd", doseRRmod="ERR")
fit.ameras.percentile <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", methods=c("FMA"), CI="percentile", 
                            doseRRmod="ERR")

summary(fit.ameras.hpd)
summary(fit.ameras.percentile)

