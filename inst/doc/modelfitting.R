## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ameras)

## -----------------------------------------------------------------------------
data(data, package="ameras")
head(data)

## ----modelfit.linreg, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(12345)
fit.ameras.linreg <- ameras(Y.gaussian~dose(V1:V10)+X1+X2, data=data, family="gaussian", 
                            methods=c("RC", "ERC", "MCML", "FMA", "BMA"), 
                            niter.BMA = 5000, nburnin.BMA = 1000)

## ----eval = identical(Sys.getenv("NOT_CRAN"), "true")-------------------------
str(fit.ameras.linreg)

## ----eval = identical(Sys.getenv("NOT_CRAN"), "true")-------------------------
fit.ameras.linreg$RC

## ----eval = identical(Sys.getenv("NOT_CRAN"), "true")-------------------------
summary(fit.ameras.linreg)

## ----eval = identical(Sys.getenv("NOT_CRAN"), "true")-------------------------
coef(fit.ameras.linreg)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.linreg)

## ----modelfit.logreg, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(33521)
fit.ameras.logreg <- ameras(Y.binomial~dose(V1:V10, deg=2, model="EXP")+X1+X2, data=data, 
                            family="binomial", methods=c("RC", "ERC", "MCML", "FMA", "BMA"), 
                            niter.BMA = 5000, nburnin.BMA = 1000)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.logreg)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.logreg)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.logreg)

## ----modelfit.logreg.lin, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(3521216)
fit.ameras.logreg.lin <- ameras(Y.binomial~dose(V1:V10, deg=1, model="EXP")+X1+X2,  data=data, 
                                family="binomial", methods=c("RC", "ERC", "MCML", "FMA", "BMA"), 
                                niter.BMA = 5000, nburnin.BMA = 1000)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.logreg.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.logreg.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.logreg.lin)

## ----modelfit.poisson, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(332101)
fit.ameras.poisson <- ameras(Y.poisson~dose(V1:V10, deg=2, model="EXP")+X1+X2, data=data, 
                             family="poisson", methods=c("RC", "ERC", "MCML", "FMA", "BMA"), 
                             niter.BMA = 5000, nburnin.BMA = 1000)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.poisson)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.poisson)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.poisson)

## ----modelfit.poisson.lin, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(24252)
fit.ameras.poisson.lin <- ameras(Y.poisson~dose(V1:V10, deg=1, model="EXP")+X1+X2, data=data, 
                                 family="poisson", methods=c("RC", "ERC", "MCML", "FMA", "BMA"), 
                                 niter.BMA = 5000, nburnin.BMA = 1000)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.poisson.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.poisson.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.poisson.lin)

## ----modelfit.prophaz, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(332120000)
fit.ameras.prophaz <- ameras(Surv(time, status)~dose(V1:V10, deg=2, model="EXP")+X1+X2, 
                             data=data, family="prophaz", 
                             methods=c("RC", "ERC", "MCML", "FMA", "BMA"), niter.BMA = 5000, 
                             nburnin.BMA = 1000)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.prophaz)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.prophaz)

## ----eval = identical(Sys.getenv("NOT_CRAN"), "true")-------------------------
fit.ameras.prophaz$BMA$prophaz.timepoints

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.prophaz)

## ----modelfit.prophaz.lin, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(24978252)
fit.ameras.prophaz.lin <- ameras(Surv(time, status)~dose(V1:V10, deg=1, model="EXP")+X1+X2, 
                             data=data, family="prophaz", 
                             methods=c("RC", "ERC", "MCML", "FMA", "BMA"), niter.BMA = 5000, 
                             nburnin.BMA = 1000)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.prophaz.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.prophaz.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.prophaz.lin)


## ----modelfit.multinomial, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(33)
fit.ameras.multinomial <- ameras(Y.multinomial~dose(V1:V10, deg=2, model="EXP")+X1+X2, data=data, 
                            family="multinomial",methods=c("RC", "ERC", "MCML", "FMA", "BMA"), 
                            niter.BMA = 5000, nburnin.BMA = 1000)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.multinomial)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.multinomial)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.multinomial)

## ----modelfit.multinomial.lin, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(44)
fit.ameras.multinomial.lin <- ameras(Y.multinomial~dose(V1:V10, deg=1, model="EXP")+X1+X2, data=data, 
                            family="multinomial",methods=c("RC", "ERC", "MCML", "FMA", "BMA"), 
                            niter.BMA = 5000, nburnin.BMA = 1000)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.multinomial.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.multinomial.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.multinomial.lin)

## ----modelfit.clogit, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(3301)
fit.ameras.clogit <- ameras(Y.clogit~dose(V1:V10, deg=2, model="EXP")+X1+X2+strata(setnr), data=data, 
                            family="clogit", methods=c("RC", "ERC", "MCML", "FMA", "BMA"), 
                            niter.BMA = 5000, nburnin.BMA = 1000)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.clogit)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.clogit)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.clogit)

## ----modelfit.clogit.lin, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(4401)
fit.ameras.clogit.lin <- ameras(Y.clogit~dose(V1:V10, deg=2, model="EXP")+X1+X2+strata(setnr), data=data, 
                            family="clogit", methods=c("RC", "ERC", "MCML", "FMA", "BMA"), 
                            niter.BMA = 5000, nburnin.BMA = 1000)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.clogit.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.clogit.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.clogit.lin)

