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

## -----------------------------------------------------------------------------
dosevars <- paste0("V", 1:10)

## ----modelfit.linreg, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(12345)
fit.ameras.linreg <- ameras(Y="Y.gaussian", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="gaussian", methods=c("RC", "ERC", "MCML", "FMA", "BMA"), 
                            niter.BMA = 5000, nburnin.BMA = 1000, CI=c("wald.orig","percentile"))

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
fit.ameras.logreg <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", deg=2, doseRRmod = "EXP", 
                            methods=c("RC", "ERC", "MCML", "FMA", "BMA"), niter.BMA = 5000, 
                            nburnin.BMA = 1000, CI=c("wald.orig","percentile"))

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.logreg)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.logreg)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.logreg)

## ----modelfit.logreg.lin, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(3521216)
fit.ameras.logreg.lin <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                                family="binomial", deg=1, doseRRmod = "EXP", 
                                methods=c("RC", "ERC", "MCML", "FMA", "BMA"), niter.BMA = 5000, 
                                nburnin.BMA = 1000, CI=c("wald.orig","percentile"))

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.logreg.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.logreg.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.logreg.lin)

## ----modelfit.poisson, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(332101)
fit.ameras.poisson <- ameras(Y="Y.poisson", dosevars=dosevars, X=c("X1","X2"), data=data, 
                             family="poisson", deg=2, doseRRmod = "EXP", 
                             methods=c("RC", "ERC", "MCML", "FMA", "BMA"), niter.BMA = 5000, 
                             nburnin.BMA = 1000, CI=c("wald.orig","percentile"))

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.poisson)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.poisson)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.poisson)

## ----modelfit.poisson.lin, cache=TRUE, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
set.seed(24252)
fit.ameras.poisson.lin <- ameras(Y="Y.poisson", dosevars=dosevars, X=c("X1","X2"), data=data, 
                                 family="poisson", deg=1, doseRRmod = "EXP", 
                                 methods=c("RC", "ERC", "MCML", "FMA", "BMA"), niter.BMA = 5000, 
                                 nburnin.BMA = 1000, CI=c("wald.orig","percentile"))

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
summary(fit.ameras.poisson.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
coef(fit.ameras.poisson.lin)

## ----fig.fullwidth=TRUE, fig.show="hold", out.width='100%', fig.width=6, fig.height=8, eval = identical(Sys.getenv("NOT_CRAN"), "true")----
traceplot(fit.ameras.poisson.lin)

