## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ameras)
library(ggplot2)
data(data, package="ameras")

## ----modelfit.exp-------------------------------------------------------------
fit.ameras.exp <- ameras(Y.binomial~dose(V1:V10, deg=2, model="EXP")+X1+X2, 
                         data=data, family="binomial", methods="RC")
summary(fit.ameras.exp)

## ----modelfit.err-------------------------------------------------------------
fit.ameras.err <- ameras(Y.binomial~dose(V1:V10, deg=2, model="ERR")+X1+X2, 
                         data=data, family="binomial", methods="RC")
summary(fit.ameras.err)

## ----modelfit.linexp----------------------------------------------------------
fit.ameras.linexp <- ameras(Y.binomial~dose(V1:V10, model="LINEXP")+X1+X2, 
                         data=data, family="binomial", methods="RC")
summary(fit.ameras.linexp)

## ----comparison, fig.width=7, fig.height=6------------------------------------
ggplot(data.frame(x=c(0, 5)), aes(x))+
  theme_light(base_size=15)+
  xlab("Exposure")+
  ylab("Relative risk")+
  labs(col="Model", lty="Model") +
  theme(legend.position = "inside", 
        legend.position.inside = c(.2,.85),
        legend.box.background = element_rect(color = "black", fill = "white", linewidth = 1))+
  stat_function(aes(col="Linear-quadratic ERR", lty="Linear-quadratic ERR" ),fun=function(x){
    1+fit.ameras.err$RC$coefficients["dose"]*x + fit.ameras.err$RC$coefficients["dose_squared"]*x^2
  }, linewidth=1.2) + 
  stat_function(aes(col="Exponential", lty="Exponential"),fun=function(x){
    exp(fit.ameras.exp$RC$coefficients["dose"]*x + fit.ameras.exp$RC$coefficients["dose_squared"]*x^2)
  }, linewidth=1.2) +
  stat_function(aes(col="Linear-exponential", lty="Linear-exponential"),fun=function(x){
    1+fit.ameras.linexp$RC$coefficients["dose_linear"]*x * exp(fit.ameras.linexp$RC$coefficients["dose_exponential"]*x)
  }, linewidth=1.2)



