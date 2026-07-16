## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE
)

## ----setup--------------------------------------------------------------------
library(ameras)
data(data, package = "ameras")

## ----contrast-fit-------------------------------------------------------------
fit_contrast <- ameras(
  Y.gaussian ~ dose(V1:V10, modifier = ~ M1) + X1 + X2,
  data = data,
  family = "gaussian",
  methods = "RC"
)

coef(fit_contrast)

## ----subgroup-fit-------------------------------------------------------------
fit_subgroup <- ameras(
  Y.gaussian ~ dose(V1:V10, modifier = ~ 0 + M1) + X1 + X2,
  data = data,
  family = "gaussian",
  methods = "RC"
)

coef(fit_subgroup)

## ----coding-equivalence-------------------------------------------------------
contrast_coef <- fit_contrast$RC$coefficients
subgroup_coef <- fit_subgroup$RC$coefficients

data.frame(
  term = c("reference group", "non-reference group"),
  from_contrast = c(
    contrast_coef["dose"],
    contrast_coef["dose"] + contrast_coef["dose:M1"]
  ),
  from_subgroup = c(
    subgroup_coef["dose[M1=0]"],
    subgroup_coef["dose[M1=1]"]
  ),
  row.names = NULL
)

## ----factor-setup-------------------------------------------------------------
data$M_factor <- factor(
  ifelse(data$M1 == 0, "unexposed modifier group", "modifier group"),
  levels = c("unexposed modifier group", "modifier group")
)

## ----factor-contrast----------------------------------------------------------
fit_factor_contrast <- ameras(
  Y.gaussian ~ dose(V1:V10, modifier = ~ M_factor) + X1 + X2,
  data = data,
  family = "gaussian",
  methods = "RC"
)

coef(fit_factor_contrast)

## ----factor-subgroup----------------------------------------------------------
fit_factor_subgroup <- ameras(
  Y.gaussian ~ dose(V1:V10, modifier = ~ 0 + M_factor) + X1 + X2,
  data = data,
  family = "gaussian",
  methods = "RC"
)

coef(fit_factor_subgroup)

## ----three-level-factor-------------------------------------------------------
data$M3 <- factor(
  rep(c("low", "middle", "high"), length.out = nrow(data)),
  levels = c("low", "middle", "high")
)

fit_three_level <- ameras(
  Y.gaussian ~ dose(V1:V10, modifier = ~ 0 + M3) + X1 + X2,
  data = data,
  family = "gaussian",
  methods = "RC"
)

coef(fit_three_level)

## ----subgroup-ci, eval = FALSE------------------------------------------------
# fit_subgroup <- confint(
#   fit_subgroup,
#   type = "proflik",
#   parm = "dose"
# )
# 
# summary(fit_subgroup)

