# Effect modification

``` r

library(ameras)
data(data, package = "ameras")
```

## Introduction

Effect modification is specified inside `dose()` with the `modifier`
argument. Formula modifiers can be binary numeric/logical variables or
factors. Character variables are not accepted as modifiers; convert them
to factors first so the level order is explicit.

There are two useful ways to code effect modification:

- `modifier = ~ M` uses reference-plus-contrast coding.
- `modifier = ~ 0 + M` or `modifier = ~ M - 1` uses subgroup-specific
  coding.

These are reparameterizations of the same model. They give the same
fitted likelihood when the same model and data are used, but the
coefficients answer different questions directly.

## Reference-plus-contrast coding

Reference-plus-contrast coding is the default formula-modifier coding.
For a binary modifier `M1`, the coefficient named `dose` is the dose
effect in the reference group, and `dose:M1` is the difference between
the non-reference group and the reference group.

``` r

fit_contrast <- ameras(
  Y.gaussian ~ dose(V1:V10, modifier = ~ M1) + X1 + X2,
  data = data,
  family = "gaussian",
  methods = "RC"
)

coef(fit_contrast)
#>                     RC
#> (Intercept) -1.3602199
#> X1           0.4788215
#> X2          -0.5192733
#> dose         1.0836911
#> dose:M1      0.1566173
#> sigma        1.1021197
```

This coding is convenient when the effect-modification contrast is the
main quantity of interest.

## Subgroup-specific coding

Subgroup-specific coding reports one dose-effect parameter per subgroup.
For the same binary modifier `M1`, the coefficient names are
`dose[M1=0]` and `dose[M1=1]`.

``` r

fit_subgroup <- ameras(
  Y.gaussian ~ dose(V1:V10, modifier = ~ 0 + M1) + X1 + X2,
  data = data,
  family = "gaussian",
  methods = "RC"
)

coef(fit_subgroup)
#>                     RC
#> (Intercept) -1.3602200
#> X1           0.4788212
#> X2          -0.5192734
#> dose[M1=0]   1.0836912
#> dose[M1=1]   1.2403085
#> sigma        1.1021195
```

The two codings are equivalent. The non-reference subgroup effect is the
reference effect plus the contrast:

``` r

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
#>                  term from_contrast from_subgroup
#> 1     reference group      1.083691      1.083691
#> 2 non-reference group      1.240308      1.240309
```

Subgroup-specific coding is usually the cleaner choice when
subgroup-specific effect estimates and confidence intervals are the main
goal. For linear ERR and linear-exponential models, it can also be more
numerically robust because default bounds are applied directly to each
subgroup dose-effect parameter. With reference-plus-contrast coding,
contrast parameters are not bounded by the default transformation.

## Factor modifiers

Factor modifiers use the existing level order of the factor. The first
level is the reference level.

``` r

data$M_factor <- factor(
  ifelse(data$M1 == 0, "unexposed modifier group", "modifier group"),
  levels = c("unexposed modifier group", "modifier group")
)
```

With `modifier = ~ M_factor`, ameras reports contrast terms for the
non-reference levels.

``` r

fit_factor_contrast <- ameras(
  Y.gaussian ~ dose(V1:V10, modifier = ~ M_factor) + X1 + X2,
  data = data,
  family = "gaussian",
  methods = "RC"
)

coef(fit_factor_contrast)
#>                                      RC
#> (Intercept)                  -1.3602199
#> X1                            0.4788215
#> X2                           -0.5192733
#> dose                          1.0836911
#> dose:M_factor=modifier group  0.1566173
#> sigma                         1.1021197
```

With `modifier = ~ 0 + M_factor`, ameras reports subgroup-specific
effects for all levels.

``` r

fit_factor_subgroup <- ameras(
  Y.gaussian ~ dose(V1:V10, modifier = ~ 0 + M_factor) + X1 + X2,
  data = data,
  family = "gaussian",
  methods = "RC"
)

coef(fit_factor_subgroup)
#>                                                 RC
#> (Intercept)                             -1.3602200
#> X1                                       0.4788212
#> X2                                      -0.5192734
#> dose[M_factor=unexposed modifier group]  1.0836912
#> dose[M_factor=modifier group]            1.2403085
#> sigma                                    1.1021195
```

Multi-level factors are supported. With reference-plus-contrast coding,
ameras reports one contrast for each non-reference level. With
subgroup-specific coding, ameras reports one dose-effect parameter for
each factor level.

``` r

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
#>                         RC
#> (Intercept)     -1.3619937
#> X1               0.4813357
#> X2              -0.5196477
#> dose[M3=low]     1.1433969
#> dose[M3=middle]  1.1715892
#> dose[M3=high]    1.1807593
#> sigma            1.1073548
```

## Confidence intervals

For subgroup-specific confidence intervals, fit the model using
subgroup-specific coding and then call
[`confint()`](https://rdrr.io/r/stats/confint.html) as usual. For `RC`,
`ERC`, and `MCML`, profile likelihood intervals treat subgroup effects
as ordinary coefficient rows.

``` r

fit_subgroup <- confint(
  fit_subgroup,
  type = "proflik",
  parm = "dose"
)

summary(fit_subgroup)
```

For `FMA` and `BMA`, sample-based confidence or credible intervals are
computed from the sampled coefficients. When subgroup-specific coding is
used, BMA priors and MCMC sampling remain on the internal
reference-plus-contrast scale, but stored samples, summaries,
diagnostics, trace plots, and sample-based intervals are reported on the
subgroup-specific scale.

## Choosing a coding

Use reference-plus-contrast coding when the contrast itself is the main
estimand. Use subgroup-specific coding when subgroup effects are the
main estimands, when subgroup-specific profile likelihood intervals are
desired, or when ERR/LINEXP optimization is sensitive to the unbounded
contrast parameterization.

Testing whether effect modification improves the model is a
model-comparison question. Fit a model without the modifier and a model
with the modifier, then compare the fitted likelihoods using the
appropriate degrees of freedom.
[`dose_lrt()`](https://ameras.sanderroberti.com/reference/dose_lrt.md)
tests dose-related parameters, but its global dose test is not the same
as a test for effect modification.
