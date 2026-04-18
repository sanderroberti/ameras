set.seed(123)
data("data", package="ameras")


# Clogit, basic model, all methods

for(method in all_methods){
  test_that(paste("clogit snapshot:", method), {
    if(method%in%c("ERC","MCML","BMA")){
      skip_on_cran()
    }
    
    fit <- fit_combination(
      family    = "clogit",
      Y         = "Y.clogit",
      doseRRmod = "EXP",
      deg       = 2,
      X         = NULL,
      M         = NULL,
      data      = data,
      methods   = method,
      niter.BMA = 1000,
      nburnin.BMA = 500
    )
    fit <- confint(fit, type=c("wald.orig","percentile"))
    expect_snapshot(fit[[method]]$coefficients)
    expect_snapshot(fit[[method]]$sd)
    expect_snapshot(fit[[method]]$CI)
  })
}

