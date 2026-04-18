set.seed(123)
data("data", package="ameras")



for(method in c("RC","ERC","MCML")){
  test_that(paste("proflik/wald.transformed snapshot:", method), {
    if(method%in%c("ERC","MCML")){
      skip_on_cran()
    }
    
    fit <- fit_combination(
      family    = "binomial",
      Y         = "Y.binomial",
      doseRRmod = "ERR",
      deg       = 2,
      X         = NULL,
      M         = NULL,
      data      = data,
      methods   = method
    )
    fit1 <- confint(fit, type=c("proflik"))
    expect_snapshot(fit1[[method]]$coefficients)
    expect_snapshot(fit1[[method]]$sd)
    expect_snapshot(fit1[[method]]$CI)
    
    fit2 <- confint(fit, type=c("wald.transformed"))
    expect_snapshot(fit2[[method]]$coefficients)
    expect_snapshot(fit2[[method]]$sd)
    expect_snapshot(fit2[[method]]$CI)
  })
}


for(method in c("FMA","BMA")){
  test_that(paste("percentile/hpd snapshot:", method), {
    skip_on_cran()
    
    fit <- fit_combination(
      family    = "binomial",
      Y         = "Y.binomial",
      doseRRmod = "ERR",
      deg       = 2,
      X         = NULL,
      M         = NULL,
      data      = data,
      methods   = method
    )
    fit1 <- confint(fit, type=c("percentile"))
    expect_snapshot(fit1[[method]]$coefficients)
    expect_snapshot(fit1[[method]]$sd)
    expect_snapshot(fit1[[method]]$CI)
    
    fit2 <- confint(fit, type=c("hpd"))
    expect_snapshot(fit2[[method]]$coefficients)
    expect_snapshot(fit2[[method]]$sd)
    expect_snapshot(fit2[[method]]$CI)
  })
}