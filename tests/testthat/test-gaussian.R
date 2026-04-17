set.seed(123)
data("data", package="ameras")
dosevars <- paste0("V",1:10)

# Gaussian, basic model, all methods 

for(method in all_methods){
  test_that(paste("gaussian snapshot:", method), {
    if(method=="BMA"){
      skip_on_cran()
    }
    
    fit <- fit_combination(
      family    = "gaussian",
      Y         = "Y.gaussian",
      deg       = 2,
      X         = NULL,
      M         = NULL,
      data      = data,
      dosevars  = dosevars,
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