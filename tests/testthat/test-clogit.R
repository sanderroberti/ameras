set.seed(123)
data("data", package="ameras")
dosevars <- paste0("V",1:10)

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
      dosevars  = dosevars,
      methods   = method,
      CI        = c("wald.orig", "percentile"),
      niter.BMA = 1000,
      nburnin.BMA = 500
    )
    expect_snapshot(fit[[method]]$coefficients)
    expect_snapshot(fit[[method]]$sd)
    expect_snapshot(fit[[method]]$CI)
  })
}

