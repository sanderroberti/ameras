set.seed(123)
data("data", package="ameras")
dosevars <- paste0("V", 1:10)


# Snapshot test 
binomial_combos <- Filter(
  function(x) x$family == "binomial", 
  snapshot_combinations
)

for (combo in binomial_combos) {
  
  label <- sprintf("%s_%s_deg%d",
                   combo$family,
                   combo$doseRRmod,
                   combo$deg)
  
  test_that(paste("snapshot:", label), {
    fit <- fit_combination(
      family    = combo$family,
      Y         = combo$Y,
      doseRRmod = combo$doseRRmod,
      deg       = combo$deg,
      X         = "X1",
      M         = NULL,
      data      = data,
      dosevars  = dosevars
    )
    fit <- confint(fit, type=c("wald.orig","percentile"))
    expect_snapshot(fit$RC$coefficients)
    expect_snapshot(fit$RC$sd)
    expect_snapshot(fit$RC$CI)
  })
}


# Test all non-RC methods with snapshot for a basic model
for(method in setdiff(all_methods,"RC")){
  test_that(paste("binomial snapshot:", method), {
    if(method%in%c("ERC","MCML","BMA")){
      skip_on_cran()
    }
    
    fit <- fit_combination(
      family    = "binomial",
      Y         = "Y.binomial",
      doseRRmod = "EXP",
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


# Basic no-error check for RC and all combinations of doseRRmod, deg, and lengths of X and M

# for(combo in binomial_combos){
#   for(cov_combo in covariate_combinations){
#     X_label <- if(is.null(cov_combo$X)) "NULL" else paste(cov_combo$X, collapse="-")
#     M_label <- if(is.null(cov_combo$M)) "NULL" else paste(cov_combo$M, collapse="-")
#     label <- sprintf("%s_%s_deg%d_X%s_M%s",
#                      combo$family,
#                      combo$doseRRmod,
#                      combo$deg,
#                      X_label,
#                      M_label)
#     
#     test_that(label, {
#       
#       expect_no_error({
#         fit_combination(
#           family    = combo$family,
#           Y         = combo$Y,
#           deg       = combo$deg,
#           doseRRmod = combo$doseRRmod,
#           X         = cov_combo$X,
#           M         = cov_combo$M,
#           data      = data,
#           dosevars  = dosevars,
#           methods   = "RC",
#           CI        = "wald.orig"
#         )
#       })
#     })
#   }
# }

# Test RC with all 45 combinations of doseRRmod, deg and lengths of X and M

test_that("ameras binomial ERR deg=1 methods=RC X=c(X1,X2) M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1, methods="RC", X=c("X1","X2"), M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:6))
})

test_that("ameras binomial ERR deg=2 methods=RC X=c(X1,X2) M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=2, methods="RC", X=c("X1","X2"), M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:9))
})

test_that("ameras binomial ERR deg=1 methods=RC X=X1 M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=1, methods="RC", X=c("X1"), M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:5))
})

test_that("ameras binomial ERR deg=2 methods=RC X=X1 M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=2, methods="RC", X=c("X1"), M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:8))
})

test_that("ameras binomial ERR deg=1 methods=RC X=NULL M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=1, methods="RC", X=NULL, M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:4))
})

test_that("ameras binomial ERR deg=2 methods=RC X=NULL M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=2, methods="RC", X=NULL, M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:7))
})






test_that("ameras binomial ERR deg=1 methods=RC X=c(X1,X2) M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=1, methods="RC", X=c("X1","X2"), M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:5))
})

test_that("ameras binomial ERR deg=2 methods=RC X=c(X1,X2) M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=2, methods="RC", X=c("X1","X2"), M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:7))
})

test_that("ameras binomial ERR deg=1 methods=RC X=X1 M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=1, methods="RC", X=c("X1"), M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:4))
})

test_that("ameras binomial ERR deg=2 methods=RC X=X1 M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=2, methods="RC", X=c("X1"), M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:6))
})

test_that("ameras binomial ERR deg=1 methods=RC X=NULL M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=1, methods="RC", X=NULL, M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:3))
})

test_that("ameras binomial ERR deg=2 methods=RC X=NULL M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=2, methods="RC", X=NULL, M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:5))
})






test_that("ameras binomial ERR deg=1 methods=RC X=c(X1,X2) M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=1, methods="RC", X=c("X1","X2"), M=NULL, doseRRmod="ERR"))
})

test_that("ameras binomial ERR deg=2 methods=RC X=c(X1,X2) M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=2, methods="RC", X=c("X1","X2"), M=NULL, doseRRmod="ERR"))
})

test_that("ameras binomial ERR deg=1 methods=RC X=X1 M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=1, methods="RC", X=c("X1"), M=NULL, doseRRmod="ERR"))
})

test_that("ameras binomial ERR deg=2 methods=RC X=X1 M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=2, methods="RC", X=c("X1"), M=NULL, doseRRmod="ERR"))
})

test_that("ameras binomial ERR deg=1 methods=RC X=NULL M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",  dosevars=dosevars, deg=1, methods="RC", X=NULL, M=NULL, doseRRmod="ERR"))
})

test_that("ameras binomial ERR deg=2 methods=RC X=NULL M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2, methods="RC", X=NULL, M=NULL, doseRRmod="ERR"))
})





#------------------------


test_that("ameras binomial LINEXP   methods=RC X=c(X1,X2) M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars,   methods="RC", X=c("X1","X2"), M=c("M1", "M2"), doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:9))
})

test_that("ameras binomial LINEXP   methods=RC X=X1 M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars,  methods="RC", X=c("X1"), M=c("M1", "M2"), doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:8))
})


test_that("ameras binomial LINEXP   methods=RC X=NULL M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars,  methods="RC", X=NULL, M=c("M1", "M2"), doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:7))
})






test_that("ameras binomial LINEXP   methods=RC X=c(X1,X2) M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars,   methods="RC", X=c("X1","X2"), M="M1", doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:7))
})


test_that("ameras binomial LINEXP   methods=RC X=X1 M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars,   methods="RC", X=c("X1"), M="M1", doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:6))
})


test_that("ameras binomial LINEXP   methods=RC X=NULL M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars,  methods="RC", X=NULL, M="M1", doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:5))
})







test_that("ameras binomial LINEXP   methods=RC X=c(X1,X2) M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars,   methods="RC", X=c("X1","X2"), M=NULL, doseRRmod="LINEXP"))
})


test_that("ameras binomial LINEXP   methods=RC X=X1 M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars,  methods="RC", X=c("X1"), M=NULL, doseRRmod="LINEXP"))
})

test_that("ameras binomial LINEXP   methods=RC X=NULL M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars,   methods="RC", X=NULL, M=NULL, doseRRmod="LINEXP"))
})


# ------------------------

test_that("ameras binomial  EXP deg=1 methods=RC X=c(X1,X2) M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,  methods="RC", X=c("X1","X2"), M=c("M1", "M2"), doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=c(X1,X2) M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,  methods="RC", X=c("X1","X2"), M=c("M1", "M2"), doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=X1 M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1, methods="RC", X=c("X1"), M=c("M1", "M2"), doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=X1 M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,  methods="RC", X=c("X1"), M=c("M1", "M2"), doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=NULL M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,  methods="RC", X=NULL, M=c("M1", "M2"), doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=NULL M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2, methods="RC", X=NULL, M=c("M1", "M2"), doseRRmod="EXP"))
})






test_that("ameras binomial  EXP deg=1 methods=RC X=c(X1,X2) M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,  methods="RC", X=c("X1","X2"), M="M1", doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=c(X1,X2) M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2, methods="RC", X=c("X1","X2"), M="M1", doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=X1 M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,methods="RC", X=c("X1"), M="M1", doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=X1 M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,  methods="RC", X=c("X1"), M="M1", doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=NULL M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,  methods="RC", X=NULL, M="M1", doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=NULL M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,  methods="RC", X=NULL, M="M1", doseRRmod="EXP"))
})






test_that("ameras binomial  EXP deg=1 methods=RC X=c(X1,X2) M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,  methods="RC", X=c("X1","X2"), M=NULL, doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=c(X1,X2) M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,  methods="RC", X=c("X1","X2"), M=NULL, doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=X1 M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,  methods="RC", X=c("X1"), M=NULL, doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=X1 M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,  methods="RC", X=c("X1"), M=NULL, doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=NULL M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,  methods="RC", X=NULL, M=NULL, doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=NULL M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,  methods="RC", X=NULL, M=NULL, doseRRmod="EXP"))
})








