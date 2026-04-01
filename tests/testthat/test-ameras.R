
# Test RC with all 45 combinations of doseRRmod, deg and lengths of X and M

test_that("ameras binomial ERR deg=1 methods=RC X=c(X1,X2) M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"), dosevars=dosevars, deg=1, methods="RC", X=c("X1","X2"), M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:6))
})

test_that("ameras binomial ERR deg=2 methods=RC X=c(X1,X2) M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=2, methods="RC", X=c("X1","X2"), M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:9))
})

test_that("ameras binomial ERR deg=1 methods=RC X=X1 M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=1, methods="RC", X=c("X1"), M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:5))
})

test_that("ameras binomial ERR deg=2 methods=RC X=X1 M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=2, methods="RC", X=c("X1"), M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:8))
})

test_that("ameras binomial ERR deg=1 methods=RC X=NULL M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=1, methods="RC", X=NULL, M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:4))
})

test_that("ameras binomial ERR deg=2 methods=RC X=NULL M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=2, methods="RC", X=NULL, M=c("M1", "M2"), doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:7))
})






test_that("ameras binomial ERR deg=1 methods=RC X=c(X1,X2) M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=1, methods="RC", X=c("X1","X2"), M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:5))
})

test_that("ameras binomial ERR deg=2 methods=RC X=c(X1,X2) M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial",CI=c("wald.orig"), dosevars=dosevars, deg=2, methods="RC", X=c("X1","X2"), M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:7))
})

test_that("ameras binomial ERR deg=1 methods=RC X=X1 M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=1, methods="RC", X=c("X1"), M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:4))
})

test_that("ameras binomial ERR deg=2 methods=RC X=X1 M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial",CI=c("wald.orig"), dosevars=dosevars, deg=2, methods="RC", X=c("X1"), M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:6))
})

test_that("ameras binomial ERR deg=1 methods=RC X=NULL M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",CI=c("wald.orig"), dosevars=dosevars, deg=1, methods="RC", X=NULL, M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:3))
})

test_that("ameras binomial ERR deg=2 methods=RC X=NULL M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=2, methods="RC", X=NULL, M="M1", doseRRmod="ERR", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:5))
})






test_that("ameras binomial ERR deg=1 methods=RC X=c(X1,X2) M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",CI=c("wald.orig"), dosevars=dosevars, deg=1, methods="RC", X=c("X1","X2"), M=NULL, doseRRmod="ERR"))
})

test_that("ameras binomial ERR deg=2 methods=RC X=c(X1,X2) M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=2, methods="RC", X=c("X1","X2"), M=NULL, doseRRmod="ERR"))
})

test_that("ameras binomial ERR deg=1 methods=RC X=X1 M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=1, methods="RC", X=c("X1"), M=NULL, doseRRmod="ERR"))
})

test_that("ameras binomial ERR deg=2 methods=RC X=X1 M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=2, methods="RC", X=c("X1"), M=NULL, doseRRmod="ERR"))
})

test_that("ameras binomial ERR deg=1 methods=RC X=NULL M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial",CI=c("wald.orig"), dosevars=dosevars, deg=1, methods="RC", X=NULL, M=NULL, doseRRmod="ERR"))
})

test_that("ameras binomial ERR deg=2 methods=RC X=NULL M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", CI=c("wald.orig"),dosevars=dosevars, deg=2, methods="RC", X=NULL, M=NULL, doseRRmod="ERR"))
})





#------------------------


test_that("ameras binomial LINEXP   methods=RC X=c(X1,X2) M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, CI=c("wald.orig"),  methods="RC", X=c("X1","X2"), M=c("M1", "M2"), doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:9))
})

test_that("ameras binomial LINEXP   methods=RC X=X1 M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, CI=c("wald.orig"),  methods="RC", X=c("X1"), M=c("M1", "M2"), doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:8))
})


test_that("ameras binomial LINEXP   methods=RC X=NULL M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, CI=c("wald.orig"),  methods="RC", X=NULL, M=c("M1", "M2"), doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:7))
})






test_that("ameras binomial LINEXP   methods=RC X=c(X1,X2) M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, CI=c("wald.orig"),  methods="RC", X=c("X1","X2"), M="M1", doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=4:7))
})


test_that("ameras binomial LINEXP   methods=RC X=X1 M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, CI=c("wald.orig"),  methods="RC", X=c("X1"), M="M1", doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=3:6))
})


test_that("ameras binomial LINEXP   methods=RC X=NULL M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, CI=c("wald.orig"),  methods="RC", X=NULL, M="M1", doseRRmod="LINEXP", transform=transform1, transform.jacobian=transform1.jacobian, index.t=2:5))
})







test_that("ameras binomial LINEXP   methods=RC X=c(X1,X2) M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, CI=c("wald.orig"),  methods="RC", X=c("X1","X2"), M=NULL, doseRRmod="LINEXP"))
})


test_that("ameras binomial LINEXP   methods=RC X=X1 M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, CI=c("wald.orig"),  methods="RC", X=c("X1"), M=NULL, doseRRmod="LINEXP"))
})

test_that("ameras binomial LINEXP   methods=RC X=NULL M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, CI=c("wald.orig"),  methods="RC", X=NULL, M=NULL, doseRRmod="LINEXP"))
})


# ------------------------

test_that("ameras binomial  EXP deg=1 methods=RC X=c(X1,X2) M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1, CI=c("wald.orig"), methods="RC", X=c("X1","X2"), M=c("M1", "M2"), doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=c(X1,X2) M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2, CI=c("wald.orig"),methods="RC", X=c("X1","X2"), M=c("M1", "M2"), doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=X1 M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,CI=c("wald.orig"), methods="RC", X=c("X1"), M=c("M1", "M2"), doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=X1 M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2, CI=c("wald.orig"),methods="RC", X=c("X1"), M=c("M1", "M2"), doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=NULL M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1, CI=c("wald.orig"),methods="RC", X=NULL, M=c("M1", "M2"), doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=NULL M=c(M1,M2)", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,CI=c("wald.orig"), methods="RC", X=NULL, M=c("M1", "M2"), doseRRmod="EXP"))
})






test_that("ameras binomial  EXP deg=1 methods=RC X=c(X1,X2) M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1, CI=c("wald.orig"),methods="RC", X=c("X1","X2"), M="M1", doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=c(X1,X2) M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,CI=c("wald.orig"), methods="RC", X=c("X1","X2"), M="M1", doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=X1 M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,CI=c("wald.orig"), methods="RC", X=c("X1"), M="M1", doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=X1 M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2, CI=c("wald.orig"),methods="RC", X=c("X1"), M="M1", doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=NULL M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,CI=c("wald.orig"), methods="RC", X=NULL, M="M1", doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=NULL M=M1", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,CI=c("wald.orig"), methods="RC", X=NULL, M="M1", doseRRmod="EXP"))
})






test_that("ameras binomial  EXP deg=1 methods=RC X=c(X1,X2) M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,CI=c("wald.orig"), methods="RC", X=c("X1","X2"), M=NULL, doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=c(X1,X2) M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,CI=c("wald.orig"), methods="RC", X=c("X1","X2"), M=NULL, doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=X1 M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1, CI=c("wald.orig"),methods="RC", X=c("X1"), M=NULL, doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=X1 M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,CI=c("wald.orig"), methods="RC", X=c("X1"), M=NULL, doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=1 methods=RC X=NULL M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=1,CI=c("wald.orig"), methods="RC", X=NULL, M=NULL, doseRRmod="EXP"))
})

test_that("ameras binomial  EXP deg=2 methods=RC X=NULL M=NULL", {
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="binomial", Y="Y.binomial", dosevars=dosevars, deg=2,CI=c("wald.orig"), methods="RC", X=NULL, M=NULL, doseRRmod="EXP"))
})




#-------------------------

# Basic test of all methods for all other outcome models

# Poisson, basic model, all methods 
test_that("ameras poisson EXP deg=2 methods=RC,ERC,MCML,FMA X=NULL M=NULL", {
  skip_on_cran()
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="poisson", Y="Y.poisson", dosevars=dosevars, deg=2, CI=c("wald.orig", "percentile"),methods=c("RC", "ERC", "MCML", "FMA", "BMA"), X=NULL, M=NULL, doseRRmod="EXP"))
})

# Gaussian, basic model, all methods 
test_that("ameras gaussian deg=2 methods=RC,ERC,MCML,FMA X=NULL M=NULL", {
  skip_on_cran()
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="gaussian", Y="Y.gaussian", dosevars=dosevars, deg=2, CI=c("wald.orig", "percentile"),methods=c("RC", "ERC", "MCML", "FMA", "BMA"), X=NULL, M=NULL))
})


# Multinomial, basic model, all methods except BMA
test_that("ameras multinomial EXP deg=2 methods=RC,ERC,MCML,FMA X=NULL M=NULL", {
  skip_on_cran()
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="multinomial", Y="Y.multinomial", dosevars=dosevars, deg=2,CI=c("wald.orig", "percentile"), methods=c("RC", "ERC", "MCML", "FMA", "BMA"), X=NULL, M=NULL, doseRRmod="EXP"))
})


# Prophaz, basic model, all methods 
test_that("ameras prophaz EXP deg=2 methods=RC,ERC,MCML,FMA X=NULL M=NULL", {
  skip_on_cran()
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_warning(ameras(data=data, family="prophaz", Y="status", exit="time", dosevars=dosevars, deg=2,CI=c("wald.orig", "percentile"), methods=c("RC", "ERC", "MCML", "FMA", "BMA"), X=NULL, M=NULL, doseRRmod="EXP"))
})


# Clogit, basic model, all methods
test_that("ameras clogit EXP deg=2 methods=RC,ERC,MCML,FMA X=NULL M=NULL", {
  skip_on_cran()
  data("data", package="ameras")
  dosevars <- paste0("V",1:10)
  expect_no_error(ameras(data=data, family="clogit", Y="Y.clogit", setnr="setnr", dosevars=dosevars, deg=2, CI=c("wald.orig", "percentile"),methods=c("RC", "ERC", "MCML", "FMA", "BMA"), X=NULL, M=NULL, doseRRmod="EXP"))
})