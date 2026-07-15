# proflik/wald.transformed snapshot: RC

    structure(list(RC = c(-0.695663582187652, 0.0358348096332653,
    0.273524431388668)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0853170991392384, dose = 0.204795009277685,
    dose_squared = 0.077517578544823)

---

    structure(list(lower = c(NA, 0.124824485619156), upper = c(0.481936149545665,
    0.388528138266399), pval.lower = c(NA, 0.0497852731937511), pval.upper = c(0.0500507010336095,
    0.0493297310517574), iter.lower = c(NA, 7L), iter.upper = c(11L,
    7L)), row.names = c("dose", "dose_squared"), class = "data.frame")

---

    Code
      print(summary(fit1)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1     RC  (Intercept)   -0.696 0.085       NA       NA
      2     RC         dose    0.036 0.205       NA     0.48
      3     RC dose_squared    0.274 0.078     0.12     0.39

---

    structure(list(RC = c(-0.695663582187652, 0.0358348096332653,
    0.273524431388668)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0853170991392384, dose = 0.204795009277685,
    dose_squared = 0.077517578544823)

---

    structure(list(lower = c(-0.862882023765993, 4.89452239688687e-07,
    0.155981926298724), upper = c(-0.528445140609312, 2623.61365895298,
    0.474675836352198)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    Code
      print(summary(fit2)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1     RC  (Intercept)   -0.696 0.085 -8.6e-01    -0.53
      2     RC         dose    0.036 0.205  4.9e-07  2623.61
      3     RC dose_squared    0.274 0.078  1.6e-01     0.47

# proflik/wald.transformed snapshot: ERC

    structure(list(ERC = c(-0.693434923733658, 5.33950999146791e-08,
    0.315000416519103)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0525498471312224, dose = 0.000228253766612278,
    dose_squared = 0.0469447229554221)

---

    structure(list(lower = c(NA, 0.23166893036439), upper = c(NA,
    0.420152375346261), pval.lower = c(NA, 0.0493096605911812), pval.upper = c(NA,
    0.0494652222864629), iter.lower = c(NA, 6L), iter.upper = c(NA,
    4)), row.names = c("dose", "dose_squared"), class = "data.frame")

---

    Code
      print(summary(fit1)$summary_table, digits = 2)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1    ERC  (Intercept) -6.9e-01 0.05255       NA       NA
      2    ERC         dose  5.3e-08 0.00023       NA       NA
      3    ERC dose_squared  3.2e-01 0.04694     0.23     0.42

---

    structure(list(ERC = c(-0.693434923733658, 5.33950999146791e-08,
    0.315000416519103)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0525498471312224, dose = 0.000228253766612278,
    dose_squared = 0.0469447229554221)

---

    structure(list(lower = c(-0.79643073150394, 0, 0.234892430905717
    ), upper = c(-0.590439115963377, Inf, 0.421391646834024)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    Code
      print(summary(fit2)$summary_table, digits = 2)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1    ERC  (Intercept) -6.9e-01 0.05255    -0.80    -0.59
      2    ERC         dose  5.3e-08 0.00023     0.00      Inf
      3    ERC dose_squared  3.2e-01 0.04694     0.23     0.42

# proflik/wald.transformed snapshot: MCML

    structure(list(MCML = c(-0.678139870284149, 0.0083653968680493,
    0.275598584484292)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0841326277237695, dose = 0.197716781910587,
    dose_squared = 0.0748177402132151)

---

    structure(list(lower = c(NA, 0.13142716402424), upper = c(0.44192112536042,
    0.380982697766377), pval.lower = c(NA, 0.049834619436305), pval.upper = c(0.0500763939119137,
    0.0494530130995442), iter.lower = c(NA, 7L), iter.upper = c(14L,
    7L)), row.names = c("dose", "dose_squared"), class = "data.frame")

---

    Code
      print(summary(fit1)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1   MCML  (Intercept)  -0.6781 0.084       NA       NA
      2   MCML         dose   0.0084 0.198       NA     0.44
      3   MCML dose_squared   0.2756 0.075     0.13     0.38

---

    structure(list(MCML = c(-0.678139870284149, 0.0083653968680493,
    0.275598584484292)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0841326277237695, dose = 0.197716781910587,
    dose_squared = 0.0748177402132151)

---

    structure(list(lower = c(-0.843036790547453, 6.37197822923014e-23,
    0.160979662105735), upper = c(-0.513242950020845, 1098243939989482624,
    0.467385155968411)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    Code
      print(summary(fit2)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1   MCML  (Intercept)  -0.6781 0.084 -8.4e-01 -5.1e-01
      2   MCML         dose   0.0084 0.198  6.4e-23  1.1e+18
      3   MCML dose_squared   0.2756 0.075  1.6e-01  4.7e-01

# percentile/hpd snapshot: FMA

    structure(list(FMA = c(-0.679420383964918, 0.0129581869892215,
    0.274062556829787)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0834856681317757, dose = 0.195303870223553,
    dose_squared = 0.0743630356097565)

---

    structure(list(lower = c(-0.844712598991456, -0.364490456842211,
    0.126625280646566), upper = c(-0.515992835789905, 0.401607615396295,
    0.419327076449348)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    Code
      print(summary(fit1)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1    FMA  (Intercept)   -0.679 0.083    -0.84    -0.52
      2    FMA         dose    0.013 0.195    -0.36     0.40
      3    FMA dose_squared    0.274 0.074     0.13     0.42

---

    structure(list(FMA = c(-0.679420383964918, 0.0129581869892215,
    0.274062556829787)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0834856681317757, dose = 0.195303870223553,
    dose_squared = 0.0743630356097565)

---

    structure(list(lower = c(-0.847790185835778, -0.366287150520233,
    0.123789499311631), upper = c(-0.519530340551063, 0.399298136536058,
    0.416242200399084)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    Code
      print(summary(fit2)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1    FMA  (Intercept)   -0.679 0.083    -0.85    -0.52
      2    FMA         dose    0.013 0.195    -0.37     0.40
      3    FMA dose_squared    0.274 0.074     0.12     0.42

# percentile/hpd snapshot: BMA

    structure(list(BMA = c(-0.72335627769738, 0.134336359543895,
    0.248122664727861)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.108373977815855, dose = 0.242082441719764,
    dose_squared = 0.0694602084341149)

---

    structure(list(lower = c(-0.99584915733811, -0.242788192084572,
    0.128595431717399), upper = c(-0.53415204288286, 0.685070077998039,
    0.399869448439078)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    Code
      print(summary(fit1)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper Rhat n.eff
      1    BMA  (Intercept)    -0.72 0.108    -1.00    -0.53  1.3    18
      2    BMA         dose     0.13 0.242    -0.24     0.69  1.4    13
      3    BMA dose_squared     0.25 0.069     0.13     0.40  1.1    58

---

    structure(list(BMA = c(-0.72335627769738, 0.134336359543895,
    0.248122664727861)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.108373977815855, dose = 0.242082441719764,
    dose_squared = 0.0694602084341149)

---

    structure(list(lower = c(-0.912917097932078, -0.242048264480045,
    0.105492052659696), upper = c(-0.510214024881486, 0.693812077832168,
    0.374409980934656)), class = "data.frame", row.names = c("(Intercept)",
    "dose", "dose_squared"))

---

    Code
      print(summary(fit2)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper Rhat n.eff
      1    BMA  (Intercept)    -0.72 0.108    -0.91    -0.51  1.3    18
      2    BMA         dose     0.13 0.242    -0.24     0.69  1.4    13
      3    BMA dose_squared     0.25 0.069     0.11     0.37  1.1    58

