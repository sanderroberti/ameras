# proflik/wald.transformed snapshot: RC

    structure(list(RC = c(-0.695663582187652, 0.0358348096332653, 
    0.273524431388668)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0853170991392384, dose = 0.204795009277685, 
    dose_squared = 0.077517578544823)

---

    structure(list(lower = c(0, 0.124824485619196), upper = c(0.481936149545829, 
    0.388528138266399), pval.lower = c(0.86103023251831, 0.049785273193913
    ), pval.upper = c(0.0500507010332296, 0.0493297310517574), iter.lower = c(NA, 
    7L), iter.upper = c(11L, 7L)), row.names = c("dose", "dose_squared"
    ), class = "data.frame")

---

    Code
      print(summary(fit1)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper pval.lower pval.upper
      1     RC  (Intercept) -0.69566 0.08532       NA       NA         NA         NA
      2     RC         dose  0.03583 0.20480   0.0000   0.4819    0.86103    0.05005
      3     RC dose_squared  0.27352 0.07752   0.1248   0.3885    0.04979    0.04933

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
      print(summary(fit2)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE   CI.lower  CI.upper
      1     RC  (Intercept) -0.69566 0.08532 -8.629e-01   -0.5284
      2     RC         dose  0.03583 0.20480  4.895e-07 2623.6137
      3     RC dose_squared  0.27352 0.07752  1.560e-01    0.4747

# proflik/wald.transformed snapshot: ERC

    structure(list(ERC = c(-0.693434923733658, 5.33950999146791e-08, 
    0.315000416519103)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0525498471312224, dose = 0.000228253766612278, 
    dose_squared = 0.0469447229554221)

---

    structure(list(lower = c(0, 0.23166893036439), upper = c(NA, 
    0.420152375346261), pval.lower = c(0.992315558085395, 0.0493096605911812
    ), pval.upper = c(NA, 0.0494652222864629), iter.lower = c(NA, 
    6L), iter.upper = c(NA, 4L)), row.names = c("dose", "dose_squared"
    ), class = "data.frame")

---

    Code
      print(summary(fit1)$summary_table, digits = 4)
    Output
        Method         Term   Estimate        SE CI.lower CI.upper pval.lower
      1    ERC  (Intercept) -6.934e-01 0.0525498       NA       NA         NA
      2    ERC         dose  5.340e-08 0.0002283   0.0000       NA    0.99232
      3    ERC dose_squared  3.150e-01 0.0469447   0.2317   0.4202    0.04931
        pval.upper
      1         NA
      2         NA
      3    0.04947

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
      print(summary(fit2)$summary_table, digits = 4)
    Output
        Method         Term   Estimate        SE CI.lower CI.upper
      1    ERC  (Intercept) -6.934e-01 0.0525498  -0.7964  -0.5904
      2    ERC         dose  5.340e-08 0.0002283   0.0000      Inf
      3    ERC dose_squared  3.150e-01 0.0469447   0.2349   0.4214

# proflik/wald.transformed snapshot: MCML

    structure(list(MCML = c(-0.678139870284149, 0.0083653968680493, 
    0.275598584484292)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0841326277237695, dose = 0.197716781910587, 
    dose_squared = 0.0748177402132151)

---

    structure(list(lower = c(0, 0.131427164024257), upper = c(0.441921125360596, 
    0.38098269776637), pval.lower = c(0.966228359151512, 0.0498346194362645
    ), pval.upper = c(0.0500763939116691, 0.0494530130996915), iter.lower = c(NA, 
    7L), iter.upper = c(14L, 7L)), row.names = c("dose", "dose_squared"
    ), class = "data.frame")

---

    Code
      print(summary(fit1)$summary_table, digits = 4)
    Output
        Method         Term  Estimate      SE CI.lower CI.upper pval.lower pval.upper
      1   MCML  (Intercept) -0.678140 0.08413       NA       NA         NA         NA
      2   MCML         dose  0.008365 0.19772   0.0000   0.4419    0.96623    0.05008
      3   MCML dose_squared  0.275599 0.07482   0.1314   0.3810    0.04983    0.04945

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
      print(summary(fit2)$summary_table, digits = 4)
    Output
        Method         Term  Estimate      SE   CI.lower   CI.upper
      1   MCML  (Intercept) -0.678140 0.08413 -8.430e-01 -5.132e-01
      2   MCML         dose  0.008365 0.19772  6.372e-23  1.098e+18
      3   MCML dose_squared  0.275599 0.07482  1.610e-01  4.674e-01

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
      print(summary(fit1)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1    FMA  (Intercept) -0.67942 0.08349  -0.8447  -0.5160
      2    FMA         dose  0.01296 0.19530  -0.3645   0.4016
      3    FMA dose_squared  0.27406 0.07436   0.1266   0.4193

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
      print(summary(fit2)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1    FMA  (Intercept) -0.67942 0.08349  -0.8478  -0.5195
      2    FMA         dose  0.01296 0.19530  -0.3663   0.3993
      3    FMA dose_squared  0.27406 0.07436   0.1238   0.4162

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
      print(summary(fit1)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper Rhat n.eff
      1    BMA  (Intercept)  -0.7234 0.10837  -0.9958  -0.5342 1.26    18
      2    BMA         dose   0.1343 0.24208  -0.2428   0.6851 1.40    13
      3    BMA dose_squared   0.2481 0.06946   0.1286   0.3999 1.15    58

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
      print(summary(fit2)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper Rhat n.eff
      1    BMA  (Intercept)  -0.7234 0.10837  -0.9129  -0.5102 1.26    18
      2    BMA         dose   0.1343 0.24208  -0.2420   0.6938 1.40    13
      3    BMA dose_squared   0.2481 0.06946   0.1055   0.3744 1.15    58

