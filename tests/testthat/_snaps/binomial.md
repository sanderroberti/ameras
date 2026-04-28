# snapshot: binomial_ERR_deg1

    structure(list(RC = c(-1.12138445681014, 0.43747816756955, 0.829152319286417
    )), class = "data.frame", row.names = c("(Intercept)", "X1", 
    "dose"))

---

    c("(Intercept)" = 0.0861480826034211, X1 = 0.0761152077180876, 
    dose = 0.141835536319699)

---

    structure(list(lower = c(-1.29023159605003, 0.288295101766313, 
    0.551159776371883), upper = c(-0.952537317570258, 0.586661233372787, 
    1.10714486220095)), class = "data.frame", row.names = c("(Intercept)", 
    "X1", "dose"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method        Term Estimate    SE CI.lower CI.upper
      1     RC (Intercept)    -1.12 0.086    -1.29    -0.95
      2     RC          X1     0.44 0.076     0.29     0.59
      3     RC        dose     0.83 0.142     0.55     1.11

# snapshot: binomial_ERR_deg2

    structure(list(RC = c(-0.930459974364865, 0.442578152424314, 
    0.0355699195390312, 0.283828439726294)), class = "data.frame", row.names = c("(Intercept)", 
    "X1", "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0958670525012646, X1 = 0.0765527212331611, 
    dose = 0.208874156520848, dose_squared = 0.0794109809669211)

---

    structure(list(lower = c(-1.11835594457135, 0.292537575888783, 
    -0.373815904543013, 0.128185777054133), upper = c(-0.742564004158376, 
    0.592618728959844, 0.444955743621076, 0.439471102398455)), class = "data.frame", row.names = c("(Intercept)", 
    "X1", "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1     RC  (Intercept)   -0.930 0.096    -1.12    -0.74
      2     RC           X1    0.443 0.077     0.29     0.59
      3     RC         dose    0.036 0.209    -0.37     0.44
      4     RC dose_squared    0.284 0.079     0.13     0.44

# snapshot: binomial_EXP_deg1

    structure(list(RC = c(-1.03751752518051, 0.442783042151235, 0.441545076940068
    )), class = "data.frame", row.names = c("(Intercept)", "X1", 
    "dose"))

---

    c("(Intercept)" = 0.069938585800817, X1 = 0.0765042646998762, 
    dose = 0.0408085033922897)

---

    structure(list(lower = c(-1.17459463447978, 0.292837438675759, 
    0.3615618800282), upper = c(-0.900440415881248, 0.592728645626711, 
    0.521528273851936)), class = "data.frame", row.names = c("(Intercept)", 
    "X1", "dose"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method        Term Estimate    SE CI.lower CI.upper
      1     RC (Intercept)    -1.04 0.070    -1.17    -0.90
      2     RC          X1     0.44 0.077     0.29     0.59
      3     RC        dose     0.44 0.041     0.36     0.52

# snapshot: binomial_EXP_deg2

    structure(list(RC = c(-1.00098240262208, 0.442428476057805, 0.363241568514543, 
    0.0223337183123334)), class = "data.frame", row.names = c("(Intercept)", 
    "X1", "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0826504750506285, X1 = 0.0765013616470964, 
    dose = 0.103699331834104, dose_squared = 0.027501634282057)

---

    structure(list(lower = c(-1.16297435702644, 0.292488562461223, 
    0.159994612898832, -0.0315684943964904), upper = c(-0.838990448217726, 
    0.592368389654388, 0.566488524130253, 0.0762359310211572)), class = "data.frame", row.names = c("(Intercept)", 
    "X1", "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1     RC  (Intercept)   -1.001 0.083   -1.163   -0.839
      2     RC           X1    0.442 0.077    0.292    0.592
      3     RC         dose    0.363 0.104    0.160    0.566
      4     RC dose_squared    0.022 0.028   -0.032    0.076

# snapshot: binomial_LINEXP_deg2

    structure(list(RC = c(-0.989261627830125, 0.442426615995325, 
    0.310382666344636, 0.352282802770536)), class = "data.frame", row.names = c("(Intercept)", 
    "X1", "dose_linear", "dose_exponential"))

---

    c("(Intercept)" = 0.0839891703082157, X1 = 0.0765152415844903, 
    dose_linear = 0.115231089136967, dose_exponential = 0.109420939226206
    )

---

    structure(list(lower = c(-1.15387737672563, 0.292459498221342, 
    0.084533881736856, 0.137821702732626), upper = c(-0.824645878934621, 
    0.592393733769307, 0.536231450952416, 0.566743902808447)), class = "data.frame", row.names = c("(Intercept)", 
    "X1", "dose_linear", "dose_exponential"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method             Term Estimate    SE CI.lower CI.upper
      1     RC      (Intercept)    -0.99 0.084   -1.154    -0.82
      2     RC               X1     0.44 0.077    0.292     0.59
      3     RC      dose_linear     0.31 0.115    0.085     0.54
      4     RC dose_exponential     0.35 0.109    0.138     0.57

# binomial snapshot: ERC

    structure(list(ERC = c(-0.742108872660136, 0.293763533354128, 
    0.0445786596682705)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0705045165265905, dose = 0.103780566801031, 
    dose_squared = 0.027108601176109)

---

    structure(list(lower = c(-0.880295185799663, 0.0903573601289543, 
    -0.0085532223081632), upper = c(-0.60392255952061, 0.497169706579302, 
    0.0977105416447042)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1    ERC  (Intercept)   -0.742 0.071  -0.8803   -0.604
      2    ERC         dose    0.294 0.104   0.0904    0.497
      3    ERC dose_squared    0.045 0.027  -0.0086    0.098

# binomial snapshot: MCML

    structure(list(MCML = c(-0.732701069038212, 0.291000601776992, 
    0.0407221187083245)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0698236734295423, dose = 0.103397739211373, 
    dose_squared = 0.0276278889020142)

---

    structure(list(lower = c(-0.869552954228401, 0.0883447568398362, 
    -0.0134275485084972), upper = c(-0.595849183848023, 0.493656446714148, 
    0.0948717859251461)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1   MCML  (Intercept)   -0.733 0.070   -0.870   -0.596
      2   MCML         dose    0.291 0.103    0.088    0.494
      3   MCML dose_squared    0.041 0.028   -0.013    0.095

# binomial snapshot: FMA

    structure(list(FMA = c(-0.733806899724863, 0.293742722421584, 
    0.0398867147023287)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0703078399632643, dose = 0.104813195539034, 
    dose_squared = 0.0281101186064282)

---

    structure(list(lower = c(-0.871991768545499, 0.0900154912248243, 
    -0.0163552402022485), upper = c(-0.596009473907393, 0.501185640210265, 
    0.0946192633548769)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1    FMA  (Intercept)    -0.73 0.070   -0.872   -0.596
      2    FMA         dose     0.29 0.105    0.090    0.501
      3    FMA dose_squared     0.04 0.028   -0.016    0.095

# binomial snapshot: BMA

    structure(list(BMA = c(-0.730173444353113, 0.287133677575896, 
    0.0421037845656069)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.061525993306574, dose = 0.0901418668188222, 
    dose_squared = 0.0242758155720054)

---

    structure(list(lower = c(-0.857291851271742, 0.107033728876186, 
    0.00426550225354544), upper = c(-0.618610154718928, 0.439121194911745, 
    0.0886645543721738)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper Rhat n.eff
      1    BMA  (Intercept)   -0.730 0.062  -0.8573   -0.619  1.2    44
      2    BMA         dose    0.287 0.090   0.1070    0.439  2.1    23
      3    BMA dose_squared    0.042 0.024   0.0043    0.089  2.2    35

