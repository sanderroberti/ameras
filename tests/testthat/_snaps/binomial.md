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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method        Term Estimate      SE CI.lower CI.upper
      1     RC (Intercept)  -1.1214 0.08615  -1.2902  -0.9525
      2     RC          X1   0.4375 0.07612   0.2883   0.5867
      3     RC        dose   0.8292 0.14184   0.5512   1.1071

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1     RC  (Intercept) -0.93046 0.09587  -1.1184  -0.7426
      2     RC           X1  0.44258 0.07655   0.2925   0.5926
      3     RC         dose  0.03557 0.20887  -0.3738   0.4450
      4     RC dose_squared  0.28383 0.07941   0.1282   0.4395

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method        Term Estimate      SE CI.lower CI.upper
      1     RC (Intercept)  -1.0375 0.06994  -1.1746  -0.9004
      2     RC          X1   0.4428 0.07650   0.2928   0.5927
      3     RC        dose   0.4415 0.04081   0.3616   0.5215

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1     RC  (Intercept) -1.00098 0.08265 -1.16297 -0.83899
      2     RC           X1  0.44243 0.07650  0.29249  0.59237
      3     RC         dose  0.36324 0.10370  0.15999  0.56649
      4     RC dose_squared  0.02233 0.02750 -0.03157  0.07624

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method             Term Estimate      SE CI.lower CI.upper
      1     RC      (Intercept)  -0.9893 0.08399 -1.15388  -0.8246
      2     RC               X1   0.4424 0.07652  0.29246   0.5924
      3     RC      dose_linear   0.3104 0.11523  0.08453   0.5362
      4     RC dose_exponential   0.3523 0.10942  0.13782   0.5667

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE  CI.lower CI.upper
      1    ERC  (Intercept) -0.74211 0.07050 -0.880295 -0.60392
      2    ERC         dose  0.29376 0.10378  0.090357  0.49717
      3    ERC dose_squared  0.04458 0.02711 -0.008553  0.09771

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1   MCML  (Intercept) -0.73270 0.06982 -0.86955 -0.59585
      2   MCML         dose  0.29100 0.10340  0.08834  0.49366
      3   MCML dose_squared  0.04072 0.02763 -0.01343  0.09487

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1    FMA  (Intercept) -0.73381 0.07031 -0.87199 -0.59601
      2    FMA         dose  0.29374 0.10481  0.09002  0.50119
      3    FMA dose_squared  0.03989 0.02811 -0.01636  0.09462

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE  CI.lower CI.upper Rhat n.eff
      1    BMA  (Intercept)  -0.7302 0.06153 -0.857292 -0.61861 1.25    44
      2    BMA         dose   0.2871 0.09014  0.107034  0.43912 2.15    23
      3    BMA dose_squared   0.0421 0.02428  0.004266  0.08866 2.22    35

