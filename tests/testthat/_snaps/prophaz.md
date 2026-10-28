# prophaz snapshot: RC

    structure(list(RC = c(0.57628615469363, -0.0339913365375492)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    c(dose = 0.0871406378193708, dose_squared = 0.0173704325187007
    )

---

    structure(list(lower = c(0.405493642977815, -0.0680367586700859
    ), upper = c(0.747078666409446, 5.40855949875754e-05)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1     RC         dose    0.576 0.087    0.405  7.5e-01
      2     RC dose_squared   -0.034 0.017   -0.068  5.4e-05

# prophaz snapshot: ERC

    structure(list(ERC = c(0.329438692494966, -0.0111115564109217
    )), class = "data.frame", row.names = c("dose", "dose_squared"
    ))

---

    c(dose = NA_real_, dose_squared = NA_real_)

---

    structure(list(lower = c(NA_real_, NA_real_), upper = c(NA_real_, 
    NA_real_)), class = "data.frame", row.names = c("dose", "dose_squared"
    ))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate SE CI.lower CI.upper
      1    ERC         dose    0.329 NA       NA       NA
      2    ERC dose_squared   -0.011 NA       NA       NA

# prophaz snapshot: MCML

    structure(list(MCML = c(0.581318732120769, -0.038895795621568
    )), class = "data.frame", row.names = c("dose", "dose_squared"
    ))

---

    c(dose = 0.0783869247717831, dose_squared = 0.0144550228208374
    )

---

    structure(list(lower = c(0.427683182709223, -0.0672271197461139
    ), upper = c(0.734954281532314, -0.0105644714970222)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1   MCML         dose    0.581 0.078    0.428    0.735
      2   MCML dose_squared   -0.039 0.014   -0.067   -0.011

# prophaz snapshot: FMA

    structure(list(FMA = c(0.582119594670597, -0.0390576222751472
    )), class = "data.frame", row.names = c("dose", "dose_squared"
    ))

---

    c(dose = 0.0793271699716439, dose_squared = 0.0147186668798125
    )

---

    structure(list(lower = c(0.427209037843427, -0.0682958725268815
    ), upper = c(0.738624278475251, -0.0103135694046743)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1    FMA         dose    0.582 0.079    0.427     0.74
      2    FMA dose_squared   -0.039 0.015   -0.068    -0.01

# prophaz snapshot: BMA

    structure(list(BMA = c(0.59279789774428, -0.0414212934705637, 
    0.437203687857773, 0.469562756212489, 0.361723426783169, 0.399297505241577, 
    0.425333707523287, 0.534964491715369, 0.363159201191068, 0.416056927906623, 
    0.364887564898759, 0.408584041024697)), class = "data.frame", row.names = c("dose", 
    "dose_squared", "h0[1]", "h0[2]", "h0[3]", "h0[4]", "h0[5]", 
    "h0[6]", "h0[7]", "h0[8]", "h0[9]", "h0[10]"))

---

    c(dose = 0.0865673624987242, dose_squared = 0.0160840746369599, 
    "h0[1]" = 0.0650190466925804, "h0[2]" = 0.0700296236233704, "h0[3]" = 0.0534606790704257, 
    "h0[4]" = 0.0590182695444318, "h0[5]" = 0.0627350246705671, "h0[6]" = 0.0783349475651827, 
    "h0[7]" = 0.0527499563241878, "h0[8]" = 0.0623892073388782, "h0[9]" = 0.0527739019161336, 
    "h0[10]" = 0.0594264808463789)

---

    structure(list(lower = c(0.415858268325095, -0.0710975251943536, 
    0.321676888497718, 0.342359632103572, 0.266385042051409, 0.292949976290781, 
    0.307543358100789, 0.393562245215681, 0.266547266821683, 0.308795516624463, 
    0.272666810004812, 0.304932474667468), upper = c(0.742507227407166, 
    -0.0102085761388503, 0.57949726921734, 0.620076779930237, 0.475326275089895, 
    0.526948482468846, 0.559164404210773, 0.696252553582887, 0.476942448795476, 
    0.550814280324657, 0.479965602753711, 0.54198698857668)), class = "data.frame", row.names = c("dose", 
    "dose_squared", "h0[1]", "h0[2]", "h0[3]", "h0[4]", "h0[5]", 
    "h0[6]", "h0[7]", "h0[8]", "h0[9]", "h0[10]"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
         Method         Term Estimate    SE CI.lower CI.upper Rhat n.eff
      1     BMA         dose    0.593 0.087    0.416     0.74    1    71
      2     BMA dose_squared   -0.041 0.016   -0.071    -0.01    1    83
      3     BMA        h0[1]    0.437 0.065    0.322     0.58    1   622
      4     BMA        h0[2]    0.470 0.070    0.342     0.62    1   667
      5     BMA        h0[3]    0.362 0.053    0.266     0.48    1   622
      6     BMA        h0[4]    0.399 0.059    0.293     0.53    1   867
      7     BMA        h0[5]    0.425 0.063    0.308     0.56    1   720
      8     BMA        h0[6]    0.535 0.078    0.394     0.70    1   551
      9     BMA        h0[7]    0.363 0.053    0.267     0.48    1   773
      10    BMA        h0[8]    0.416 0.062    0.309     0.55    1   592
      11    BMA        h0[9]    0.365 0.053    0.273     0.48    1   480
      12    BMA       h0[10]    0.409 0.059    0.305     0.54    1   574

# prophaz snapshot: FMA 1-par

    structure(list(FMA = 1.06798681864581), class = "data.frame", row.names = "dose")

---

    c(dose = 0.188879612817662)

---

    structure(list(lower = 0.69966449201711, upper = 1.43913843180591), class = "data.frame", row.names = "dose")

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method Term Estimate   SE CI.lower CI.upper
      1    FMA dose      1.1 0.19      0.7      1.4

