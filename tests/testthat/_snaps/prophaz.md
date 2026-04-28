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

    structure(list(BMA = c(0.593733608470395, -0.0414268433298425, 
    0.435292940645859, 0.468312631020603, 0.361145286966187, 0.398655056485894, 
    0.423494435311201, 0.530014491158026, 0.363376146103184, 0.416304878824973, 
    0.36500182771853, 0.408801899562158)), class = "data.frame", row.names = c("dose", 
    "dose_squared", "h0[1]", "h0[2]", "h0[3]", "h0[4]", "h0[5]", 
    "h0[6]", "h0[7]", "h0[8]", "h0[9]", "h0[10]"))

---

    c(dose = 0.0858714637186771, dose_squared = 0.0159958734208164, 
    "h0[1]" = 0.0631247900751881, "h0[2]" = 0.069704174169237, "h0[3]" = 0.0534302447360716, 
    "h0[4]" = 0.0578064805636316, "h0[5]" = 0.0621258115697789, "h0[6]" = 0.0809826783480894, 
    "h0[7]" = 0.0527224498942387, "h0[8]" = 0.0619635623464885, "h0[9]" = 0.0524023190936904, 
    "h0[10]" = 0.0596156582935786)

---

    structure(list(lower = c(0.426056112543774, -0.0705717755449735, 
    0.326168770410745, 0.341929480718954, 0.263644597840207, 0.300288351640787, 
    0.306923599676445, 0.388890305589162, 0.266535997514222, 0.311903611451052, 
    0.268238268689863, 0.306578383936537), upper = c(0.739280585296266, 
    -0.00924914504354091, 0.568679601913492, 0.614810335470171, 0.473086060067528, 
    0.519107642884882, 0.555222316677207, 0.697255352685562, 0.471068595526541, 
    0.548109448491731, 0.474296446186198, 0.544839307789043)), class = "data.frame", row.names = c("dose", 
    "dose_squared", "h0[1]", "h0[2]", "h0[3]", "h0[4]", "h0[5]", 
    "h0[6]", "h0[7]", "h0[8]", "h0[9]", "h0[10]"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
         Method         Term Estimate    SE CI.lower CI.upper Rhat n.eff
      1     BMA         dose    0.594 0.086    0.426   0.7393    1    72
      2     BMA dose_squared   -0.041 0.016   -0.071  -0.0092    1    73
      3     BMA        h0[1]    0.435 0.063    0.326   0.5687    1   283
      4     BMA        h0[2]    0.468 0.070    0.342   0.6148    1   235
      5     BMA        h0[3]    0.361 0.053    0.264   0.4731    1   254
      6     BMA        h0[4]    0.399 0.058    0.300   0.5191    1   546
      7     BMA        h0[5]    0.423 0.062    0.307   0.5552    1   362
      8     BMA        h0[6]    0.530 0.081    0.389   0.6973    1   270
      9     BMA        h0[7]    0.363 0.053    0.267   0.4711    1   306
      10    BMA        h0[8]    0.416 0.062    0.312   0.5481    1   328
      11    BMA        h0[9]    0.365 0.052    0.268   0.4743    1   365
      12    BMA       h0[10]    0.409 0.060    0.307   0.5448    1   356

