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

    structure(list(BMA = c(0.533696101638339, -0.0314646605284873, 
    0.468050459824929, 0.497144184808074, 0.378122531576888, 0.422077059362406, 
    0.446882198491802, 0.563255604067662, 0.380876072980271, 0.441334868176305, 
    0.383620100011815, 0.430494554551972)), class = "data.frame", row.names = c("dose", 
    "dose_squared", "h0[1]", "h0[2]", "h0[3]", "h0[4]", "h0[5]", 
    "h0[6]", "h0[7]", "h0[8]", "h0[9]", "h0[10]"))

---

    c(dose = 0.115174360092467, dose_squared = 0.0207808664027294, 
    "h0[1]" = 0.0753052758072568, "h0[2]" = 0.0796729523810991, "h0[3]" = 0.0627608244251918, 
    "h0[4]" = 0.068585155809537, "h0[5]" = 0.068208183828608, "h0[6]" = 0.0902695323480563, 
    "h0[7]" = 0.0618813135478922, "h0[8]" = 0.0692661665000544, "h0[9]" = 0.0623140202440979, 
    "h0[10]" = 0.0715894859433258)

---

    structure(list(lower = c(0.206868962883771, -0.0613299720039562, 
    0.348353206297311, 0.364447740645051, 0.282459753623268, 0.314571773733772, 
    0.332133845497682, 0.398292175891734, 0.273382219140136, 0.328663879301948, 
    0.282174042338372, 0.306617893531167), upper = c(0.68951980566149, 
    0.0277577014007621, 0.64272193782141, 0.676860851554093, 0.543121285321583, 
    0.57058605730801, 0.615354672640777, 0.770140142768364, 0.512860427619233, 
    0.595768836840876, 0.525075869402902, 0.595138591536962)), class = "data.frame", row.names = c("dose", 
    "dose_squared", "h0[1]", "h0[2]", "h0[3]", "h0[4]", "h0[5]", 
    "h0[6]", "h0[7]", "h0[8]", "h0[9]", "h0[10]"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
         Method         Term Estimate    SE CI.lower CI.upper Rhat n.eff
      1     BMA         dose    0.534 0.115    0.207    0.690  1.3    21
      2     BMA dose_squared   -0.031 0.021   -0.061    0.028  1.2    25
      3     BMA        h0[1]    0.468 0.075    0.348    0.643  1.1    80
      4     BMA        h0[2]    0.497 0.080    0.364    0.677  1.2    63
      5     BMA        h0[3]    0.378 0.063    0.282    0.543  1.1    74
      6     BMA        h0[4]    0.422 0.069    0.315    0.571  1.0   123
      7     BMA        h0[5]    0.447 0.068    0.332    0.615  1.0    61
      8     BMA        h0[6]    0.563 0.090    0.398    0.770  1.2   122
      9     BMA        h0[7]    0.381 0.062    0.273    0.513  1.1   115
      10    BMA        h0[8]    0.441 0.069    0.329    0.596  1.1   114
      11    BMA        h0[9]    0.384 0.062    0.282    0.525  1.0   128
      12    BMA       h0[10]    0.430 0.072    0.307    0.595  1.0    73

