# gaussian snapshot: RC

    structure(list(RC = c(-0.859644829141104, 0.412396141643607, 
    0.192379759719246, 1.10074804303578)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    c("(Intercept)" = 0.0354442346838672, dose = 0.0482456043446934, 
    dose_squared = 0.011320561162419, sigma = 0.0142104511465155)

---

    structure(list(lower = c(-0.929114252581069, 0.317836494715639, 
    0.170191867556122, 1.07289607058454), upper = c(-0.790175405701139, 
    0.506955788571575, 0.214567651882371, 1.12860001548701)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1     RC  (Intercept)    -0.86 0.035    -0.93    -0.79
      2     RC         dose     0.41 0.048     0.32     0.51
      3     RC dose_squared     0.19 0.011     0.17     0.21
      4     RC        sigma     1.10 0.014     1.07     1.13

# gaussian snapshot: ERC

    structure(list(ERC = c(-0.857590768586508, 0.407372444065381, 
    0.19375103179669, 1.09986471318972)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    c("(Intercept)" = 0.0354221950080613, dose = 0.0482492720914936, 
    dose_squared = 0.0113501136639654, sigma = 0.0141843297477061
    )

---

    structure(list(lower = c(-0.927016995055663, 0.31280560848578, 
    0.171505217794882, 1.07206393773938), upper = c(-0.788164542117354, 
    0.501939279644982, 0.215996845798498, 1.12766548864006)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1    ERC  (Intercept)    -0.86 0.035    -0.93    -0.79
      2    ERC         dose     0.41 0.048     0.31     0.50
      3    ERC dose_squared     0.19 0.011     0.17     0.22
      4    ERC        sigma     1.10 0.014     1.07     1.13

# gaussian snapshot: MCML

    structure(list(MCML = c(-0.863286201325093, 0.48796776533509, 
    0.148289235834366, 1.15104527636048)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    c("(Intercept)" = 0.0367088203386334, dose = 0.0497968776967233, 
    dose_squared = 0.0113010207378352, sigma = 0.014859941178711)

---

    structure(list(lower = c(-0.935234167103766, 0.390367678506966, 
    0.126139642199668, 1.12192032683783), upper = c(-0.79133823554642, 
    0.585567852163213, 0.170438829469063, 1.18017022588314)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1   MCML  (Intercept)    -0.86 0.037    -0.94    -0.79
      2   MCML         dose     0.49 0.050     0.39     0.59
      3   MCML dose_squared     0.15 0.011     0.13     0.17
      4   MCML        sigma     1.15 0.015     1.12     1.18

# gaussian snapshot: FMA

    structure(list(FMA = c(-0.863778905799158, 0.490749307581192, 
    0.146947464269109, 1.1511772699684)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    c("(Intercept)" = 0.0370220947588556, dose = 0.0519113151824632, 
    dose_squared = 0.0129043157144241, sigma = 0.0148845339274785
    )

---

    structure(list(lower = c(-0.936517795605244, 0.390260380881284, 
    0.118305280366457, 1.12194626692474), upper = c(-0.791068923763542, 
    0.594056623308337, 0.170739665607116, 1.18030645823)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1    FMA  (Intercept)    -0.86 0.037    -0.94    -0.79
      2    FMA         dose     0.49 0.052     0.39     0.59
      3    FMA dose_squared     0.15 0.013     0.12     0.17
      4    FMA        sigma     1.15 0.015     1.12     1.18

# gaussian snapshot: BMA

    structure(list(BMA = c(-0.861039230169073, 0.487861545224281, 
    0.14629264262881, 1.15221514057702)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    c("(Intercept)" = 0.0418836048547209, dose = 0.0483052036046024, 
    dose_squared = 0.0134316186907068, sigma = 0.0152618089835416
    )

---

    structure(list(lower = c(-0.932181495784226, 0.399744498981121, 
    0.116521528523675, 1.12572933691506), upper = c(-0.772456861961563, 
    0.589715462876947, 0.168357634915395, 1.18362357773664)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper Rhat n.eff
      1    BMA  (Intercept)    -0.86 0.042    -0.93    -0.77 0.99   100
      2    BMA         dose     0.49 0.048     0.40     0.59 1.05    56
      3    BMA dose_squared     0.15 0.013     0.12     0.17 1.17    32
      4    BMA        sigma     1.15 0.015     1.13     1.18 1.00   210

