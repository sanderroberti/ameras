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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1     RC  (Intercept)  -0.8596 0.03544  -0.9291  -0.7902
      2     RC         dose   0.4124 0.04825   0.3178   0.5070
      3     RC dose_squared   0.1924 0.01132   0.1702   0.2146
      4     RC        sigma   1.1007 0.01421   1.0729   1.1286

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1    ERC  (Intercept)  -0.8576 0.03542  -0.9270  -0.7882
      2    ERC         dose   0.4074 0.04825   0.3128   0.5019
      3    ERC dose_squared   0.1938 0.01135   0.1715   0.2160
      4    ERC        sigma   1.0999 0.01418   1.0721   1.1277

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1   MCML  (Intercept)  -0.8633 0.03671  -0.9352  -0.7913
      2   MCML         dose   0.4880 0.04980   0.3904   0.5856
      3   MCML dose_squared   0.1483 0.01130   0.1261   0.1704
      4   MCML        sigma   1.1510 0.01486   1.1219   1.1802

# gaussian snapshot: FMA

    structure(list(FMA = c(-0.863778905799158, 0.490749307581192, 
    0.146947464269109, 1.1511772699684)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    c("(Intercept)" = 0.0370273042047476, dose = 0.0519364859403298, 
    dose_squared = 0.0129118107700265, sigma = 0.0148852541053793
    )

---

    structure(list(lower = c(-0.936517795605244, 0.390260380881284, 
    0.118305280366457, 1.12194626692474), upper = c(-0.791068923763542, 
    0.594056623308337, 0.170739665607116, 1.18030645823)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1    FMA  (Intercept)  -0.8638 0.03703  -0.9365  -0.7911
      2    FMA         dose   0.4907 0.05194   0.3903   0.5941
      3    FMA dose_squared   0.1469 0.01291   0.1183   0.1707
      4    FMA        sigma   1.1512 0.01489   1.1219   1.1803

# gaussian snapshot: BMA

    structure(list(BMA = c(-0.86075388734466, 0.483858407793307, 
    0.147825002240133, 1.15099747564059)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    c("(Intercept)" = 0.0375232386622422, dose = 0.0452318592245169, 
    dose_squared = 0.0111250825137037, sigma = 0.0139880095194784
    )

---

    structure(list(lower = c(-0.937124709210879, 0.395445532983898, 
    0.122642690868872, 1.12289146189322), upper = c(-0.790254797080576, 
    0.566691532270053, 0.16606777228885, 1.17564663726968)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared", "sigma"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper Rhat n.eff
      1    BMA  (Intercept)  -0.8608 0.03752  -0.9371  -0.7903 1.00    81
      2    BMA         dose   0.4839 0.04523   0.3954   0.5667 1.03    73
      3    BMA dose_squared   0.1478 0.01113   0.1226   0.1661 1.12    73
      4    BMA        sigma   1.1510 0.01399   1.1229   1.1756 1.23    91

