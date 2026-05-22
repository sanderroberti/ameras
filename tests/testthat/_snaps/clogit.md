# clogit snapshot: RC

    structure(list(RC = c(0.636933079842829, -0.0469518747182343)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    c(dose = 0.0994212957461848, dose_squared = 0.0220451378778829
    )

---

    structure(list(lower = c(0.442070920884002, -0.0901595509931046
    ), upper = c(0.831795238801657, -0.00374419844336402)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1     RC         dose    0.637 0.099     0.44   0.8318
      2     RC dose_squared   -0.047 0.022    -0.09  -0.0037

# clogit snapshot: ERC

    structure(list(ERC = c(0.211216746541785, 0.0529938610513979)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    c(dose = 0.0353015895183663, dose_squared = 0.00892501075843195
    )

---

    structure(list(lower = c(0.142026902488767, 0.0355011614032402
    ), upper = c(0.280406590594797, 0.0704865606995584)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate     SE CI.lower CI.upper
      1    ERC         dose    0.211 0.0353    0.142     0.28
      2    ERC dose_squared    0.053 0.0089    0.036     0.07

# clogit snapshot: MCML

    structure(list(MCML = c(0.668032212804173, -0.056784395590818
    )), class = "data.frame", row.names = c("dose", "dose_squared"
    ))

---

    c(dose = 0.0981315804650418, dose_squared = 0.0208508161639185
    )

---

    structure(list(lower = c(0.475697849346696, -0.0976512443203638
    ), upper = c(0.860366576261649, -0.0159175468612722)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1   MCML         dose    0.668 0.098    0.476    0.860
      2   MCML dose_squared   -0.057 0.021   -0.098   -0.016

# clogit snapshot: FMA

    structure(list(FMA = c(0.665505653115859, -0.0562187934345348
    )), class = "data.frame", row.names = c("dose", "dose_squared"
    ))

---

    c(dose = 0.102484025982308, dose_squared = 0.0221437635004394
    )

---

    structure(list(lower = c(0.460406659549313, -0.0992885679303401
    ), upper = c(0.863557680009724, -0.01171218044798)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1    FMA         dose    0.666 0.102    0.460    0.864
      2    FMA dose_squared   -0.056 0.022   -0.099   -0.012

# clogit snapshot: BMA

    structure(list(BMA = c(0.62957288938639, -0.0491256940956094)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    c(dose = 0.137109688051081, dose_squared = 0.0304239879548994
    )

---

    structure(list(lower = c(0.363967101239792, -0.0986573570228772
    ), upper = c(0.858500724121195, 0.00677148649021496)), class = "data.frame", row.names = c("dose", 
    "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate   SE CI.lower CI.upper Rhat n.eff
      1    BMA         dose    0.630 0.14    0.364   0.8585  3.3    49
      2    BMA dose_squared   -0.049 0.03   -0.099   0.0068  3.6    29

# clogit snapshot: FMA 1-par

    structure(list(FMA = 1.09177057785782), class = "data.frame", row.names = "dose")

---

    c(dose = 0.211208164015115)

---

    structure(list(lower = 0.679044763960655, upper = 1.50712257637798), class = "data.frame", row.names = "dose")

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method Term Estimate   SE CI.lower CI.upper
      1    FMA dose      1.1 0.21     0.68      1.5

