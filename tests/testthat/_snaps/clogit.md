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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower  CI.upper
      1     RC         dose  0.63693 0.09942  0.44207  0.831795
      2     RC dose_squared -0.04695 0.02205 -0.09016 -0.003744

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate       SE CI.lower CI.upper
      1    ERC         dose  0.21122 0.035302   0.1420  0.28041
      2    ERC dose_squared  0.05299 0.008925   0.0355  0.07049

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1   MCML         dose  0.66803 0.09813  0.47570  0.86037
      2   MCML dose_squared -0.05678 0.02085 -0.09765 -0.01592

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper
      1    FMA         dose  0.66551 0.10248  0.46041  0.86356
      2    FMA dose_squared -0.05622 0.02214 -0.09929 -0.01171

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
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower CI.upper Rhat n.eff
      1    BMA         dose  0.62957 0.13711  0.36397 0.858501 3.33    49
      2    BMA dose_squared -0.04913 0.03042 -0.09866 0.006771 3.57    29

