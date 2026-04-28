# Poisson snapshot: RC

    structure(list(RC = c(-0.884511877326043, 0.62049067304047, -0.0400924621455127
    )), class = "data.frame", row.names = c("(Intercept)", "dose", 
    "dose_squared"))

---

    c("(Intercept)" = 0.0407525835288347, dose = 0.0398084600441328, 
    dose_squared = 0.0073719270436653)

---

    structure(list(lower = c(-0.96438547331952, 0.542467525073967, 
    -0.0545411736477535), upper = c(-0.804638281332567, 0.698513821006972, 
    -0.0256437506432719)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate       SE CI.lower CI.upper
      1     RC  (Intercept) -0.88451 0.040753 -0.96439 -0.80464
      2     RC         dose  0.62049 0.039808  0.54247  0.69851
      3     RC dose_squared -0.04009 0.007372 -0.05454 -0.02564

# Poisson snapshot: ERC

    structure(list(ERC = c(-0.879900561003282, 0.611249720235731, 
    -0.0378744548428055)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0404559831820766, dose = 0.0387383055934341, 
    dose_squared = 0.00701626766896617)

---

    structure(list(lower = c(-0.95919283099931, 0.535324036450494, 
    -0.051626086779872), upper = c(-0.800608291007254, 0.687175404020969, 
    -0.024122822905739)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate       SE CI.lower CI.upper
      1    ERC  (Intercept) -0.87990 0.040456 -0.95919 -0.80061
      2    ERC         dose  0.61125 0.038738  0.53532  0.68718
      3    ERC dose_squared -0.03787 0.007016 -0.05163 -0.02412

# Poisson snapshot: MCML

    structure(list(MCML = c(-0.85434971404576, 0.593745669672536, 
    -0.0368254864743646)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0385394894778423, dose = 0.0345033491570035, 
    dose_squared = 0.00580625643204904)

---

    structure(list(lower = c(-0.929885725404891, 0.526120347978798, 
    -0.0482055399661847), upper = c(-0.778813702686629, 0.661370991366273, 
    -0.0254454329825444)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate       SE CI.lower CI.upper
      1   MCML  (Intercept) -0.85435 0.038539 -0.92989 -0.77881
      2   MCML         dose  0.59375 0.034503  0.52612  0.66137
      3   MCML dose_squared -0.03683 0.005806 -0.04821 -0.02545

# Poisson snapshot: FMA

    structure(list(FMA = c(-0.85774451193125, 0.599610478857139, 
    -0.0382776866437891)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0400221044672458, dose = 0.0391061771073553, 
    dose_squared = 0.0072962645111157)

---

    structure(list(lower = c(-0.937639840406445, 0.527981485004092, 
    -0.0562889521857883), upper = c(-0.780338412722336, 0.685022704940847, 
    -0.0259222692708027)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate       SE CI.lower CI.upper
      1    FMA  (Intercept) -0.85774 0.040022 -0.93764 -0.78034
      2    FMA         dose  0.59961 0.039106  0.52798  0.68502
      3    FMA dose_squared -0.03828 0.007296 -0.05629 -0.02592

# Poisson snapshot: BMA

    structure(list(BMA = c(-0.816042010970857, 0.548296349764407, 
    -0.030468230010967)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.0848227674909227, dose = 0.102316694541126, 
    dose_squared = 0.0174258321558769)

---

    structure(list(lower = c(-0.936287928910797, 0.327186090023886, 
    -0.0575862008232684), upper = c(-0.628417092320246, 0.703246734058374, 
    0.00519617749309753)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method         Term Estimate      SE CI.lower  CI.upper Rhat n.eff
      1    BMA  (Intercept) -0.81604 0.08482 -0.93629 -0.628417 1.39     7
      2    BMA         dose  0.54830 0.10232  0.32719  0.703247 1.77     4
      3    BMA dose_squared -0.03047 0.01743 -0.05759  0.005196 1.82     4

