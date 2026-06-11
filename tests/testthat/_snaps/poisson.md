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
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate     SE CI.lower CI.upper
      1     RC  (Intercept)    -0.88 0.0408   -0.964   -0.805
      2     RC         dose     0.62 0.0398    0.542    0.699
      3     RC dose_squared    -0.04 0.0074   -0.055   -0.026

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
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper
      1    ERC  (Intercept)   -0.880 0.040   -0.959   -0.801
      2    ERC         dose    0.611 0.039    0.535    0.687
      3    ERC dose_squared   -0.038 0.007   -0.052   -0.024

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
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate     SE CI.lower CI.upper
      1   MCML  (Intercept)   -0.854 0.0385   -0.930   -0.779
      2   MCML         dose    0.594 0.0345    0.526    0.661
      3   MCML dose_squared   -0.037 0.0058   -0.048   -0.025

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
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate     SE CI.lower CI.upper
      1    FMA  (Intercept)   -0.858 0.0400   -0.938   -0.780
      2    FMA         dose    0.600 0.0391    0.528    0.685
      3    FMA dose_squared   -0.038 0.0073   -0.056   -0.026

# Poisson snapshot: BMA

    structure(list(BMA = c(-0.85217931252014, 0.596796037428263, 
    -0.0394542111389367)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    c("(Intercept)" = 0.114340566756727, dose = 0.150770058437058, 
    dose_squared = 0.0285392945914736)

---

    structure(list(lower = c(-1.08057982786042, 0.273738071459514, 
    -0.098169394133088), upper = c(-0.60058016404825, 0.892978104234397, 
    0.0199593865839747)), class = "data.frame", row.names = c("(Intercept)", 
    "dose", "dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method         Term Estimate    SE CI.lower CI.upper Rhat n.eff
      1    BMA  (Intercept)   -0.852 0.114   -1.081    -0.60  1.0     7
      2    BMA         dose    0.597 0.151    0.274     0.89  1.2     3
      3    BMA dose_squared   -0.039 0.029   -0.098     0.02  1.5     3

