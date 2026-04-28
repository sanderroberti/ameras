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

    structure(list(BMA = c(0.520579799976228, -0.0285484884309152, 
    0.463482885760436, 0.498404740623583, 0.370229311923104, 0.417225072243013, 
    0.449181764399242, 0.577391390177076, 0.386764497528575, 0.446891063673647, 
    0.38199629571387, 0.437090166634038)), class = "data.frame", row.names = c("dose", 
    "dose_squared", "h0[1]", "h0[2]", "h0[3]", "h0[4]", "h0[5]", 
    "h0[6]", "h0[7]", "h0[8]", "h0[9]", "h0[10]"))

---

    c(dose = 0.0948692519103796, dose_squared = 0.0161707181298364, 
    "h0[1]" = 0.063622653260542, "h0[2]" = 0.0744718238518843, "h0[3]" = 0.0533447443201762, 
    "h0[4]" = 0.0567069263630293, "h0[5]" = 0.0663650919635696, "h0[6]" = 0.0878271516650539, 
    "h0[7]" = 0.0573336708539561, "h0[8]" = 0.0625057364697876, "h0[9]" = 0.0589596913443982, 
    "h0[10]" = 0.0699751018505027)

---

    structure(list(lower = c(0.269253823025207, -0.0492321679561386, 
    0.352683645849793, 0.378809246161618, 0.273185987131722, 0.327821940165829, 
    0.326084298655844, 0.448311698191613, 0.285298248905047, 0.333546541495982, 
    0.291872836414828, 0.307741397326048), upper = c(0.642099680525839, 
    0.0132146682850106, 0.61753745597517, 0.670675938182796, 0.476405988430706, 
    0.543585338991664, 0.581351996735346, 0.769697697935526, 0.498051734039191, 
    0.592503154462393, 0.51234606744854, 0.608480393401336)), class = "data.frame", row.names = c("dose", 
    "dose_squared", "h0[1]", "h0[2]", "h0[3]", "h0[4]", "h0[5]", 
    "h0[6]", "h0[7]", "h0[8]", "h0[9]", "h0[10]"))

---

    Code
      print(summary(fit)$summary_table, digits = 2)
    Output
         Method         Term Estimate    SE CI.lower CI.upper Rhat n.eff
      1     BMA         dose    0.521 0.095    0.269    0.642  1.5     7
      2     BMA dose_squared   -0.029 0.016   -0.049    0.013  1.6    12
      3     BMA        h0[1]    0.463 0.064    0.353    0.618  1.1    43
      4     BMA        h0[2]    0.498 0.074    0.379    0.671  1.0    55
      5     BMA        h0[3]    0.370 0.053    0.273    0.476  1.2    47
      6     BMA        h0[4]    0.417 0.057    0.328    0.544  1.1    42
      7     BMA        h0[5]    0.449 0.066    0.326    0.581  1.0    64
      8     BMA        h0[6]    0.577 0.088    0.448    0.770  1.1    41
      9     BMA        h0[7]    0.387 0.057    0.285    0.498  1.3    63
      10    BMA        h0[8]    0.447 0.063    0.334    0.593  1.0   107
      11    BMA        h0[9]    0.382 0.059    0.292    0.512  1.3    35
      12    BMA       h0[10]    0.437 0.070    0.308    0.608  1.1    71

