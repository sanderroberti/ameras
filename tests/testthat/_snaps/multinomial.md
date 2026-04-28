# multinomial snapshot: RC

    structure(list(RC = c(-0.872532246181349, 0.544231834982583, 
    -0.0313499823285036, -0.109368671249109, 0.57635409412743, -0.041641603047693
    )), class = "data.frame", row.names = c("(1)_(Intercept)", "(1)_dose", 
    "(1)_dose_squared", "(2)_(Intercept)", "(2)_dose", "(2)_dose_squared"
    ))

---

    c("(1)_(Intercept)" = 0.0910621676999024, "(1)_dose" = 0.130761271583389, 
    "(1)_dose_squared" = 0.0332067423567577, "(2)_(Intercept)" = 0.0745044468806468, 
    "(2)_dose" = 0.113600566594363, "(2)_dose_squared" = 0.0301801567232621
    )

---

    structure(list(lower = c(-1.0510108152273, 0.287944452106479, 
    -0.0964340013916494, -0.255394703823254, 0.353701074979134, -0.100793623273061
    ), upper = c(-0.694053677135394, 0.800519217858686, 0.0337340367346422, 
    0.0366573613250359, 0.799007113275726, 0.0175104171776751)), class = "data.frame", row.names = c("(1)_(Intercept)", 
    "(1)_dose", "(1)_dose_squared", "(2)_(Intercept)", "(2)_dose", 
    "(2)_dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method             Term Estimate      SE CI.lower CI.upper
      1     RC  (1)_(Intercept) -0.87253 0.09106 -1.05101 -0.69405
      2     RC         (1)_dose  0.54423 0.13076  0.28794  0.80052
      3     RC (1)_dose_squared -0.03135 0.03321 -0.09643  0.03373
      4     RC  (2)_(Intercept) -0.10937 0.07450 -0.25539  0.03666
      5     RC         (2)_dose  0.57635 0.11360  0.35370  0.79901
      6     RC (2)_dose_squared -0.04164 0.03018 -0.10079  0.01751

# multinomial snapshot: ERC

    structure(list(ERC = c(-0.874400585060071, 0.557601599677255, 
    -0.0322519686441131, -0.110267556601327, 0.582400710001095, -0.0391875439518789
    )), class = "data.frame", row.names = c("(1)_(Intercept)", "(1)_dose", 
    "(1)_dose_squared", "(2)_(Intercept)", "(2)_dose", "(2)_dose_squared"
    ))

---

    c("(1)_(Intercept)" = 0.0891371445468678, "(1)_dose" = 0.123432405861906, 
    "(1)_dose_squared" = 0.029799250545539, "(2)_(Intercept)" = 0.0725942435845731, 
    "(2)_dose" = 0.105480906960947, "(2)_dose_squared" = 0.0263306217704128
    )

---

    structure(list(lower = c(-1.04910617805667, 0.315678529662789, 
    -0.090657426479655, -0.252549659512018, 0.375661931301018, -0.0907946143124343
    ), upper = c(-0.699694992063469, 0.79952466969172, 0.0261534891914288, 
    0.0320145463093646, 0.789139488701172, 0.0124195264086766)), class = "data.frame", row.names = c("(1)_(Intercept)", 
    "(1)_dose", "(1)_dose_squared", "(2)_(Intercept)", "(2)_dose", 
    "(2)_dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method             Term Estimate      SE CI.lower CI.upper
      1    ERC  (1)_(Intercept) -0.87440 0.08914 -1.04911 -0.69969
      2    ERC         (1)_dose  0.55760 0.12343  0.31568  0.79952
      3    ERC (1)_dose_squared -0.03225 0.02980 -0.09066  0.02615
      4    ERC  (2)_(Intercept) -0.11027 0.07259 -0.25255  0.03201
      5    ERC         (2)_dose  0.58240 0.10548  0.37566  0.78914
      6    ERC (2)_dose_squared -0.03919 0.02633 -0.09079  0.01242

# multinomial snapshot: MCML

    structure(list(MCML = c(-0.858804858198922, 0.52789110098132, 
    -0.0284676930064826, -0.104061905992065, 0.580767475654764, -0.0443267873734486
    )), class = "data.frame", row.names = c("(1)_(Intercept)", "(1)_dose", 
    "(1)_dose_squared", "(2)_(Intercept)", "(2)_dose", "(2)_dose_squared"
    ))

---

    c("(1)_(Intercept)" = 0.0890760787896232, "(1)_dose" = 0.124621303806615, 
    "(1)_dose_squared" = 0.0301293702372263, "(2)_(Intercept)" = 0.0727887409089817, 
    "(2)_dose" = 0.10711493856028, "(2)_dose_squared" = 0.0266682668222654
    )

---

    structure(list(lower = c(-1.03339076451064, 0.283637833813931, 
    -0.0875201735483191, -0.246725216653687, 0.370826053870395, -0.0965956298751931
    ), upper = c(-0.684218951887208, 0.772144368148709, 0.030584787535354, 
    0.0386014046695559, 0.790708897439134, 0.00794205512829594)), class = "data.frame", row.names = c("(1)_(Intercept)", 
    "(1)_dose", "(1)_dose_squared", "(2)_(Intercept)", "(2)_dose", 
    "(2)_dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method             Term Estimate      SE CI.lower  CI.upper
      1   MCML  (1)_(Intercept) -0.85880 0.08908 -1.03339 -0.684219
      2   MCML         (1)_dose  0.52789 0.12462  0.28364  0.772144
      3   MCML (1)_dose_squared -0.02847 0.03013 -0.08752  0.030585
      4   MCML  (2)_(Intercept) -0.10406 0.07279 -0.24673  0.038601
      5   MCML         (2)_dose  0.58077 0.10711  0.37083  0.790709
      6   MCML (2)_dose_squared -0.04433 0.02667 -0.09660  0.007942

# multinomial snapshot: FMA

    structure(list(FMA = c(-0.860734562162769, 0.532383338057357, 
    -0.0299074150772933, -0.103626336376565, 0.580004407582068, -0.0442824724335984
    )), class = "data.frame", row.names = c("(1)_(Intercept)", "(1)_dose", 
    "(1)_dose_squared", "(2)_(Intercept)", "(2)_dose", "(2)_dose_squared"
    ))

---

    c("(1)_(Intercept)" = 0.0895735088513275, "(1)_dose" = 0.127222774539065, 
    "(1)_dose_squared" = 0.0311233594341049, "(2)_(Intercept)" = 0.0729724815016855, 
    "(2)_dose" = 0.107933226214456, "(2)_dose_squared" = 0.0270365031956162
    )

---

    structure(list(lower = c(-1.03623344740624, 0.2832882096991, 
    -0.0910528321303053, -0.247042472166993, 0.369769439674991, -0.0980862262304044
    ), upper = c(-0.686453470307389, 0.782562545619539, 0.031407527453333, 
    0.0395264618081422, 0.793582522875995, 0.00826657847150827)), class = "data.frame", row.names = c("(1)_(Intercept)", 
    "(1)_dose", "(1)_dose_squared", "(2)_(Intercept)", "(2)_dose", 
    "(2)_dose_squared"))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method             Term Estimate      SE CI.lower  CI.upper
      1    FMA  (1)_(Intercept) -0.86073 0.08957 -1.03623 -0.686453
      2    FMA         (1)_dose  0.53238 0.12722  0.28329  0.782563
      3    FMA (1)_dose_squared -0.02991 0.03112 -0.09105  0.031408
      4    FMA  (2)_(Intercept) -0.10363 0.07297 -0.24704  0.039526
      5    FMA         (2)_dose  0.58000 0.10793  0.36977  0.793583
      6    FMA (2)_dose_squared -0.04428 0.02704 -0.09809  0.008267

# multinomial snapshot: BMA

    structure(list(BMA = c(-0.587718745014249, -0.122836686034902, 
    0.177646762797876, 0.137220659635992, -0.00983993483463549, 0.147174985137446
    )), class = "data.frame", row.names = c("(1)_(Intercept)", "(1)_dose", 
    "(1)_dose_squared", "(2)_(Intercept)", "(2)_dose", "(2)_dose_squared"
    ))

---

    c("(1)_(Intercept)" = 0.16884744034517, "(1)_dose" = 0.427307016088735, 
    "(1)_dose_squared" = 0.155092265686904, "(2)_(Intercept)" = 0.153138980744979, 
    "(2)_dose" = 0.39067109573887, "(2)_dose_squared" = 0.146645948536132
    )

---

    structure(list(lower = c(-0.899314610991756, -0.841229407356729, 
    -0.0233626602065779, -0.131990787086232, -0.713951712373224, 
    -0.0518910784512704), upper = c(-0.267715870383468, 0.476122000646149, 
    0.484626926664673, 0.422959754764158, 0.579695497596276, 0.43280446623031
    )), class = "data.frame", row.names = c("(1)_(Intercept)", "(1)_dose", 
    "(1)_dose_squared", "(2)_(Intercept)", "(2)_dose", "(2)_dose_squared"
    ))

---

    Code
      print(summary(fit)$summary_table, digits = 4)
    Output
        Method             Term Estimate     SE CI.lower CI.upper Rhat n.eff
      1    BMA  (1)_(Intercept) -0.58772 0.1688 -0.89931  -0.2677 1.97     7
      2    BMA         (1)_dose -0.12284 0.4273 -0.84123   0.4761 2.70     3
      3    BMA (1)_dose_squared  0.17765 0.1551 -0.02336   0.4846 2.72     3
      4    BMA  (2)_(Intercept)  0.13722 0.1531 -0.13199   0.4230 1.88     5
      5    BMA         (2)_dose -0.00984 0.3907 -0.71395   0.5797 2.43     2
      6    BMA (2)_dose_squared  0.14717 0.1466 -0.05189   0.4328 2.43     2

