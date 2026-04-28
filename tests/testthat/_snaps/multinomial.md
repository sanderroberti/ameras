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
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method             Term Estimate    SE CI.lower CI.upper
      1     RC  (1)_(Intercept)   -0.873 0.091   -1.051   -0.694
      2     RC         (1)_dose    0.544 0.131    0.288    0.801
      3     RC (1)_dose_squared   -0.031 0.033   -0.096    0.034
      4     RC  (2)_(Intercept)   -0.109 0.075   -0.255    0.037
      5     RC         (2)_dose    0.576 0.114    0.354    0.799
      6     RC (2)_dose_squared   -0.042 0.030   -0.101    0.018

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
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method             Term Estimate    SE CI.lower CI.upper
      1    ERC  (1)_(Intercept)   -0.874 0.089   -1.049   -0.700
      2    ERC         (1)_dose    0.558 0.123    0.316    0.800
      3    ERC (1)_dose_squared   -0.032 0.030   -0.091    0.026
      4    ERC  (2)_(Intercept)   -0.110 0.073   -0.253    0.032
      5    ERC         (2)_dose    0.582 0.105    0.376    0.789
      6    ERC (2)_dose_squared   -0.039 0.026   -0.091    0.012

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
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method             Term Estimate    SE CI.lower CI.upper
      1   MCML  (1)_(Intercept)   -0.859 0.089   -1.033  -0.6842
      2   MCML         (1)_dose    0.528 0.125    0.284   0.7721
      3   MCML (1)_dose_squared   -0.028 0.030   -0.088   0.0306
      4   MCML  (2)_(Intercept)   -0.104 0.073   -0.247   0.0386
      5   MCML         (2)_dose    0.581 0.107    0.371   0.7907
      6   MCML (2)_dose_squared   -0.044 0.027   -0.097   0.0079

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
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method             Term Estimate    SE CI.lower CI.upper
      1    FMA  (1)_(Intercept)   -0.861 0.090   -1.036  -0.6865
      2    FMA         (1)_dose    0.532 0.127    0.283   0.7826
      3    FMA (1)_dose_squared   -0.030 0.031   -0.091   0.0314
      4    FMA  (2)_(Intercept)   -0.104 0.073   -0.247   0.0395
      5    FMA         (2)_dose    0.580 0.108    0.370   0.7936
      6    FMA (2)_dose_squared   -0.044 0.027   -0.098   0.0083

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
      print(summary(fit)$summary_table, digits = 2)
    Output
        Method             Term Estimate   SE CI.lower CI.upper Rhat n.eff
      1    BMA  (1)_(Intercept)  -0.5877 0.17   -0.899    -0.27  2.0     7
      2    BMA         (1)_dose  -0.1228 0.43   -0.841     0.48  2.7     3
      3    BMA (1)_dose_squared   0.1776 0.16   -0.023     0.48  2.7     3
      4    BMA  (2)_(Intercept)   0.1372 0.15   -0.132     0.42  1.9     5
      5    BMA         (2)_dose  -0.0098 0.39   -0.714     0.58  2.4     2
      6    BMA (2)_dose_squared   0.1472 0.15   -0.052     0.43  2.4     2

