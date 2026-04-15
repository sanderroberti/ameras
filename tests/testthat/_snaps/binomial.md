# snapshot: binomial_ERR_deg1

    Code
      fit$RC$coefficients
    Output
      (Intercept)          X1        dose 
       -1.1213845   0.4374782   0.8291523 

---

    Code
      fit$RC$sd
    Output
      (Intercept)          X1        dose 
       0.08614808  0.07611521  0.14183554 

---

    Code
      fit$RC$CI
    Output
                       lower      upper
      (Intercept) -1.2902347 -0.9525342
      X1           0.2882924  0.5866640
      dose         0.5511547  1.1071500

# snapshot: binomial_ERR_deg2

    Code
      fit$RC$coefficients
    Output
       (Intercept)           X1         dose dose_squared 
       -0.93045997   0.44257815   0.03556992   0.28382844 

---

    Code
      fit$RC$sd
    Output
       (Intercept)           X1         dose dose_squared 
        0.09586705   0.07655272   0.20887416   0.07941098 

---

    Code
      fit$RC$CI
    Output
                        lower      upper
      (Intercept)  -1.1183594 -0.7425606
      X1            0.2925348  0.5926215
      dose         -0.3738234  0.4449633
      dose_squared  0.1281829  0.4394740

# snapshot: binomial_EXP_deg1

    Code
      fit$RC$coefficients
    Output
      (Intercept)          X1        dose 
       -1.0375175   0.4427830   0.4415451 

---

    Code
      fit$RC$sd
    Output
      (Intercept)          X1        dose 
       0.06993859  0.07650426  0.04080850 

---

    Code
      fit$RC$CI
    Output
                       lower      upper
      (Intercept) -1.1745972 -0.9004379
      X1           0.2928347  0.5927314
      dose         0.3615604  0.5215297

# snapshot: binomial_EXP_deg2

    Code
      fit$RC$coefficients
    Output
       (Intercept)           X1         dose dose_squared 
       -1.00098240   0.44242848   0.36324157   0.02233372 

---

    Code
      fit$RC$sd
    Output
       (Intercept)           X1         dose dose_squared 
        0.08265048   0.07650136   0.10369933   0.02750163 

---

    Code
      fit$RC$CI
    Output
                         lower       upper
      (Intercept)  -1.16297733 -0.83898747
      X1            0.29248581  0.59237114
      dose          0.15999088  0.56649226
      dose_squared -0.03156948  0.07623692

# snapshot: binomial_LINEXP_deg2

    Code
      fit$RC$coefficients
    Output
           (Intercept)               X1      dose_linear dose_exponential 
            -0.9892616        0.4424266        0.3103827        0.3522828 

---

    Code
      fit$RC$sd
    Output
           (Intercept)               X1      dose_linear dose_exponential 
            0.08398917       0.07651524       0.11523109       0.10942094 

---

    Code
      fit$RC$CI
    Output
                             lower      upper
      (Intercept)      -1.15388040 -0.8246429
      X1                0.29245674  0.5923965
      dose_linear       0.08452973  0.5362356
      dose_exponential  0.13781776  0.5667478

# binomial snapshot: ERC

    Code
      fit[[method]]$coefficients
    Output
       (Intercept)         dose dose_squared 
       -0.74210887   0.29376353   0.04457866 

---

    Code
      fit[[method]]$sd
    Output
       (Intercept)         dose dose_squared 
        0.07050452   0.10378057   0.02710860 

---

    Code
      fit[[method]]$CI
    Output
                          lower       upper
      (Intercept)  -0.880297725 -0.60392002
      dose          0.090353622  0.49717344
      dose_squared -0.008554199  0.09771152

# binomial snapshot: MCML

    Code
      fit[[method]]$coefficients
    Output
       (Intercept)         dose dose_squared 
       -0.73270107   0.29100060   0.04072212 

---

    Code
      fit[[method]]$sd
    Output
       (Intercept)         dose dose_squared 
        0.06982367   0.10339774   0.02762789 

---

    Code
      fit[[method]]$CI
    Output
                         lower       upper
      (Intercept)  -0.86955547 -0.59584667
      dose          0.08834103  0.49366017
      dose_squared -0.01342854  0.09487278

# binomial snapshot: FMA

    Code
      fit[[method]]$coefficients
    Output
       (Intercept)         dose dose_squared 
       -0.73380690   0.29374272   0.03988671 

---

    Code
      fit[[method]]$sd
    Output
       (Intercept)         dose dose_squared 
        0.07030784   0.10481320   0.02811012 

---

    Code
      fit[[method]]$CI
    Output
                         lower       upper
      (Intercept)  -0.87199177 -0.59600947
      dose          0.09001549  0.50118564
      dose_squared -0.01635524  0.09461926

# binomial snapshot: BMA

    Code
      fit[[method]]$coefficients
    Output
       (Intercept)         dose dose_squared 
       -0.73017344   0.28713368   0.04210378 

---

    Code
      fit[[method]]$sd
    Output
       (Intercept)         dose dose_squared 
        0.06152599   0.09014187   0.02427582 

---

    Code
      fit[[method]]$CI
    Output
                          lower       upper
      (Intercept)  -0.857291851 -0.61861015
      dose          0.107033729  0.43912119
      dose_squared  0.004265502  0.08866455

