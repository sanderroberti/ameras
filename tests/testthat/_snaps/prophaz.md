# prophaz snapshot: RC

    Code
      fit[[method]]$coefficients
    Output
              dose dose_squared 
        0.57628615  -0.03399134 

---

    Code
      fit[[method]]$sd
    Output
              dose dose_squared 
        0.08714064   0.01737043 

---

    Code
      fit[[method]]$CI
    Output
                         lower        upper
      dose          0.40549364 7.470787e-01
      dose_squared -0.06803676 5.408559e-05

# prophaz snapshot: ERC

    Code
      fit[[method]]$coefficients
    Output
              dose dose_squared 
        0.32943869  -0.01111156 

---

    Code
      fit[[method]]$sd
    Output
              dose dose_squared 
                NA           NA 

---

    Code
      fit[[method]]$CI
    Output
                   lower upper
      dose            NA    NA
      dose_squared    NA    NA

# prophaz snapshot: MCML

    Code
      fit[[method]]$coefficients
    Output
              dose dose_squared 
         0.5813187   -0.0388958 

---

    Code
      fit[[method]]$sd
    Output
              dose dose_squared 
        0.07838692   0.01445502 

---

    Code
      fit[[method]]$CI
    Output
                         lower       upper
      dose          0.42768318  0.73495428
      dose_squared -0.06722712 -0.01056447

# prophaz snapshot: FMA

    Code
      fit[[method]]$coefficients
    Output
              dose dose_squared 
        0.58211959  -0.03905762 

---

    Code
      fit[[method]]$sd
    Output
              dose dose_squared 
        0.07932717   0.01471867 

---

    Code
      fit[[method]]$CI
    Output
                         lower       upper
      dose          0.42720904  0.73862428
      dose_squared -0.06829587 -0.01031357

# prophaz snapshot: BMA

    Code
      fit[[method]]$coefficients
    Output
              dose dose_squared        h0[1]        h0[2] 
        0.52057980  -0.02854849   0.46348289   0.49840474 
             h0[3]        h0[4]        h0[5]        h0[6] 
        0.37022931   0.41722507   0.44918176   0.57739139 
             h0[7]        h0[8]        h0[9]       h0[10] 
        0.38676450   0.44689106   0.38199630   0.43709017 

---

    Code
      fit[[method]]$sd
    Output
              dose dose_squared        h0[1]        h0[2] 
        0.09486925   0.01617072   0.06362265   0.07447182 
             h0[3]        h0[4]        h0[5]        h0[6] 
        0.05334474   0.05670693   0.06636509   0.08782715 
             h0[7]        h0[8]        h0[9]       h0[10] 
        0.05733367   0.06250574   0.05895969   0.06997510 

---

    Code
      fit[[method]]$CI
    Output
                         lower      upper
      dose          0.26925382 0.64209968
      dose_squared -0.04923217 0.01321467
      h0[1]         0.35268365 0.61753746
      h0[2]         0.37880925 0.67067594
      h0[3]         0.27318599 0.47640599
      h0[4]         0.32782194 0.54358534
      h0[5]         0.32608430 0.58135200
      h0[6]         0.44831170 0.76969770
      h0[7]         0.28529825 0.49805173
      h0[8]         0.33354654 0.59250315
      h0[9]         0.29187284 0.51234607
      h0[10]        0.30774140 0.60848039

