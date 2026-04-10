# Example data

Data includes outcomes of all six supported types in the appropriately
named columns. For proportional hazards regression, the observed exit
time is `time` and event status is `event`. For conditional logistic
regression, the matched set variable is `setnr`. The data has 10
exposure replicates in columns `V1`-`V10`.

## Examples

``` r
 data(data, package="ameras")

 # Display a few rows of the data
 data[1:5, ]
#>    Y.gaussian Y.binomial Y.poisson      time status setnr Y.clogit
#> 1 -0.32647093          0         0 0.3027656      0     1        0
#> 2 -0.18734036          1         0 0.1973514      1     1        1
#> 3  0.08404044          0         2 0.3027656      0     1        0
#> 4  0.22432504          0         0 0.2360258      1     1        0
#> 5 -0.46317255          0         0 0.3027656      0     2        0
#>   Y.multinomial X1 X2 M1 M2         V1         V2         V3         V4
#> 1             3  0  0  0  1 0.42868043 0.61542487 0.41960219 0.49265549
#> 2             2  1  0  1  0 0.73321154 0.35512449 0.41876478 0.49235658
#> 3             2  0  0  1  0 0.70369712 0.43407408 1.04115924 0.79882088
#> 4             3  0  0  1  0 0.01845324 0.01373367 0.02733303 0.01912686
#> 5             3  1  0  0  0 0.39389441 0.40087181 0.61932032 0.51715526
#>           V5         V6         V7         V8         V9        V10
#> 1 0.31363762 0.42218455 0.42464021 0.29630858 0.38211182 0.45751570
#> 2 0.49515815 0.56837639 0.61126842 0.67723449 0.53361810 0.49393510
#> 3 0.66613754 0.72346942 0.64077434 0.79894534 0.98278177 1.06068250
#> 4 0.01917956 0.03056413 0.01536966 0.02135999 0.01548655 0.01596626
#> 5 0.36440322 0.60255525 0.47512525 0.52567606 0.53391825 0.56026531
```
