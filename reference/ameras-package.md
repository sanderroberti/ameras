# Analyze multiple exposure realizations in association studies

Analyze association studies with multiple realizations of a noisy or
uncertain exposure. These can be obtained from e.g. a two-dimensional
Monte Carlo dosimetry system (Simon et al 2015
\<doi:10.1667/RR13729.1\>) to characterize exposure uncertainty. Methods
include regression calibration (Carroll et al. 2006
[doi:10.1201/9781420010138](https://doi.org/10.1201/9781420010138) ),
extended regression calibration (Little et al. 2023
[doi:10.1038/s41598-023-42283-y](https://doi.org/10.1038/s41598-023-42283-y)
), Monte Carlo maximum likelihood (Stayner et al. 2007
[doi:10.1667/RR0677.1](https://doi.org/10.1667/RR0677.1) ), frequentist
model averaging (Kwon et al. 2023
[doi:10.1371/journal.pone.0290498](https://doi.org/10.1371/journal.pone.0290498)
), and Bayesian model averaging (Kwon et al. 2016
[doi:10.1002/sim.6635](https://doi.org/10.1002/sim.6635) ). Supported
model families are Gaussian, binomial, multinomial, Poisson,
proportional hazards, and conditional logistic.

## Details

The main function is
[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

## Author

Sander Roberti \<sander.roberti@nih.gov\>, William Wheeler
\<WheelerB@imsweb.com\>, Ruth Pfeiffer \<pfeiffer@mail.nih.gov\>, and
Deukwoo Kwon \<DKwon@uams.edu

## References

Roberti, S., Kwon D., Wheeler W., Pfeiffer R. (in preparation). ameras:
An R Package to Analyze Multiple Exposure Realizations in Association
Studies
