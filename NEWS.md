# survHE 1.1.2
This is a summary of the changes made to the development version 1.1.2 of `survHE`. The changes are updated to 12 October 2020 and are relative to the version published on CRAN (1.1.1).

- The default position of the legend in `plot` is moved a bit further up.    
- There is a fix to a minor issue with the `make.surv` function (which would spill on the `plot` method) for models fitted using HMC and with an intercept only (no covariates). With this minor fix, `make.surv` completes the PSA estimation for such models and in turns, `plot` shows the survival curves correctly.
