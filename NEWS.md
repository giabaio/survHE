# survHE 1.1.3

* Added a `NEWS.md` file to track changes to the package.

_September 2024_

Patch code in support of plotting (`R/utils_fit_models.R`)

* In `make_surv_curve_plot` replace 
```
geom_step(data = datakm, aes(x = time, y = S, group=as.factor(strata)),
                color="darkgrey") + 
      geom_ribbon(data = datakm,
                  aes(x = time, y = S, ymin=lower, ymax=upper, group=as.factor(strata)),
                  alpha = 0.2) 
```
with
```
geom_step(data = datakm, aes(x = time, y = S, group=as.factor(strata:object_name)),
                color="darkgrey") + 
      geom_ribbon(data = datakm,
                  aes(x = time, y = S, ymin=lower, ymax=upper, group=as.factor(strata:object_name)),
                  alpha = 0.2) 
```
This means that when plotting two or more `survHE` objects, the KM is added and displayed correctly

* Adds a utility function `make_newdata` that can be used to generate profiles of covariates, to then plot specific groups of individuals' survival or hazard curves.