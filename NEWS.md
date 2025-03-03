# 2.0.3

* Fixes dependency to R >= 4.1.0 because of the use of natural pipes
* Fixes tests --- needs to use `select(!contains("(Intercept)"))` instead of `select(-matches("Intercept"),everything())` in several parts of the code. This was previously used to determine the names of the strata, used to plot the survival curves

# 2.0.2

* Added a `NEWS.md` file to track changes to the package (to replace the old `.Rd` version).

_October 2024_

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

_September 2024_

# 2.0.1

* Some refactoring (mainly thanks to @n8thangreen) with tidying up of the underlying code. One major change is that now 'make.surv' outputs an object named 'time' instead of 't', which was ambigious

_November 2022_

# 2.0

* @Philip-Cooney found a little mistake in how the print method works for hmc objects. The utility functions were accessing the data for the     first model of the list of possible models (instead of the specific 'mod' one). Also, RPS would produce bizarre results if no covariates included.
    \item @Philip-Cooney also found that 'make_sim_hmc' would break in the case of a RPS model with no covariates, because the matrix of beta coefficients would be turned into a vector, essentially, but the code would try to still subset a column. This has now been fixed so it's OK to make simulations off an RPS with no intercept model and that propagates to plots too.
    
_April 2022_

* This is a *major* change. In this version, the package is restructured to only perform, in its basic version, MLE estimates using flexsurv as 
    inferential engine. All the backbone functionalities are unchanged and the user can also expand (to revert to the "full" survHE including Bayesian modelling), by simply also adding the new packages survHEinla and/or survHEhmc. These now only contain the INLA and rstan calls and functionalities.
    
_January 2022_

# 1.1.4

* Contribution by Andrew Jones to update compatibility with the newer version of stan. Changes on StanHeaders + stan models to avoid complaints by 
the compiler because of declared variables with the same name of a function that was being defined. None of these are directly "visible" to the final user, though...

* Changes to '.Rbuildignore' to allow 'rstantool' to automatically configuring on package install. Also improves compatibility across versions of
'rstan'. See https://github.com/giabaio/survHE/pull/42

_September 2021_

# 1.1.3

* Adds an option 'what' to plot so that 'survival', 'hazard' and 'cumhazard' can be specified (and the plot is modified to the various different scales)

* Updates the Gamma, GenGamma and GenF models in HMC to include for the possibility that the data contain no censoring. Also fix a small typo in the print method for Gamma/HMC models.

* Updates in INLA means now the Gompertz model is also available for survival modelling. fit.models(), make.surv() and print() have now been updated so that the Gompertz model can be run under 'survHE'. In order to run the Gompertz model using INLA, the *testing* version (>=21.03.21) needs to be installed (see instructions here: https://www.r-inla.org/download-install).

* Related to this, to improve computational stability, *all* the INLA models are now run by 'survHE' using the following trick: first the times are rescaled (on the fly) in the interval [0-1] (by simply recomputing 'time=time/max(time'). The resulting models are *not* directly comparable to other inferential engines (because they are fitted to different data), but 'survHE' automatically rescales the estimates and model fitting statistics (eg *ICs) so that the 'plot' and 'print' methods give the correct answers. make.surv() has also been updated to reflect this.

* There are also changes to make.transition.probs(), which has been updated and streamlined to compute transition probabilities off the survival curves fitted in a 'survHE' object. The computation is quicker and now based on the more robust relationship between the cumulative hazard function
and the transition probabilities.

* A new function make_data_multi_state() to create a dataset in the format required to analyse data in a multi-state framework.

* A new function three_state_mm() added to fit a standard 3-states Markov model, based on survival curves that are then mapped onto transition probabilities (this function is under testing, though).

* Adds two new datasets (TA174 and msmdata), both from the MDM paper by Williams et al (2017) that can be used for analysis of multi-state data.

_June 2021_
