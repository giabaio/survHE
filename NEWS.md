# 2.0.53 (development version)

## November 2026
* `make.surv()`: when `newdata` is given, the design matrix `des.mat` now keeps the profiles in the order in which they are listed in `newdata` (PR [#74](https://github.com/giabaio/survHE/pull/74)). Previously, the rows of `des.mat` followed the order in which the profiles first appear in the data, and repeated profiles were collapsed to a single row. For MLE models the survival curves `S` are computed directly from `newdata`; thus, when the two orders differed (e.g. with only categorical covariates and `newdata` listing the profiles in a different order from the data), the curves were attached to the wrong profiles. This affected the labels in `plot()`, `psa.plot()`, `summary()` and `make.transition.probs()`. For INLA and HMC models, the curves and `des.mat` were consistent with each other but in data order; they now follow the order of `newdata` too (thanks to [@chenyangsu](https://github.com/chenyangsu) for this).
* `make.surv()` (MLE): for the Royston-Parmar model (`distr="rps"`), the survival curves are now computed using the `scale` and `timescale` of the fitted model (PR [#75](https://github.com/giabaio/survHE/pull/75)), in the same way as was already done for the `knots`. Previously, `make.surv()` always called `flexsurv::psurvspline()` with the defaults `scale="hazard"` and `timescale="log"`; thus, models fitted with `scale="odds"`, `scale="normal"` or `timescale="identity"` produced wrong survival curves (e.g. for the `bc` data with `scale="odds"`, $S(5)$ for the group "Poor" was 0.04 instead of the 0.23 given by `flexsurv`). This affected `plot()`, `psa.plot()`, `summary()`, `make.transition.probs()` and `write.surv()`. The problem with `timescale` only became visible after version 2.0.52, which fixed `runMLE()` so that a user-specified `timescale` is actually used in the fit (thanks to [@chenyangsu](https://github.com/chenyangsu) for this).
* `make.surv()` (HMC): for the Royston-Parmar model, the survival curves are now always computed with `scale="hazard"` and `timescale="log"`, which are the only settings implemented in `survHEhmc`. Previously, values of `scale` or `timescale` passed to `make.surv()` through `...` were forwarded to `flexsurv::psurvspline()`; thus, the curves did not match the fitted model.

# 2.0.52 

## November 2026
* `fit.models()` (MLE): user-defined initial values can now be passed to `flexsurv` through the option `inits` (PRs [#65](https://github.com/giabaio/survHE/pull/65) and [#66](https://github.com/giabaio/survHE/pull/66)). For a single distribution, `inits` can be a numeric vector or a list with one numeric vector. For several distributions, it must be a list with one element per distribution (in the same order as `distr`); each element is a numeric vector or `NULL` (to use the `flexsurv` default for that model). Initial values are on the natural scale of the parameters, in the order used by `flexsurv`; covariate effects that are not given are initialised at 0 (thanks to [@yansi5](https://github.com/yansi5) for this).
* `inits` also works for the Royston-Parmar model (`distr="rps"`), where it is passed to `flexsurv::flexsurvspline`.
* Fixes a bug in the hidden function `runMLE()`: for `distr="rps"`, the option `timescale` was read from `scale`, so a user-specified `timescale` was ignored.
* Updates the documentation of `fit.models()`: the list of available distributions now matches `load_availables()` (adding `"rps"` for MLE, `"weibullPH"` for INLA and `"rps"` and `"polyweibull"` for HMC) and the MLE-specific options (`inits`, `weights`, `subset` and the `"rps"` options `k`, `knots`, `bknots`, `scale` and `timescale`) are now described.
* `make.surv()` (MLE): all covariate profiles now use the same draws of the model parameters, so the simulations are paired across profiles, as needed for PSA (PR [#68](https://github.com/giabaio/survHE/pull/68)). Previously, when `newdata` was given, each profile had its own call to `normboot.flexsurvreg()`, so the parameter draws were independent across profiles. In addition, the sub-sampling of `nsim` draws picked different rows for each profile, which broke the pairing even when `newdata` was not given. Now there is a single call to `normboot.flexsurvreg()` with all the profiles and one set of sampled rows, shared by all profiles (thanks to [@bitowaqr](https://github.com/bitowaqr) for this).
* `make_state_occupancy()`: the computation of the state occupancy now restarts from `start` at the beginning of each simulation (PR [#70](https://github.com/giabaio/survHE/pull/70)). Previously it only restarted when `profile` changed; thus, with a single profile and `nsim>1`, each simulation carried on from the last row of the previous one, rather than from `start`. With several profiles the bug was hidden, because the profile changes at each simulation boundary  (thanks to [@bitowaqr](https://github.com/bitowaqr) for this).
* `make_data_multi_state()`: deaths without progression (transition 2, Pre-progression -> Death) are now identified by `death==1 & prog==0` (PR [#72](https://github.com/giabaio/survHE/pull/72)). Previously the condition was `death==1 & prog_t==0`, which never happens in typical data (e.g. for TA174 it gives 0 direct deaths instead of 48); thus, the Pre-progression -> Death hazard was badly wrong. The documentation of `prog_t` and `death_t` now explains that they contain the censoring time (or, for `prog_t`, the death time) when the event has not happened  (thanks to [@bitowaqr](https://github.com/bitowaqr) for this).

## April 2026
* Adds support for the handling and plotting of Gamma models fitted in INLA (which will soon be available in `survHEinla` too).
Specifically, adds the `gam` support in the hidden function `load_availables()`. It also adds the function `rescale_stats_inla_gam()`,
which takes care of the rescaling for the `print` method.

## March 2026
* Fixes a small typo in the utility function `survHE:::make_surv_curve_plot` -- replaces `size` with `linewidth` in a call to `geom_line`. This avoids a warning.
* Adds facility to implement the newly available Gamma survival model from INLA.
* Fix slight imprecision -- INLA can't do RPS as a built-in model, so `survHE:::load_availables()` has been modified accordingly.

# 2.0.51 
## January 2026

* Aligns with changes in `dplyr` and adds `id` to the global functions.

## July 2025

* Slightly modifies the `print` method, so that the `original=TRUE` option can also be compounded with an additional option `print_priors=TRUE`. For HMC models, this also adds to the table of results a summary of the prior modelling assumptions.

# 2.0.5 (June 2025) - submitted to CRAN

* Adds a `ggplot` theme `theme_survHE` to be used with combined plots (eg for a `survextrap` element)
* Fix the tests, due to changes in the upcoming new release of `ggplot2`, which is now transitioning to S7 class (see [here](https://github.com/tidyverse/ggplot2/issues/6498))

# 2.0.4 (May 2025)

* Fixes a slight mistake in the table for the Gen F. Now works properly on the beta coefficients
* Updates DESCRIPTION with details of installation for the add-ons `survHEinla` and `survHEhmc`

# 2.0.3 (March 2025)

* Fixes dependency to R >= 4.1.0 because of the use of natural pipes
* Fixes tests --- needs to use `select(!contains("(Intercept)"))` instead of `select(-matches("Intercept"),everything())` in several parts of the code. This was previously used to determine the names of the strata, used to plot the survival curves

# 2.0.2 (October 2024)

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

* Some refactoring (mainly thanks to [@n8thangreen](https://github.com/n8thangreen)) with tidying up of the underlying code. One major change is that now 'make.surv' outputs an object named 'time' instead of 't', which was ambigious

_November 2022_

# 2.0

* [@Philip-Cooney](https://github.com/Philip-Cooney) found a little mistake in how the print method works for hmc objects. The utility functions were accessing the data for the     first model of the list of possible models (instead of the specific 'mod' one). Also, RPS would produce bizarre results if no covariates included.
    \item [@Philip-Cooney](https://github.com/Philip-Cooney) also found that 'make_sim_hmc' would break in the case of a RPS model with no covariates, because the matrix of beta coefficients would be turned into a vector, essentially, but the code would try to still subset a column. This has now been fixed so it's OK to make simulations off an RPS with no intercept model and that propagates to plots too.
    
_April 2022_

* This is a *major* change. In this version, the package is restructured to only perform, in its basic version, MLE estimates using flexsurv as 
    inferential engine. All the backbone functionalities are unchanged and the user can also expand (to revert to the "full" survHE including Bayesian modelling), by simply also adding the new packages survHEinla and/or survHEhmc. These now only contain the INLA and rstan calls and functionalities.
    
_January 2022_

# 1.1.4

* Contribution by Andrew Jones to update compatibility with the newer version of stan. Changes on StanHeaders + stan models to avoid complaints by 
the compiler because of declared variables with the same name of a function that was being defined. None of these are directly "visible" to the final user, though...

* Changes to '.Rbuildignore' to allow 'rstantool' to automatically configuring on package install. Also improves compatibility across versions of
'rstan'. See PR [#42](https://github.com/giabaio/survHE/pull/42)

_September 2021_

# 1.1.3

* Adds an option 'what' to plot so that 'survival', 'hazard' and 'cumhazard' can be specified (and the plot is modified to the various different scales)

* Updates the Gamma, GenGamma and GenF models in HMC to include for the possibility that the data contain no censoring. Also fix a small typo in the print method for Gamma/HMC models.

* Updates in INLA means now the Gompertz model is also available for survival modelling. fit.models(), make.surv() and print() have now been updated so that the Gompertz model can be run under 'survHE'. In order to run the Gompertz model using INLA, the *testing* version (>=21.03.21) needs to be installed (see instructions [here](https://www.r-inla.org/download-install)).

* Related to this, to improve computational stability, *all* the INLA models are now run by 'survHE' using the following trick: first the times are rescaled (on the fly) in the interval [0-1] (by simply recomputing 'time=time/max(time'). The resulting models are *not* directly comparable to other inferential engines (because they are fitted to different data), but 'survHE' automatically rescales the estimates and model fitting statistics (eg *ICs) so that the 'plot' and 'print' methods give the correct answers. make.surv() has also been updated to reflect this.

* There are also changes to make.transition.probs(), which has been updated and streamlined to compute transition probabilities off the survival curves fitted in a 'survHE' object. The computation is quicker and now based on the more robust relationship between the cumulative hazard function
and the transition probabilities.

* A new function make_data_multi_state() to create a dataset in the format required to analyse data in a multi-state framework.

* A new function three_state_mm() added to fit a standard 3-states Markov model, based on survival curves that are then mapped onto transition probabilities (this function is under testing, though).

* Adds two new datasets (TA174 and msmdata), both from the MDM paper by Williams et al (2017) that can be used for analysis of multi-state data.

_June 2021_
