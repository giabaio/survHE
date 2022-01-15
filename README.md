# survHEhmc [![Travis-CI Build Status](https://travis-ci.org/giabaio/survHE.svg?branch=hmc)](https://travis-ci.org/giabaio/survHE)[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/giabaio/survHE?branch=hmc&svg=true)](https://ci.appveyor.com/project/giabaio/survHE)
## Survival analysis in health economic evaluation using HMC/stan

This is a module to complement the package `survHE` and expand its functionalities to run survival analysis in health economic evaluation from a Bayesian perspective, using Hamiltonian Monte Carlo (via the R package [rstan](https://mc-stan.org/users/interfaces/rstan)). `survHEhmc` "depends" on the main installation of `survHE`. This means that you shouldn't use `survHEhmc` as a standalone package --- rather you use all the functions of `survHE` (to fit the models and the post-process the results); installing `survHEhmc` basically opens up a new option in the `survHE` function `fit.models`, which allow the use of INLA to run the underlying survival analysis.

## Installation
`survHEhmc` can be installed from this GitHub repository using the package `remotes`:
```R
remotes::install_github("giabaio/survHE", ref="hmc")
```

## Usage
Once `survHEhmc` is available, then you can refer to the whole manual/instructions for `survHE`. For instance, to fit a model using HMC, the following code would work:
```R
# Load survHE
library(survHE)

# Loads an example dataset from 'flexsurv'
data(bc)

# Fits the same model using HMC
hmc = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
   distr="exp",method="hmc")
   
# Prints the results using the survHE method
print(hmc)

# Or visualises the results using the original package methods
print(hmc,original=TRUE)

# Plots the survival curves and estimates
plot(hmc)
```

Basically, the user doesn't even "see" that `survHEhmc` is being used...