# survHEinla [![Travis-CI Build Status](https://travis-ci.org/giabaio/survHE.svg?branch=inla)](https://travis-ci.org/giabaio/survHE)[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/giabaio/survHE?branch=inla&svg=true)](https://ci.appveyor.com/project/giabaio/survHE)
## Survival analysis in health economic evaluation using INLA

This is a module to complement the package `survHE` and expand its functionalities to run survival analysis in health economic evaluation from a Bayesian perspective, using Integrated Nested Laplace Integration (via the R package [INLA](http://www.r-inla.org/)). `survHEinla` "depends" on the main installation of `survHE`. This means that you shouldn't use `survHEinla` as a standalone package --- rather you use all the functions of `survHE` (to fit the models and the post-process the results); installing `survHEinla` basically opens up a new option in the `survHE` function `fit.models`, which allow the use of INLA to run the underlying survival analysis.

## Installation
`survHEinla` can be installed from this GitHub repository using the package `remotes`:
```R
remotes::install_github("giabaio/survHE", ref="inla")
```

## Usage
Once `survHEinla` is available, then you can refer to the whole manual/instructions for `survHE`. For instance, to fit a model using INLA, the following code would work:
```R
# Load survHE
library(survHE)

# Loads an example dataset from 'flexsurv'
data(bc)
     
# Fits the same model using INLA
# NB if survHEinla is installed, then the option 'method="inla"' automatically
#    loads it up in the background
inla = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
         distr="exp",method="inla")

# Prints the results in comparable fashion using the survHE method
print(inla)

# Or visualises the results using the original package methods
print(inla,original=TRUE)

# Or plots the survival curves and estimates
plot(inla)
```

Basically, the user doesn't even "see" that `survHEinla` is being used...