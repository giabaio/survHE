# survHE
## Survival analysis in health economic evaluation

Contains a suite of functions to systematise the workflow involving survival analysis in health economic evaluation. survHE can fit a large range of survival models using both a frequentist approach (by calling the R package [flexsurv](https://CRAN.R-project.org/package=flexsurv)) and a Bayesian perspective. For a selected range of models, both Integrated Nested Laplace Integration (via the R package [INLA](http://www.r-inla.org/)) and Hamiltonian Monte Carlo (via the R package [rstan](https://CRAN.R-project.org/package=rstan)) are possible. HMC models are pre-compiled so that they can run in a very efficient and fast way. In addition to model fitting, survHE provides a set of specialised functions, for example to perform Probabilistic Sensitivity Analysis, export the results of the modelling to a spreadsheet, plotting survival curves and uncertainty around the mean estimates.

## Installation
Under windows you need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first.

To ensure that all dependencies are available, run
```R
install.packages(c("flexsurv",
                   "Rcpp",
                   "rms",
                   "xlsx",
                   "data.table",
                   "rstan",
                   "StanHeaders",
                   "BH",
                   "RcppEigen"),
                 dependencies = "Depends")
install.packages("INLA",
                 repos = "https://www.math.ntnu.no/inla/R/stable",
                 dependencies = "Depends")
```

before installing the package using `devtools`:

```R
devtools::install_github("giabaio/survHE")
```
