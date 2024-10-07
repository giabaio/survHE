# survHE 

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/survHE)](https://cran.r-project.org/package=survHE)
[![CRAN_Download_Badge](https://cranlogs.r-pkg.org/badges/survHE)](https://cran.r-project.org/package=survHE)
[![CRAN_Download_Badge](https://cranlogs.r-pkg.org:443/badges/grand-total/survHE?color=orange)](https://cranlogs.r-pkg.org:443/badges/grand-total/survHE?color=orange)

## Survival analysis in health economic evaluation

Contains a suite of functions to systematise the workflow involving survival analysis in health economic evaluation. survHE can fit a large range of survival models using both a frequentist approach (by calling the R package [flexsurv](https://CRAN.R-project.org/package=flexsurv)) and a Bayesian perspective. For a selected range of models, both Integrated Nested Laplace Integration (via the R package [INLA](https://www.r-inla.org/)) and Hamiltonian Monte Carlo (via the R package [rstan](https://CRAN.R-project.org/package=rstan)) are possible. HMC models are pre-compiled so that they can run in a very efficient and fast way. In addition to model fitting, survHE provides a set of specialised functions, for example to perform Probabilistic Sensitivity Analysis, export the results of the modelling to a spreadsheet, plotting survival curves and uncertainty around the mean estimates.

**NB**: To run the Bayesian models, as of version 2.0 of `survHE`, it is necessary to install the additional packages [`survHEinla`](https://github.com/giabaio/survHEinla) and/or [`survHEhmc`](https://github.com/giabaio/survHEhmc), which are available from this GitHub repository. The reason for this structural change is that in this way, the basic backbone of `survHE` (available from this `main` branch of the repo) becomes a very lean package, whose installation is very quick. More details [here](https://gianluca.statistica.it/blog/2022-01-18-survhe-light/). All the functionalities are in place for `survHE` to easily extend to the Bayesian versions, once one or both of the additional "modules" is also installed.

## Installation
The most updated version can be installed using the following code.
```R
install.packages("remotes")
remotes::install_github("giabaio/survHE")
```

To run the Bayesian versions of the models, you also need to install the ancillary packages
```R
# Bayesian models using HMC/Stan
remotes::install_github("giabaio/survHEhmc")
# or alternatively
install.packages("survHEhmc",repos=c("https://giabaio.github.io/drat/","https://www.stats.bris.ac.uk/R/"),type="source")

# Bayesian models using INLA
remotes::install_github("giabaio/survHEinla")
# or alternatively
install.packages("survHEinla",repos=c("https://giabaio.github.io/drat/","https://www.stats.bris.ac.uk/R/"),type="source")
```
(these two are optional, in some sense, so you don't *have* to, unless you want to do the right thing and be Bayesian about it... :wink:)
