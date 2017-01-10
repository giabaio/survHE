# survHE [![Travis-CI Build Status](https://travis-ci.org/giabaio/survHE.svg?branch=master)](https://travis-ci.org/giabaio/survHE)[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/giabaio/survHE?branch=master&svg=true)](https://ci.appveyor.com/project/giabaio/survHE)
## Survival analysis in health economic evaluation

Contains a suite of functions to systematise the workflow involving survival analysis in health economic evaluation. survHE can fit a large range of survival models using both a frequentist approach (by calling the R package [flexsurv](https://CRAN.R-project.org/package=flexsurv)) and a Bayesian perspective. For a selected range of models, both Integrated Nested Laplace Integration (via the R package [INLA](http://www.r-inla.org/)) and Hamiltonian Monte Carlo (via the R package [rstan](https://CRAN.R-project.org/package=rstan)) are possible. HMC models are pre-compiled so that they can run in a very efficient and fast way. In addition to model fitting, survHE provides a set of specialised functions, for example to perform Probabilistic Sensitivity Analysis, export the results of the modelling to a spreadsheet, plotting survival curves and uncertainty around the mean estimates.

## Installation
There are two ways of installing `survHE`. A "stable" version is packaged and binary files are available for Windows and as source. To install the stable version on a Windows machine with x64 architecture, run the following commands
```R
install.packages("http://www.statistica.it/gianluca/survHE/x64/survHE_1.0.4.zip",repos=NULL)
deps <- tools::package_dependencies("survHE", db = installed.packages())[[1]]
pkgs <- deps[!deps %in% installed.packages()[,1]]
if (length(pkgs) > 0) 
  install.packages(pkgs = pkgs,
    repos = c("https://cran.rstudio.com", "https://www.math.ntnu.no/inla/R/stable"),
	dependencies="Depends")
rm(deps, installed.packages(), pkgs)
```
which first install `survHE` and then the relevant "dependencies". This process can be quite lengthy. On a Windows machine with a i386 architecture, simply replace the first line with
```R
install.packages("http://www.statistica.it/gianluca/survHE/i386/survHE_1.0.4.zip",repos=NULL)
```
Finally, to install from source (e.g. on a Linux machine), run
```R
install.packages("http://www.statistica.it/gianluca/survHE/src/survHE_1.0.4.tar.gz",type="source",repos=NULL)
```

The second way involves using the "development" version of `survHE` - this will usually be updated more frequently and may be continuously tested. On Windows machines, you need to install a few dependencies, including [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first, e.g. by running
```R
pkgs <- c("flexsurv","Rcpp","rms","xlsx","rstan","INLA","Rtools","devtools")
repos <- c("https://cran.rstudio.com", "https://www.math.ntnu.no/inla/R/stable") 
install.packages(pkgs,repos=repos,dependencies = "Depends")
```
before installing the package using `devtools`:
```R
devtools::install_github("giabaio/survHE")
```
Under Linux or MacOS, it is sufficient to install the package via `devtools`:
```R
devtools:install_github("giabaio/survHE")
```
