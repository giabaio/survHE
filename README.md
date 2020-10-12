# survHE [![Travis-CI Build Status](https://travis-ci.org/giabaio/survHE.svg?branch=master)](https://travis-ci.org/giabaio/survHE)[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/giabaio/survHE?branch=master&svg=true)](https://ci.appveyor.com/project/giabaio/survHE)[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/survHE)](https://cran.r-project.org/package=survHE)[![CRAN_Download_Badge](http://cranlogs.r-pkg.org/badges/survHE)](https://cran.r-project.org/package=survHE)[![CRAN_Download_Badge](https://cranlogs.r-pkg.org:443/badges/grand-total/survHE?color=orange)](https://cranlogs.r-pkg.org:443/badges/grand-total/survHE?color=orange)
## Survival analysis in health economic evaluation

Contains a suite of functions to systematise the workflow involving survival analysis in health economic evaluation. survHE can fit a large range of survival models using both a frequentist approach (by calling the R package [flexsurv](https://CRAN.R-project.org/package=flexsurv)) and a Bayesian perspective. For a selected range of models, both Integrated Nested Laplace Integration (via the R package [INLA](http://www.r-inla.org/)) and Hamiltonian Monte Carlo (via the R package [rstan](https://CRAN.R-project.org/package=rstan)) are possible. HMC models are pre-compiled so that they can run in a very efficient and fast way. In addition to model fitting, survHE provides a set of specialised functions, for example to perform Probabilistic Sensitivity Analysis, export the results of the modelling to a spreadsheet, plotting survival curves and uncertainty around the mean estimates.

## Installation
There are two ways of installing `survHE`. A "stable" version is packaged and binary files are available for Windows and as source. To install the stable version on a Windows machine, run the following commands
```R
install.packages("survHE",
	repos=c("http://www.statistica.it/gianluca/R",
		"https://cran.rstudio.org",
                "https://inla.r-inla-download.org/R/stable"),
	dependencies=TRUE
)
```
Note that you need to specify a vector of repositories - the first one hosts `survHE`, while the second one should be an official [CRAN mirror](https://cran.r-project.org/index.html). You can select whichever one you like, but a CRAN mirror must be provided, so that `install.packages()` can also install the "dependencies" (e.g. other packages that are required for `survHE` to work). The third one is used to install the package [`INLA`](http://www.r-inla.org/), which is used to do one version of the Bayesian analysis. This process can be quite lengthy, if you miss many of the relevant packages.

To install from source (e.g. on a Linux machine), run
```R
install.packages("survHE",
	repos=c("http://www.statistica.it/gianluca/R",
		"https://cran.rstudio.org",
		"https://inla.r-inla-download.org/R/stable"),
	type="source",
	dependencies=TRUE
)
```

The second way involves using the GitHub version of `survHE` - this will usually be updated more frequently and may be continuously tested. On Windows machines, you need to install a few dependencies, including [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first, e.g. by running
```R
pkgs <- c("flexsurv","Rcpp","rms","xlsx","rstan","INLA","Rtools","devtools")
repos <- c("https://cran.rstudio.com", "https://inla.r-inla-download.org/R/stable") 
install.packages(pkgs,repos=repos,dependencies = "Depends")
```
before installing the package using `devtools`:
```R
devtools::install_github("giabaio/survHE")
```
Under Linux or MacOS, it is sufficient to install the package via `devtools`:
```R
install.packages("devtools")
devtools:install_github("giabaio/survHE")
```

The current GitHub `master` version of `survHE` is aligned with the CRAN release. We'll update it to fix minor issues *if and when* we discover them, so we recommend that users install this to ensure a smoother workflow. 

There is also a **development** version of `survHE`, which is available under [https://github.com/giabaio/survHE/tree/devel](https://github.com/giabaio/survHE/tree/devel). This will be continously updated, but it will be less stable than the `master` version, as we will use it to test new functionalities. We recommend you use it only if you are an experienced `R` user.
