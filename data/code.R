# Provisionally set the working directory. Change this as required or delete entirely
setwd("~/Dropbox/UCL/Mapi/Projects/Survival/Documentation/2ndreview/Code/")

# Sets up width of output
options(width = 76)

# Load library
library("survHE")
# NB: This is based on the current development version of the package, which is available for now on the GitHub repo
# so you need to install the current version of the package with
# devtools:install_github("giabaio/survHE",ref="devel")

# Set seed for replication
set.seed(14081973)

# Reads the data in & shows a selection of records
data <- as_tibble(read.csv("data.csv"))
data$arm <- factor(data$arm)
data$sex <- factor(data$sex)

rbind(head(data), tail(data))

# Run a selection of models using the fit.models function in survHE
# First, defines the vector of models to be used
mods <- c("exp", "weibull", "gamma", "lnorm", "llogis", "gengamma") 
# Then the formula specifying the linear predictor
formula <- Surv(time, censored) ~ arm
# And then runs the models using MLE via flexsurv & shows the elements of the
# resulting survHE object
m1 <- fit.models(formula = formula, data = data, distr = mods)
# Explores the output
class(m1)
names(m1)
names(m1$models)
names(m1$models[[1]])
names(m1$models$Exponential)
m1$model.fitting
m1$method

# Or Bayesian analysis via INLA
# NB: INLA can only do a few models, so defines a specialized version of this
# vector

distr.inla <- c("exp", "weibull", "lnorm", "llogis")

# # 30 May 2019: Temporary hack to circumvent INLA throwing an error message.
# # This will be fixed in the new release of INLA 
# assign("enable.model.likelihood.loglogisticsurv", TRUE, 
#         envir  = INLA:::inla.get.inlaEnv())
m2 <- fit.models(formula, data, distr.inla, method = "inla",
                 control.family = list(
                     exp = list(),
                     weibull = list(prior="gamma", param = c(.1,.1)),
                     llogis = list(),
                     lnorm = list(initial = 1)
                 ))

# Or Bayesian analysis using HMC/Stan
m3 <- fit.models(formula, data, distr = mods, "hmc", cores = 2)

# This produces some diagnostic plots for the HMC model, using rstan functions
rstan::traceplot(m3$models[[2]])
rstan::stan_ac(m3$models[[2]])

# This runs the Bayesian model via HMC using the Generalised Gamma and with a specific setting for the prior of the linear predictor 
# naming convention is consistent with survHE (as in Table 1)
m4 <- fit.models(formula, data, distr = "gengamma", method = "hmc", 
                 priors = list(gengamma=list(sigma_beta = rep(5, 2))))

# But can also define priors on a set of models
priors <- list(exp=list(sigma_beta = rep(4, 2)), wei=list(mu_beta = rep(2, 2)))
m4 <- fit.models(formula, data, distr = mods, method = "hmc", priors = priors)
# Or construct a full list of priors & then re-run the models using fit.models
priors <- vector("list", 6)
# But still needs to name the list (at least for the non-null elements)
names(priors)=mods
priors[[6]] <- list(a_sigma = 2, b_sigma = 4)
m4 <- fit.models(formula, data, distr = mods, method = "hmc", priors = priors)

## Visualising the results
# Use the print method for survHE objects. This uses the first model stored in m1, by default
print(m1)
# But can choose to show the fifth model
print(m1, mod = 5)
# Or use the original flexsurv table instead
print(m1, mod = 6, original = TRUE)
# And it works just as well for other inferential engines - this is for HMC, using survHE notation
print(m3, 6)
# Or the original rstan one
print(m3, mod = 2, original = TRUE)

# Can also plot the survival curves for each object, eg.
plot(m1)

# This plots the survival curves for some of the models included in m1 and m3
plot(m1, m3, mods = c(2, 3, 8, 9), colors = c("blue", "green", "red", "yellow"), 
     labs = c("Weibull (MLE)", "Gamma (MLE)", "Weibull (HMC)", "Gamma (HMC)"))

# With more complex models, things are not as straightforward: for example, we could add a covariate:
m5 <- fit.models(Surv(time, censored) ~ arm + sex, 
                 data = data, distr = "exp")
print(m5)
plot(m5)

# In cases like this, better use 'newdata' and selects only a given profile to plot:
m5  = fit.models(Surv(time, censored) ~ as.factor(arm) + age, 
                 data = data, distr = "exp")
newdata <- list(list(arm = 0, age = mean(data$age)), list(arm = 1, age = mean(data$age)))
plot(m5, newdata = newdata)

# The summary method computes the mean survival time
summary(m1)
# We need to make sure, however, that the range of times is long enough as to characterise the full survival curve
summary(m1, t = seq(0, 60))
# And we can also include a specific profile to compute the mean
summary(m1, t = seq(0, 60), newdata = list(list(arm = "1")))

# Another important issue is to do model fitting - we can use survHE function 'model.fit.plot'
model.fit.plot(m1)
# Selecting the type of statistic ("aic" vs "bic" vs "dic")
model.fit.plot(m1, type = "bic")
# And we can fully customise the plot, mixing different objects and selecting which models we want
model.fit.plot(m1, m3, mods=c(1, 2, 3, 7, 8, 9),
               models = c("Exponential (MLE)", "Weibull (MLE)", "Gamma (MLE)", 
                        "Exponential (HMC)", "Weibull (HMC)", "Gamma (HMC)")
)
model.fit.plot(MLE=m1, HMC=m3, stacked = T, mods=c(1, 2, 3, 7, 8, 9),
               models = c("Exponential (MLE)", "Weibull (MLE)", "Gamma (MLE)", 
                        "Exponential (HMC)", "Weibull (HMC)", "Gamma (HMC)"),
               name_legend = "Inferential method"
)

model.fit.plot(m1, m3, stacked = T, mods=c(1, 2, 3, 7, 8, 9), type = "dic",
               models=c("Exponential (MLE)", "Weibull (MLE)", "Gamma (MLE)", 
                        "Exponential (HMC)", "Weibull (HMC)", "Gamma (HMC)")
)

#model.fit.plot(m1, m3, mods = c(1, 2, 3, 7, 8, 9), type = "dic", 
#               mar = c(4, 7, 3, 0.5), xlim = c(1200, 1290), 
#               models = c("Exponential (MLE)", "Weibull (MLE)", "Gamma (MLE)", 
#                          "Exponential (HMC)", "Weibull (HMC)", "Gamma (HMC)")
#)


# PSA - that's crucial to economic evaluation and survHE as well. We can use the function 'make.surv'
psa <- make.surv(fit = m3, nsim = 1000, t = seq(.1, 63))
# Show its elements
names(psa)
# And then visualise the results
psa.plot(psa)
# And fully customise the call to this function
psa.plot(psa, xlab = "Extrapolated time", ylab = "Estimation of the survival curves", 
         alpha = 0.2, col = c("dark grey", "black"), main = "PSA to survival curves",
         xpos = 30, ypos = 1, cex.txt = .95, offset = 1.5, nsmall = 0, digits = 2
)
# Can also use the standard "plot" function to get the same picture...
plot(m3, mods = 1, nsim = 1000, t = seq(.1, 63))

# And finally write the results to an Excel file
write.surv(psa, file = "temp.xls")

# We can also do "advanced" models, ie Royston-Parmar splines & Poly-Weibull
# RPS is available under MLE & Bayesian HMC modelling:
formula <- Surv(time, censored) ~ arm
m6 <- fit.models(formula = formula, data = data, distr = "rps", k = 2)
m7 <- fit.models(formula = formula, data = data, distr = "rps", k = 2, method = "hmc")
# We can use the methods print, plot and summary, eg
print(m6)
print(m7)

# For the Poly-Weibull, we need to specify a list of formulae before running the Bayesian/HMC model
formula.pw <- list(Surv(time, censored) ~ 1, Surv(time, censored) ~ arm)
m8  <- poly.weibull(formula.pw, data, cores = 2)
# But the results aren't great...
print(m8)
# And possibly use more advanced options in the background call to rstan
m9 <- poly.weibull(formula.pw, data, cores = 2, 
                   control = list(adapt_delta = .9, stepsize = .01, max_treedepth = 100)
)
# Now the results are better
print(m9)

# Digitising --- based on some input txt files with data obtained using DigitizeIt
surv_inp <- "survival.txt"
nrisk_inp <- "nrisk.txt"
km_out <- "KMdata.txt"
ipd_out <- "IPDdata.txt"
digitise(surv_inp = surv_inp, nrisk_inp = nrisk_inp, 
         km_output = km_out, ipd_output = ipd_out)
