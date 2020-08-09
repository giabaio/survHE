#' Helper function to create data in the correct format for rstan
#' 
#' @param formula a formula specifying the model to be used, in the form
#' \code{Surv(time,event)~treatment[+covariates]} in flexsurv terms, or
#' \code{inla.surv(time,event)~treatment[+covariates]} in INLA terms.
#' @param data A data frame containing the data to be used for the analysis.
#' This must contain data for the 'event' variable. In case there is no
#' censoring, then \code{event} is a column of 1s.
#' @return \item{data.stan}{A list containing the variables needed to pass
#' to 'stan' when calling \code{fit.models} with \code{method="hmc"}}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models Bayesian inference via Hamiltonian
#' Monte Carlo 
make_data_stan=function(formula,data,distr3) {
  # Modifies the original formula to separate 'time' and 'event'
  formula_temp <- update(formula,paste(all.vars(formula,data)[1],"~",all.vars(formula,data)[2],"+."))
  # Creates a model.frame + renames the variables to conform with stan's expectations later
  mf <- as_tibble(model.frame(formula_temp,data)) %>% 
    rename(time=1, event=2) %>% rename_if(is.factor,.funs=~gsub("as.factor[( )]","",.x)) %>% 
    rename_if(is.factor,.funs=~gsub("[( )]","",.x)) %>% 
    bind_cols(as_tibble(model.matrix(formula_temp,data)) %>% select(contains("Intercept"))) %>% 
    select(time,event,contains("Intercept"),everything()) %>% tibble::rownames_to_column("ID")

  # Now arrange data in list to pass to 'stan'
  # NB: Need different formatting depending on the underlying sampling distribution
  if (distr3 %in% c("gam","gga","gef")) {
    # If model is Gamma, GenGamma or GenF, then use the "obs vs" censored format
    data.stan <- list(t=(mf %>% filter(event==1))$time,
                      d=(mf %>% filter(event==0))$time,
                      n_obs=mf %>% filter(event==1) %>% with(nrow(.)),
                      n_cens=mf %>% filter(event==0) %>% with(nrow(.))
    )
    data.stan$X_obs <- matrix(model.matrix(formula,data)[(mf %>% filter(event==1))$ID,],nrow=data.stan$n_obs)
    data.stan$X_cens <- matrix(model.matrix(formula,data)[(mf %>% filter(event==0))$ID,],nrow=data.stan$n_cens)
    data.stan$H=ncol(data.stan$X_obs)

    # NB: Stan doesn't allow vectors of size 1, so if there's only one covariate (eg intercept only), needs a little trick
    if (data.stan$H==1) {
      data.stan$X_obs <- cbind(data.stan$X_obs,rep(0,data.stan$n_obs))
      data.stan$X_cens <- cbind(data.stan$X_cens,rep(0,data.stan$n_cens))
      data.stan$H <- ncol(data.stan$X_obs)
    }
  }

  if (distr3 %in% c("exp", "gom", "wei", "wph", "llo", "lno")) {
    # If it's one of the others (except polyweibull), use the "h,S" format
    data.stan <- list(t=(mf$time),
                      d=mf$event,
                      n=nrow(mf),
                      X=matrix(model.matrix(formula,data),nrow=nrow(mf)),
                      H=ncol(model.matrix(formula,data))
    )
    # NB: Stan doesn't allow vectors of size 1, so if there's only one covariate (eg intercept only), needs a little trick
    if (data.stan$H==1) {
      data.stan$X <- cbind(data.stan$X,rep(0,data.stan$n))
      data.stan$H <- ncol(data.stan$X)
    }
  }
  
  if (distr3=="rps"){
    # If it's Royston-Parmar splines, then gets the correct data 
    knots <- quantile(log((mf %>% filter(event==1))$time), seq(0,1,length=k+2))
    # Uses flexsurv to compute the basis and derivatives of the basis
    B <- flexsurv::basis(knots,log(mf$time))
    DB <- flexsurv::dbasis(knots,log(mf$time))
    # Now checks to see whether the user wants to specify covariates and removes the intercept from the formula (for identifiability)
    mm <- model.matrix(formula,data)[,-1]
    # a. if the formula is ~ 1, then adds two fictional covariates of all 0s
    if (length(mm)<1) {
      mm <- matrix(rep(0,nrow(mf)),nrow=nrow(mf),ncol=2)
    }
    # b. in case there's only one covariate, then adds another fake covariate of all 0s
    if (is.null(dim(mm))) {
      mm <- cbind(mm,rep(0,length(mm)))
    }
    data.stan <- list(t=mf$time,
                      d=mf$event,
                      n=nrow(df),
                      M=k,
                      X=mm,
                      H=ncol(mm),
                      B=B,
                      DB=DB,
                      mu_gamma=rep(0,k+2),
                      sigma_gamma=rep(5,k+2),
                      knots=knots
    )
  }
  data.stan
}
