#' Helper function to create the covariates profile to use in the 
#' computation of the survival curves
#' 
#' @param formula a formula specifying the model to be used, in the form
#' \code{Surv(time,event)~treatment[+covariates]} in flexsurv terms, or
#' \code{inla.surv(time,event)~treatment[+covariates]} in INLA terms.
#' @param data A data frame containing the data to be used for the analysis.
#' This must contain data for the 'event' variable. In case there is no
#' censoring, then \code{event} is a column of 1s.
#' @param newdata a list (of lists), specifiying the values of the covariates
#' at which the computation is performed. For example
#' \code{list(list(arm=0),list(arm=1))} will create two survival curves, one
#' obtained by setting the covariate \code{arm} to the value 0 and the other by
#' setting it to the value 1. In line with \code{flexsurv} notation, the user
#' needs to either specify the value for *all* the covariates or for none (in
#' which case, \code{newdata=NULL}, which is the default). If some value is
#' specified and at least one of the covariates is continuous, then a single
#' survival curve will be computed in correspondence of the average values of
#' all the covariates (including the factors, which in this case are expanded
#' into indicators).
#' @return \item{X}{A matrix with the covariates profile selected to compute
#' the survival curves}.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso make.surv
#' @references Baio (2020). survHE
#' @keywords Parametric survival models
make_profile_surv <- function(formula,data,newdata) {
  # Checks how many elements are given in 'newdata'
  n.elements <- ifelse(is.null(newdata),0,length(newdata))
  n.provided <- unlist(lapply(newdata,function(x) length(x)))
  
  # Temporarily re-writes the model formula to avoid issues with naming
  formula_temp <- update(formula,paste(all.vars(formula,data)[1],"~",all.vars(formula,data)[2],"+."))
  # Creates a tibble with all the covariates in their original format
  covs <- data %>% model.frame(formula_temp,.) %>% as_tibble(.) %>%  select(-c(1:2)) %>% 
    rename_if(is.factor,.funs=~gsub("as.factor[( )]","",.x)) %>% 
    rename_if(is.factor,.funs=~gsub("[( )]","",.x)) %>% 
    bind_cols(as_tibble(model.matrix(formula_temp,data)) %>% select(contains("Intercept"))) %>% 
    select(`(Intercept)`,everything())
  ncovs <- ncol(covs) - 1
  # Selects the subset of categorical covariates
  is.fac <- covs %>% select(where(is.factor))
  nfacts <- covs %>% select(where(is.factor)) %>% with(ncol(.))
  
  # If formula is in 'inla' terms now change it back to 'flexsurv' terms
  formula_temp <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
  # Computes the "average" profile of the covariates
  X <- data %>% model.matrix(formula_temp,.) %>% as_tibble(.) %>% summarise_all(mean) 
  colnames(X)=colnames(covs)
  
  # The way the object X *must* be formatted depends on which way it's been generated.
  # The point is that it *always* has to be a matrix for other functions to process
  if(n.elements==0){
    # If all the covariates are factors, then get survival curves for *all* the combinations
    if(nfacts==ncovs & nfacts>0) {
      X <- apply(covs %>% unique(.),2,as.numeric)
    } else {
      X <- as.matrix(X,nrow=nrow(X),ncol=ncol(X))
    }
  }
  
  # If 'newdata' provides a specific (set of) profile(s), then use that
  if (n.elements>=1) {
    # This means that n.provided will also be > 0 but needs to check it's the right number and if not, stop
    if (!all(n.provided==ncovs)) {
      stop("You need to provide data for *all* the covariates specified in the model, in the list 'newdata'")
    } else {
      # Creates a design matrix containing the information provided in 'newdata'
      X <- bind_rows(lapply(newdata,function(x) as_tibble(x)))
      # But if there's an intercept in the model matrix, then add it
      if ("(Intercept)" %in% colnames(model.matrix(formula,data))) {
        X <- X %>% mutate(`(Intercept)`=1) %>% select(`(Intercept)`,everything())
      }
      X <- as.matrix(X)
    }
  }
  # Returns the output (the matrix with the covariates profile)
  return(X)
}
