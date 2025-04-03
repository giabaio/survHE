#' make.transition.probs
#' 
#' Computes the transition probabilities (to be passed to a Markov model) from
#' the cumulative hazard curves obtained using \code{fit.models}, using the formula 
#' p(t)=1-exp(H(t-k)/H(t)), where k is the Markov model cycle length (or the 
#' difference across two consecutive times) and t is a generic time
#' 
#' @aliases make.transition.probs
#' @param fit an object obtained as output of the call to \code{fit.models}
#' @param labs a vector with labels to identify the 'profiles' ie the 
#' combination of covariates that have been passed onto the model formula.
#' If 'NULL' (default), then figures it out from the 'survHE' object.
#' @param ...  additional arguments. Includes the standard inputs to the 
#' call to \code{make.surv}, so \code{mod} (the index of the possibly many
#' models stored in the 'survHE' object), \code{t} (the vector of times 
#' over which to compute the survival curves), \code{newdata} (a list that
#' defines the profile of covariates) and \code{nsim} (the number of 
#' simulations to use - default is \code{nsim}=1)
#' @return A tibble 'lambda' with an indicator for the treatment arm,
#' the times at which the probabilities have been computed and \code{nsim}
#' columns each with a simulation of the transition probabilities for 
#' all the times specified by the user
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso \code{\link{make.surv}}
#' @references Something will go here
#' @keywords Transition probabilities Markov models
#' @examples
#' \dontrun{
#' # Something will go here
#' }
#' 
#' @export make.transition.probs
make.transition.probs <- function(fit,labs=NULL,...) {
  exArgs <- list(...)
  
  # Makes default for parameters to the call to 'make.surv' (which are overwritten if the user has specified them
  # separately)
  if(exists("mod",exArgs)) {mod <- exArgs$mod} else {mod <- 1}
  if(exists("t",exArgs)) {t <- exArgs$t} else {t <- NULL}
  if(is.null(t)) {
    t <- sort(unique(fit$misc$km$time))
    # Add an extra time=0 at the beginning. This ensures the computation can be done for all the actual times
    # specified by the user (because the first one has no lag and so the ratio gives NA...)
    #if(t[1]>0) {t=c(0,t)}
  }
  # Add an extra time=0 at the beginning. This ensures the computation can be done for all the actual times
  # specified by the user (because the first one has no lag and so the ratio gives NA...)
  t <- c(0,t)
  if(exists("newdata",exArgs)) {newdata <- exArgs$newdata} else {newdata <- NULL}
  if(exists("nsim",exArgs)) {nsim <- exArgs$nsim} else {nsim <- 1}
  # Now computes the simulations using 'make.surv'
  s <- make.surv(fit,mod=mod,t=t,newdata=newdata,nsim=nsim)
  
  # Get labels of the 'strata' 
  strata <- lapply(1:nrow(s$des.mat),function(x){
    s$des.mat %>%
      as_tibble() %>%
      select(!contains("(Intercept)")) %>%
      slice(x) %>% 
      round(digits=2) %>%
      mutate(strata=paste0(names(.),"=",.,collapse=","))
    }) %>%
    bind_rows(.) %>%
    select(strata) %>%
    pull(strata)
  
  # Now retrieves the transition probabilities 'lambda' applying the formula
  # lambda(t)=1-S(t+k)/S(t) where k is the MM cycle length and t is a generic time
  lambda <- s$mat %>% 
    # First stacks together all the matrices with the simulation(s) for the survival curves
    bind_rows() %>% 
    # Then creates a treatment indicator (1,2,...) & place it as the first column
    mutate(profile=rep(strata,each=nrow(s$mat[[1]]))) %>% select(profile,everything()) %>% 
    # Then group by treatment and computes the cumulative hazards
    group_by(profile) %>% mutate(across(starts_with("S"),~-log(.))) %>% 
    # Then computes the probabilities using 1-exp(H(t-k)-H(t)), where k is the cycle length (or difference across times)
    ####mutate(across(starts_with("S"),~case_when(.==0~0,TRUE~1-exp(lag(.)-.)))) %>% 
    mutate(across(starts_with("S"), ~ 1-exp(lag(.)-.))) %>% 
    # Then removes the first row (in each treatment) - that was just artificially added anyway...
    slice(-1) %>% ungroup() 

  # And now renames the columns from S(_1,S_2,...,S_nsim) to lambda(_1,lambda_2,...,lambda_nsim)
  if (nsim==1) {
    lambda <- lambda %>%
      rename_with(starts_with("S"), .fn=~"lambda")
  } else {
    lambda <- lambda %>%
      rename_with(starts_with("S"), .fn=~paste0("lambda_",1:nsim))
  }

  lambda
}


#' make_data_multi_state
#' 
#' Takes as input an individual-level dataset including data on both
#' progression and death time (**jointly**) and manipulates it using
#' \code{dplyr} functions to create a full "multi-state" dataset, in
#' which all the transitions are tracked. This can then be used
#' to fit survival models and compute all the estimates for the
#' whole set of transition probabilities
#' 
#' @aliases make_data_multi_state
#' @param data dataset containing the full ILD with information on both
#' progression and death. Can be a data.frame or a tibble
#' @param id The column with the individual identifier. Can be NULL (in
#' which case, it will be created from scratch)
#' @param prog The progression indicator: takes value 1 if the individual
#' has progressed and 0 otherwise. Defaults to the column named 'prog' in
#' the dataset
#' @param death The death indicator: takes value 1 if the individual
#' has died and 0 otherwise. Defaults to the column named 'death' in
#' the dataset
#' @param prog_t The progression time. Defaults to the column named 
#' 'prog_t' in the dataset
#' @param death_t The death time. Defaults to the column named 
#' 'death_t' in the dataset
#' @param keep A vector of strings with the names of the additional 
#' variables from the original dataset to keep into the multistate
#' dataset. If 'NULL' (default), then keeps all
#' @param ...  additional arguments. 
#' @return A tibble containing the event history for each individual
#' and with the following variables: id = Patients ID; from = Initial 
#' state (1=Pre-progression, 2=Progression, 3=Death); to = End state 
#' (1=Pre-progression, 2=Progression, 3=Death); trans = Transition ID:
#' 1=Pre-progression -> Progression; 2=Pre-Progression -> Death;
#' 3=Progression -> Death; Tstart = Entry time (either entry or 
#' progression); Tstop = Exit time (time of event or censoring time); 
#' status = Event indicator (1=yes, 0=censored), **for the specific
#' event under consideration**; treat = Treatment indicator
#' All the other original variables are appended to these, but can be
#' removed
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso Something will go here
#' @references Something will go here
#' @keywords Transition probabilities Markov models Multistate models
#' @examples
#' \dontrun{
#' # Something will go here
#' }
#' 
#' @export make_data_multi_state
make_data_multi_state=function(data,id="id",prog="prog",death="death",prog_t="prog_t",
                               death_t="death_t",keep=NULL,...) { 
  
  # If ID is not passed (=NULL) then create it
  if(is.null(id)) {
    data=data %>% mutate(id=row_number())
    id="id"
  } 
  
  # Uses dplyr to manipulate the original dataset and create the mstate version
  # NB: uses the notation '!!sym(name_variable)' to address a specific column in the dataset
  msmdata=
    # Transition Pre to Post
    data %>% mutate(
      id=!!sym(id),                         # patient ID
      from=1,                               # starting state
      to=2,                                 # arriving state 
      trans=1,                              # transition code (1 = Pre -> Progression) 
      Tstart=0,                             # entry time
      Tstop=!!sym(prog_t),                  # exit time
      time=Tstop-Tstart,                    # time-to-event = Tstop-Tstart
      status=case_when(                     # censoring indicator: 
        !!sym(prog)==1~1,                   #  1 if progressed; 0 otherwise
        TRUE~0
      )
    ) %>% select(id,from,to,trans,Tstart,Tstop,time,status,treat,everything()) %>% 
    bind_rows(
      # Transition Pre to Death
      data %>% mutate(
        id=!!sym(id),                         # patient ID
        from=1,                               # starting state
        to=3,                                 # arriving state 
        trans=2,                              # transition code (2 = Pre -> Death)  
        Tstart=0,                             # entry time
        Tstop=!!sym(death_t),                 # exit time
        time=Tstop-Tstart,                    # time-to-event = Tstop-Tstart
        status=case_when(                     # censoring indicator:
                                              #  1 if died at progression; 0 otherwise
          (!!sym(death)==1 & !!sym(prog_t)==!!sym(death_t))~1,     
          TRUE~0
        )
      )  %>% select(id,from,to,trans,Tstart,Tstop,time,status,treat,everything()) 
    ) %>% 
    bind_rows(
      # Transition Post to Death
      data %>% filter(!!sym(prog)==1) %>% mutate(
        id=!!sym(id),                         # patient ID
        from=2,                               # starting state
        to=3,                                 # arriving state 
        trans=3,                              # transition code (2 = Pre -> Death)  
        Tstart=!!sym(prog_t),                 # entry time
        Tstop=!!sym(death_t),                 # exit time
        time=Tstop-Tstart,                    # time-to-event = Tstop-Tstart
        status=case_when(                     # censoring indicator: 
          !!sym(death)==1~1,                  #  1 if died; 0 otherwise
          TRUE~0
        )
      ) %>% select(id,from,to,trans,Tstart,Tstop,time,status,treat,everything()) 
    ) %>% arrange(id,trans)
  
  if (!is.null(keep)) {
    msmdata=msmdata %>% select(id,from,to,trans,Tstart,Tstop,time,status,keep)
  }
  return(msmdata)
}


#' three_state_mm
#' 
#' General purpose function to run a standard three-state Markov model
#' (typically used in cancer modelling). The states are typically
#' 'Pre-progression', 'Progressed' and 'Death'. No backward transition
#' from 'Progressed' to 'Pre-progression' is allowed and 'Death' is 
#' obviously an absorbing state. All other transitions are possible.
#' The crucial assumption is that *individual-level data* are available
#' recording an indicator and the time of progression and death for each
#' individual. The function returns the full transition matrix
#' 
#' @aliases three_state_mm
#' @param m_12 A 'survHE' object (output to a call to \code{fit.models})
#' estimating the parameters of a model for the transition from 
#' 'Pre-progression' (state 1) to 'Progressed' (state 2). Given the 
#' individual level data with the complete event history (in the object
#' 'data'), can be done with a call like 'x=make_data_multi_state(data)'
#' and then \code{fit.models(Surv(time,status)~...,data=x \%>\% filter(trans==1),...)}
#' @param m_13 A 'survHE' object (output to a call to \code{fit.models})
#' estimating the parameters of a model for the transition from 
#' 'Pre-progression' (state 1) to 'Death' (state 3).  Given the 
#' individual level data with the complete event history (in the object
#' 'data'), can be done with a call like 'x=make_data_multi_state(data)'
#' and then \code{fit.models(Surv(time,status)~...,data=x \%>\% filter(trans==2),...)}
#' @param m_23 A 'survHE' object (output to a call to \code{fit.models})
#' estimating the parameters of a model for the transition from 
#' 'Progressed' (state 2) to 'Death' (state 3).  Given the 
#' individual level data with the complete event history (in the object
#' 'data'), can be done with a call like 'x=make_data_multi_state(data)'
#' and then \code{fit.models(Surv(time,status)~...,data=x \%>\% filter(trans==3),...)}
#' @param nsim The number of simulations for the model parameters that are 
#' used to compute the survival curves. Defaults to \code{nsim}=1,
#' which simply creates one survival curve for each treatment arm.
#' @param start A vector of initial state occupancy. By default assumes 1000 
#' individuals, all initially allocated to 'Pre-progression'
#' @param basecase Should the base case be computed as well, based on the 
#' point estimate of the underlying model parameters? (Default=FALSE)
#' @param ...  additional arguments. 
#' @return A list including the state occupancy simulations in an object 'm'.
#' This is a tibble with the number of individuals in each of the 3 states 
#' at each of the times specified by the user. If \code{nsim}>1, then the tibble
#' also contains a simulation index to keep track of that. The list also 
#' includes the computation time to obtain the state occupancy tibble (in the
#' object 'running_time'). If \code{basecase==TRUE}, then the function also
#' computes the "base case scenario" (based on 1 simulation from of the 
#' underlying survival curves, i.e. the point estimate of the model parameters)
#' and stores it in the object 'base_case'
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso make.transition.probs make_data_multi_state
#' @references Something will go here
#' @keywords Transition probabilities Markov models Three-state cancer model
#' @examples
#' \dontrun{
#' # Something will go here
#' }
#' 
#' @export three_state_mm
three_state_mm = function(m_12,m_13,m_23,nsim=1,start=c(1000,0,0),basecase=FALSE,...){
  
  exArgs <- list(...)
  
  # Initialises the base_case object
  base_case=NULL
  
  # Makes default for parameters to the call to 'make.surv' (which are overwritten if the user has specified them
  # separately)
  if(exists("mod",exArgs)) {mod=exArgs$mod} else {mod=1}
  if(exists("t",exArgs)) {t=exArgs$t} else {t=NULL}
  if(is.null(t)) {
    t <- sort(unique(m_12$misc$km$time))
  }
  if(exists("newdata",exArgs)) {newdata=exArgs$newdata} else {newdata=NULL}
  
  # Computes the transition probabilities for the transitions that are directly
  # identifiable from the observed data
  lambda_12=make.transition.probs(m_12,mod=mod,t=t,newdata=newdata,nsim=nsim)
  lambda_13=make.transition.probs(m_13,mod=mod,t=t,newdata=newdata,nsim=nsim)
  lambda_23=make.transition.probs(m_23,mod=mod,t=t,newdata=newdata,nsim=nsim)
  
  # Derives lambda_11 by subtraction (as all transition probs out of state 1 must sum to 1)
  # NB: Without further checks, it is possible that (lambda_12+lambda_13)>1 and so lambda_11<0
  #     This should be carefully checked!
  lambda_11=(lambda_12 %>% select(starts_with("lambda")) + lambda_13 %>% select(starts_with("lambda"))) %>% 
    as_tibble() %>% bind_cols(lambda_12 %>% select(profile,time)) %>% select(profile,time,everything()) %>% 
    mutate(across(starts_with("lambda"),~1-.))
  # Derives lambda_22 by subtraction (as all transition probs out of state 2 must sum to 1)
  lambda_22=(1-lambda_23 %>% select(starts_with("lambda"))) %>% as_tibble() %>% 
    bind_cols(lambda_23 %>% select(profile,time)) %>% select(profile,time,everything())

  # Computes the state occupancy for all the simulations
  tic=Sys.time()
  m=make_state_occupancy(nsim=nsim,lambda_11,lambda_12,lambda_13,lambda_22,lambda_23,start)
  toc=Sys.time()
  running_time=toc-tic

  if(basecase) {
    # Makes base-case Markov model (by considering the point estimate of the model parameters)
    # Computes the transition probabilities for the transitions that are directly
    # identifiable from the observed data
    lambda_12=make.transition.probs(m_12,mod=mod,t=t,newdata=newdata,nsim=1)
    lambda_13=make.transition.probs(m_13,mod=mod,t=t,newdata=newdata,nsim=1)
    lambda_23=make.transition.probs(m_23,mod=mod,t=t,newdata=newdata,nsim=1)
    
    # Derives lambda_11 by subtraction (as all transition probs out of state 1 must sum to 1)
    lambda_11=(lambda_12 %>% select(starts_with("lambda"))+lambda_13 %>% select(starts_with("lambda"))) %>% 
      as_tibble() %>% bind_cols(lambda_12 %>% select(profile,time)) %>% select(profile,time,everything()) %>% 
      mutate(across(starts_with("lambda"),~1-.))
    # Derives lambda_22 by subtraction (as all transition probs out of state 2 must sum to 1)
    lambda_22=(1-lambda_23 %>% select(starts_with("lambda"))) %>% as_tibble() %>% 
      bind_cols(lambda_23 %>% select(profile,time)) %>% select(profile,time,everything())
    
    base_case=make_state_occupancy(nsim=1,lambda_11,lambda_12,lambda_13,lambda_22,lambda_23,start)
  }
  
  
  # Outputs of the function
  list(m=m,running_time=running_time,base_case=base_case)
}


#' make_state_occupancy
#' 
#' Utility function to compute the state occupancy in a three state MM
#' 
#' @param nsim The number of simulations for the model parameters that are 
#' used to compute the survival curves. 
#' @param lambda_11 the tibble containing the transition probabilities for
#' the transition "Pre-progression -> Pre-progression" 
#' @param lambda_12 the tibble containing the transition probabilities for
#' the transition "Pre-progression -> Progression" 
#' @param lambda_13 the tibble containing the transition probabilities for
#' the transition "Pre-progression -> Death" 
#' @param lambda_22 the tibble containing the transition probabilities for
#' the transition "Progression -> Progression" 
#' @param lambda_23 the tibble containing the transition probabilities for
#' the transition "Progression -> Death" 
#' @return A list including the state occupancy simulations in an object 'm'.
#' This is a tibble with the number of individuals in each of the 3 states 
#' at each of the times specified by the user. If \code{nsim}>1, then the tibble
#' also contains a simulation index to keep track of that.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso make.transition.probs make_data_multi_state
#' @references Something will go here
#' @keywords Transition probabilities Markov models Three-state cancer model
#' @examples
#' \dontrun{
#' # Something will go here
#' }
#' 
#' @noRd 
make_state_occupancy=function(nsim,lambda_11,lambda_12,lambda_13,lambda_22,lambda_23,start) {

  m=list()
  # Initialises the lists
  m=lapply(1:nsim,function(i) {
    # Creates the tibbles (one per each of the nsim simulations)
    m[[i]]=tibble(profile=lambda_11$profile,time=lambda_11$time,`Pre-progressed`=NA,Progressed=NA,Death=NA)
    # Now adds in the relevant transition probabilities in the correct rows
    m[[i]]=m[[i]] %>% left_join(lambda_11 %>% select(profile,time,starts_with("lambda")[[i]]) %>% rename("lambda_11"=starts_with("lambda")), by=c("profile","time")) %>% 
      left_join(lambda_12 %>% select(profile,time,starts_with("lambda")[[i]]) %>% rename("lambda_12"=starts_with("lambda")), by=c("profile","time")) %>% 
      left_join(lambda_13 %>% select(profile,time,starts_with("lambda")[[i]]) %>% rename("lambda_13"=starts_with("lambda")), by=c("profile","time")) %>% 
      left_join(lambda_22 %>% select(profile,time,starts_with("lambda")[[i]]) %>% rename("lambda_22"=starts_with("lambda")), by=c("profile","time")) %>% 
      left_join(lambda_23 %>% select(profile,time,starts_with("lambda")[[i]]) %>% rename("lambda_23"=starts_with("lambda")), by=c("profile","time"))
    # Initialise the tibbles with the start values
    m[[i]]=m[[i]] %>% group_by(profile) %>% mutate(
      `Pre-progressed`=replace(`Pre-progressed`,row_number()==1,start[1]),
      Progressed=replace(Progressed,row_number()==1,start[2]),
      Death=replace(Death,row_number()==1,start[3]),
    ) %>% ungroup() 

    # Re-compute the probabilities to *always* sum up to 1 (this gets rid of weird cases where two
    # probs sum to something bigger than 1, which generates negative values)
    m[[i]]=m[[i]] %>% mutate(check=(lambda_12+lambda_13)>1) %>% mutate(
      l_12=case_when(
        check==TRUE~lambda_12/(lambda_12+lambda_13),
        check==FALSE~lambda_12
      ),
      l_13=case_when(
        check==TRUE~lambda_13/(lambda_12+lambda_13),
        check==FALSE~lambda_13
      )
    ) %>% mutate(l_11=1-(l_12+l_13)) %>% mutate(lambda_11=l_11,lambda_12=l_12,lambda_13=l_13) %>%
      select(-c(l_11,l_12,l_13,check))
  })
  
  # Convert the list into a massive tibble with a simulation index
  m=m %>% bind_rows() %>% mutate(sim_idx=rep(1:nsim,each=nrow(m[[1]])))
  
  # Loops over times to compute the state occupancy
  for (j in 2:nrow(m)) {
    # This is a simple trick to restart the computation when the treatment indicator changes
    if (m$profile[j]==m$profile[j-1]) {
      m$`Pre-progressed`[j]=sum(c(m$`Pre-progressed`[j-1]*m$lambda_11[j]),na.rm=T)
      m$Progressed[j]=sum(c(m$`Pre-progressed`[j-1]*m$lambda_12[j], m$Progressed[j-1]*m$lambda_22[j]),na.rm=T)
      m$Death[j]=sum(c(m$`Pre-progressed`[j-1]*m$lambda_13[j], m$Progressed[j-1]*m$lambda_23[j], m$Death[j-1]),na.rm=T)
    }
  }
  m=m %>% select(profile,time,`Pre-progressed`,Progressed,Death,sim_idx,everything())
  # Need to set to 0 all negative values!
  m=m %>% mutate(
    `Pre-progressed`=if_else(`Pre-progressed`<0,0,`Pre-progressed`),
    Progressed=if_else(Progressed<0,0,Progressed),
    Death=if_else(Death<0,0,Death)
  )
  
  return(m)
}


#' Markov trace
#' 
#' Plots the Markov Trace from an object generated using \code{three_state_mm}
#' 
#' @aliases markov_trace
#' @param mm  The output of a call to \code{three_state_mm}
#' @param interventions A vector of labels for the interventions
#' @param ...  additional arguments. 
#' @return Plot
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso make.surv, three_state_mm
#' @references Something will go here
#' @keywords Transition probabilities Markov models Markov trace
#' @examples
#' \dontrun{
#' # Something will go here
#' }
#' 
#' @export markov_trace
markov_trace <- function(mm, interventions=NULL,...) {
  # First reshape the data
  if(!is.null(interventions)) {
    # Figures out how many observations there are in each treatment & replaces the values passed 
    # as arguments in 'interventions'
    mm$m <- mm$m %>%
      mutate(profile = rep(interventions,each=mm$m %>%
                             count(profile) %>%
                             slice(1) %>%
                             pull(n)))
  }
  
  pl=mm$m |> group_by(profile,time) |> 
    summarise(
      `Pre-progressed`=mean(`Pre-progressed`),Progressed=mean(Progressed),
      Death=mean(Death)
    ) |> 
    ungroup()
  pl <- pl |>
    select(profile,time,`Pre-progressed`) %>%
    rename(npeople=`Pre-progressed`) %>%
    mutate(group="Pre-progressed") %>%
    bind_rows(pl %>% select(profile,time,`Progressed`) %>%
                rename(npeople=Progressed) %>%
                mutate(group="Progressed")) %>%
    bind_rows(pl %>% select(profile,time,Death) %>%
                rename(npeople=Death) %>%
                mutate(group="Death")) %>%
    # Create a numeric/factor group label to help manage the appearance of the graph
    mutate(
      grp_lab=as.factor(case_when(
        group=="Pre-progressed"~3,
        group=="Progressed"~2,
        TRUE~1
      ))
    ) %>% 
    ggplot(aes(x=time, y=npeople, fill=grp_lab)) +
    geom_bar(position="stack",stat="identity") +
    labs(x="Cycle",y="Number of people",title="Markov trace",fill="State") +
    facet_wrap(~profile) +
    theme_bw() +
    # Add control to the legend label
    scale_fill_discrete(breaks=c(3,2,1),labels=c("Pre-progressed","Progressed","Death"))

  return(pl)
}


### TO DO --- THIS NICELY COMPUTES THE MEAN SURVIVAL TIME OVER THE MM SIMULATIONS BUT CURRENTLY ONLY
###           DOES IT FOR WHEN nsim=1!!!
compute_mean_time=function(x,disc.rate=0) {
  mean_time=x %>% mutate(Alive=(1000-Death)/1000,disc=1/(1+disc.rate)^time) %>% group_by(profile) %>% 
    mutate(dt=time-lag(time)) %>%  summarise(mean_survival_time=sum(Alive*dt*disc,na.rm=T)) %>% ungroup()
  
  return(mean_time)
}

