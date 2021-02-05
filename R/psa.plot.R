#' Graphical depiction of the probabilistic sensitivity analysis for the
#' survival curves
#' 
#' Plots the survival curves for all the PSA simulations. The function is 
#' actually deprecated - similar graphs can be obtained directly using 
#' the \code{plot} method (with options), which allows a finer depiction
#' of the results. 
#' 
#' 
#' @param psa the result of the call to the function \code{make.surv}
#' @param ...  Optional graphical parameters, such as: \code{xlab} = label for
#' the x-axis \code{ylab} = label for the y-axis \code{col} = (vector) of
#' colors for the lines to be plotted \code{alpha} = the level of transparency
#' for the curves (default = 0.2)
#' @author Gianluca Baio
#' @template refs
#' @seealso \code{make.surv}, \code{write.surv}
#' @keywords Survival models Bootstrap Probabilistic sensitivity analysis
#' @examples
#' \dontrun{ 
#' data(bc)
#' 
#' # Fits the same model using the 3 inference methods
#' mle = fit.models(formula=Surv(recyrs,censrec)~group,data=bc,
#'     distr="exp",method="mle")
#' p.mle = make.surv(mle,nsim=100)
#' psa.plot(p.mle)
#' }
#' 
#' @export psa.plot
psa.plot <- function(psa,...) {
  # Plots the survival curves for all the PSA simulations
  # psa = the result of the call to the function make.surv
  # ... = additional arguments
  # xlab = label for the x-axis
  # ylab = label for the y-axis
  # col = vector of colours with which to plot the curves
  # alpha = parameter to determine the transparency (default = 0.2)
  # main = a string to write the title
  # labs = a vector with non-standard names for the legend values
  # name_labs = the non-standard title for the legend
  # xlim = a vector of limits for the times
  
  exArgs=list(...)

  # Creates the dataset to plot with the survival curves for all profiles  
  strata=lapply(1:nrow(psa$des.mat),function(x) {
    psa$des.mat %>% as_tibble() %>% select(-matches("(Intercept)",everything())) %>% slice(x) %>% 
      round(digits=2) %>% mutate(strata=paste0(names(.),"=",.,collapse=","))
  }) %>% bind_rows(.) %>% select(strata)
  toplot=lapply(1:length(psa$S),function(i) {
    psa$S[[i]] %>% bind_cols(strata=as.factor(strata[i,]))
  }) %>% bind_rows(.)

  if(exists("alpha",where=exArgs)){alpha=exArgs$alpha} else {alpha=0.2}
  if(exists("name_labs",where=exArgs)){name_labs=exArgs$name_labs} else {name_labs="Profile"}
  
  psa.plot=ggplot(data=toplot,aes(x=t,y=S,colour=strata))+
    geom_line(size=.9) +
    theme_bw() + 
    theme(axis.text.x = element_text(color="black",size=12,angle=0,hjust=.5,vjust=.5),
          axis.text.y = element_text(color="black",size=12,angle=0,hjust=.5,vjust=.5),
          axis.title.x = element_text(color="black",size=14,angle=0,hjust=.5,vjust=.5),
          axis.title.y = element_text(color="black",size=14,angle=90,hjust=.5,vjust=.5)) +
    theme(axis.line = element_line(colour = "black"),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size=18, face="bold")) +
    theme(legend.position=c(.75,.9),
          legend.title=element_text(size=15,face="bold"),
          #legend.title = element_blank(),
          legend.text = element_text(colour="black", size=14, face="plain"),
          legend.background=element_blank()) +
    labs(y="Survival",x="Time",title=NULL,
         color=name_labs) 

    # If there are more than 1 simulation, then there also are the low and upp extremes and plots them too
  if(any(grepl("low",names(toplot)))) {
    psa.plot=psa.plot+geom_ribbon(data=toplot,aes(x=t,y=S,ymin=low,ymax=upp,fill=strata),alpha=alpha,show.legend=F)
  }
  
  # Optional arguments
  if(exists("col",where=exArgs)) {
    psa.plot=psa.plot+scale_color_manual(values=exArgs$col)+scale_fill_manual(values=exArgs$col)
  }
  if(exists("xlab",where=exArgs)){
    psa.plot=psa.plot+labs(x=exArgs$xlab)
  }
  if(exists("ylab",where=exArgs)){
    psa.plot=psa.plot+labs(y=exArgs$ylab)
  }
  if(exists("main",where=exArgs)) {
    psa.plot=psa.plot+labs(title=exArgs$main)+theme(plot.title=element_text(size=18,face="bold"))
  }
  if(exists("labs",where=exArgs)) {
    psa.plot=psa.plot+scale_color_discrete(labels=exArgs$labs)
  }
  if(exists("xlim",where=exArgs)){
    psa.plot=psa.plot+xlim(exArgs$xlim)
  }
  
  # Finally renders the plot
  psa.plot
}
