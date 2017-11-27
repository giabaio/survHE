#' Create an individual level dataset from digitised data
#' 
#' Piles in the simulated IPD resulting from running digitise for more than one
#' treatment arm
#' 
#' 
#' @param ipd_files a list including the names of the IPD files created as
#' output of digitise
#' @param ctr the index of the file associated with the control arm (default,
#' the first file). This will be coded as 0
#' @param var.labs a vector of labels for the column of the resulting data
#' matrix. NB these should match the arguments to the formula specified for
#' fit.models. The user can specify values. These should be 4 elements (ID,
#' TIME, EVENT, ARM)
#' @author Gianluca Baio
#' @seealso Something will go here
#' @references Something will go here
#' @keywords Digitized Kaplan Meier curve
#' @examples
#' 
#' # Something will go here
#' 
#' @export make.ipd
make.ipd <- function(ipd_files,ctr=1,var.labs=c("time","event","arm")) {
  ## Piles in the simulated IPD resulting from running digitise for more than one treatment arm  
  ## ipd_files = a list including the names of the IPD files created as output of digitise
  ## ctr = the index of the file associated with the control arm (default, the first file).
  ##       This will be coded as 0
  ## var.labs = a vector of labels for the column of the resulting data matrix. NB these
  ##            should match the arguments to the formula specified for fit.models. The
  ##            user can specify values. These should be 3 elements (TIME, EVENT, ARM)
  
  # Identifies the number of arms (= number of IPD files)
  n_arms <- length(ipd_files)
  index <- 1:n_arms
  active <- index[-ctr]
  data <- read.table(ipd_files[[ctr]],header=TRUE,row.names=NULL)
  data[,"arm"] <- 0 # sets the value of "arm" to 0, for the control group
  arm.ind <- 1
  for (i in active) {
    tmp <- read.table(ipd_files[[index[i]]],header=TRUE,row.names=NULL)
    tmp[,"arm"] <- arm.ind
    data <- rbind(data,tmp)
    arm.ind <- arm.ind+1
  }
  colnames(data) <- var.labs
  return(data)
}
