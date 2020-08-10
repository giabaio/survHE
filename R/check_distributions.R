#' Helper function to check that the distribution(s) provided by the user are
#' consistent with the method chosen for inference.
#' 
#' \code{'inla'} or \code{'hmc'}). If \code{method} is set to \code{'hmc'},
#' then \code{survHE} will write suitable model code in the Stan language
#' (according to the specified distribution), prepare data and initial values
#' and then run the model.
#' @param distr3 A vector of distribution labels (as created by 'fit.models'). 
#' It's a 3-letters label to identify the distributions
#' @param availables A list with the distributions available for each method.
#' @note Something will go here
#' @author Gianluca Baio
#' @seealso fit.models
#' @references Baio (2020). survHE
#' @keywords Parametric survival models 
check_distributions <- function(method,distr) {
  # Loads in the available models in each method
  availables <- load_availables()
  # Uses the helper 'manipulated_distributions' to create the vectors distr, distr3 and labs
  distr3 <- manipulate_distributions(distr)$distr3
  
  # If 'method' is either 'inla' or 'hmc but we're trying to run a model that is not available, then
  # falls back to 'mle'
  if(method %in% c("inla","hmc")) {
    if(!all(distr3 %in% availables[[method]])) {
      modelsString <- unname(labelTable[availables[[method]]])
      modelsString[length(modelsString)] = paste0("or ", modelsString[length(modelsString)])
      message(paste0(
        "NB: ",toupper(method)," can only fit ",
        paste(modelsString, collapse = ", "),
        " parametric survival models. Falling back on MLE analysis")
      )
    }
    method <- "mle"
  }
  
  # 'mle' can implement all the possible models, excpet the PolyWeibull
  # In this case, I choose to *stop* execution, rather than falling back to 'hmc'!
  if (method == "mle") {
    if(!all(distr3 %in% availables[[method]])) {
      stop(paste0("The Poly-Weibull model is only implemented under method='hmc'.
       Please set this option in your call to 'fit.models'"))
    }
  }
}