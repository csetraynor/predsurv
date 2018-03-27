#' Calculate Penalised Likelihood
#'
#' This function calculates the penalised likelihood of the Cox model
#' @param
#' train datset for training  \cr
#' formula \cr
#' cutoff an optional parameter to integrate the Brier score in a different time, default is max observed time \cr
#' @return pp_bried the posterior predicted integrated brier score
#' @export
penalised_likelihood_fit <- function(formula=as.formula(~1) , y, family = "cox", plot=FALSE){
  fit = glmnet::glmnet(x, y, family = family)
  if(plot == FALSE){
    return(fit)
  }else{
    plot(fit)
  }
}
