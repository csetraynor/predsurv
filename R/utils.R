#' Load Data to environment
#'
#' @param
#' rdata \cr
#' @return a coxph fit object
#' @export
LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}





