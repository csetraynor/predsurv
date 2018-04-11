#' Simulate data, function 1 with Weibull distribution.
#'
#' This function simulates a survival dataset given the parameters and covariates.
#' 
#' @param
#' rate rate parameter for the exponential distribution \cr
#' beta beta regression coefficient\cr
#' @return simulated data
#' @export
#' @importFrom magrittr %>%
sim_data <- function(mu, alpha, beta) {
  lambda = mu +  X %*% beta ;
  
  data <- data.frame(surv_months = rweibull(n = n, alpha, exp(-(mu + pi )/alpha)),
                     censor_months = rexp(n = n, rate = 1/100),
                     stringsAsFactors = F
  ) %>%
    dplyr::mutate(os_status = ifelse(surv_months < censor_months,
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months,
                       surv_months, censor_months
    )
    )
  
  return(data)
}