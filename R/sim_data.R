#' Simulate data, function 1 with Weibull distribution.
#'
#' This function simulates a survival dataset given the parameters and covariates.
#'
#' @param
#' n number of individuals default 80\cr
#' features number of features in the experiment default 1000 \cr
#' @return simulated data
#' @export
#' @importFrom magrittr %>%
surv_sim_data <- function( alpha = 0.8, mu = -3, n = 80, features = 1000, seed = 111) {
  set.seed(seed)
  x<-matrix(rnorm(features*n),ncol= n)

  survdata <- data.frame(surv_months = rweibull(n = n, alpha, exp(-(mu + svd(x[1:60,] )$v[,1] + + .1*rnorm(n))/alpha)),
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
  featurenames <- paste("feature",as.character(1:1000),sep="")
  x <- t(x)
  colnames(x) <- featurenames

  survdata <- cbind(survdata, x)

  return(survdata)
}


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
