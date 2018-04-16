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
surv_sim_data <- function(  mu = -3, alpha = 5,  N = 80, features = 1000, seed = 111, p = 0.25,  cens_rate = 1/10) {
  set.seed(seed)


  x<-matrix(rnorm(features*N),ncol= N)
  beta <- rnorm(round(features * p, 0 ))

  survdata <- data.frame(surv_months = rweibull(n = N, alpha, exp(-(mu +  t(x)[, sample(1:features, size = length(beta) )] %*% beta )/alpha)),
                         censor_months = rexp(n = N, rate = cens_rate),
                         stringsAsFactors = F
  ) %>%
    dplyr::mutate(os_status = ifelse(surv_months < censor_months,
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months,
                       surv_months, censor_months
    ),
    os_deceased = (os_status == "DECEASED")
    ) %>%
    dplyr::select(os_months, os_status, os_deceased)
  featurenames <- paste("feature",as.character(1:features),sep="")
  x <- t(x)
  colnames(x) <- featurenames

  survdata <- cbind(survdata, x)

  survdata <- survdata %>% dplyr::select(-os_status)

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

  data <- data.frame(surv_months = rweibull(n = N, alpha, exp(-(mu + svd(x[sample(1:N, size = round(N * p, 0 )),])$v[,1]  +  .1*rnorm(N))/alpha)),
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
