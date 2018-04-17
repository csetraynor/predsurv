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
surv_sim_data <- function(  N = 80, features = 100, seed = 111, p = 0.25,  rate = 1/100) {
  set.seed(seed)
  x<-matrix(rnorm(features*N, mean = 0, sd = 1),ncol=  features)
  beta <- rnorm(round(features * p, 0 ), mean = 0, sd = 10)
  #actual x
  x_true <- x[,sample(1:ncol(x), length(beta) )]
  alpha = 1 + runif(1)
  outdata <- data.frame(surv_months = rweibull(n = N,
                 alpha,exp(-rate *(  x_true %*% beta )/alpha)),
                         censor_months = rexp(n = N, rate = rate),
                         stringsAsFactors = F
  ) %>%
    dplyr::mutate(os_status = ifelse(surv_months < censor_months,
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months,
                       surv_months, censor_months
    ),
    os_deceased = (os_status == "DECEASED")
    )
  outdata <- outdata %>% dplyr::select(os_months, os_status, os_deceased)
  featurenames <- paste("feature",as.character(1:features),sep="")
  colnames(x) <- featurenames

  outdata <- cbind(outdata, x)

  outdata <- outdata %>% dplyr::select(-os_status) %>% dplyr::mutate(os_months = os_months*12)

  return(outdata)
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
