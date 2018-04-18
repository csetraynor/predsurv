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
surv_sim_data <- function(  N = 80, features = 100, seed = 111, p = 0.02, CenRate = 1/100, init = FALSE) {

  if(!init){
    X = matrix(rnorm(features*N, mean = 0, sd = 1),ncol=  features)
    beta = rep(0, features)
    #actual beta
    active <- round(features * p, 0 )
    beta_active <- rcauchy(active, 0,1)
    beta[1:active] <- beta_active
    beta <- sample(beta)
    #set alpha
    alpha = 1 + runif(1)
  }else{
    beta <- list(...)
  }

  outdata <- data.frame(surv_months = rweibull(n = N,
                 alpha,exp(- (X %*% beta) /alpha)),
                         censor_months = rexp(n = N, rate = CenRate),
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
  colnames(X) <- featurenames

  outdata <- cbind(outdata, X)

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

Sim_data <- function(N, features, p = 0.02, horizon = 48, init = FALSE, seed = 111, CenRate = 1/100, smooth = 10, ...){
  # set.seed(seed)
    #Initiate vars
  if(!init){
    tau = c(seq(0, horizon, length.out = round(N/smooth,0)+1 ))
    baseline <- invgamma::rinvgamma(length(tau)-1 , 2, 0.1)
    mean = rcauchy(features, 0, 2)
    X = sapply(mean, function(m) rnorm(N, mean = m))
    beta = rep(0, features)
    #actual beta
    active <- round(features * p, 0 )
    beta_active <- rcauchy(active, 0,2)
    beta[1:active] <- beta_active
    beta <- sample(beta)
  }else{
    beta <- list(...)
  }
  #format check
  beta <- as.vector(as.numeric(beta))
  baseline<- as.vector(as.numeric(baseline))
  X <- array(matrix(as.numeric(X)), dim = c(N, features ))
  #prognostic index
  mu = exp(X %*% beta)
  #extract first interval baseline hazard
  baseline0 <- baseline[1]
  #compute relative hazard for each interval respect to the first
  rel_base_risk <- baseline/baseline0
  rel_risk = lapply(mu, "*" , rel_base_risk)
  #caculate duration
  dt = diff(tau)
  assertthat::assert_that(length(dt) == length(baseline))
  #create a helping matrix for finding event time
  LD <- matrix(0, nrow = length(tau), ncol = length(rel_base_risk))
  LD[lower.tri(LD)] <- 1;
  #compute log survival
  logsurv <- log(1-runif(N))
  #compute log survival for each interval tau
  lsm = lapply(rel_risk, function(x) -baseline0 * as.vector(LD %*% (x*dt)))
  t <- (rep(NA,N))
  #find appropiate time interval
  t = mapply(function(x, y, z) {
    for (i in seq_along(baseline)) {
      t = ifelse(x[i]>=z & z>x[i+1], tau[i] + (x[i] - z)/baseline0/y[i], t)
    }
    return(t)
  } , x = lsm, y = rel_risk , z = as.list(logsurv)
  )
  #create data set
  sim.data <- data_frame(surv_months = t) %>%
    mutate(os_status = ifelse(is.na(surv_months), 'LIVING', 'DECEASED'),
           surv_months = ifelse(is.na(surv_months), tau[length(tau)], surv_months),
           id = seq(n),
           censor_months = rexp(n = N, rate = CenRate))   %>% #censoring rate
    dplyr::mutate(os_status = ifelse(surv_months < censor_months & os_status != 'LIVING',
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months  & os_status != 'LIVING',
                       surv_months, censor_months
    )
    ) %>%   cbind(X) #joint covariates
  return(sim.data)
}

#' Get high dimensional survival data.
#'
#' This function gets a survival data experiment
#'
#' @param
#' rate rate parameter for the exponential distribution \cr
#' beta beta regression coefficient\cr
#' @return simulated data
#' @export

get_data <- function(dataset = "survival"){
  data("beer.survival")
  data("beer.exprs")

  out <- cbind(beer.survival, t(as.matrix(beer.exprs)))
  colnames(out) <- gsub("-", "_", colnames(out))
  colnames(out) <- gsub("/", "_", colnames(out))
  rm(list = c("beer.exprs", "beer.survival"))
  return(out)
}


