#' This function generates the long data neeeded for a piecewise exponential model
#'
#' @param
#' d a dataset \cr
#' formula a formula to pass covariates \cr
#' k the number of intervals that want to be done
#' c an optional parameter to specify lenghts of intervals
#' @return a long format (count model) dataset
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
gen_long_data <- function(d, time = os_months, status = os_status, sample_id = patient_id, k, c=1, event_type = "DECEASED") {

  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  sample_id <- dplyr::enquo(sample_id)

  #add small error to time 0
  add_epsilon <- function(x) {
      ifelse(x == 0, x + 1e-3, x)
  }
  d <- d %>% dplyr::mutate_at(rlang::quo_text(time), dplyr::funs(add_epsilon))

  #set the tau interval times
  tau <- d %>% dplyr::filter(!!status == event_type) %>%
    dplyr::select(!!time) %>% unlist %>% unique %>% sort()
  tau <- tau[seq(1, k, c)]
  print(tau)

  longdata <- survival::survSplit(Surv(time2event, os_event ) ~.,
                        data = d %>%
                        dplyr::mutate(os_event = ( !!status == event_type),
                                      time2event = !!time) ,
                        cut = tau)
  #Create time ID and calculate t_dur
  longdata <- longdata %>%
    dplyr::group_by(!!sample_id) %>%
    dplyr::mutate(t_id = seq(n())) %>%
    dplyr::mutate(t_dur = time2event - tstart) %>%
    dplyr::ungroup()


  return(longdata)
}

#' This function generates the stan data for a clinical PEM model, the rason for not #' introducing genomic covariates is that it takes too much memory, for instance the #' genomic covariates will be handled ina module in Python.
#'
#' @param
#' d a dataset \cr
#' formula a formula to pass covariates \cr
#' time = time to event
#' status = status covariate
#' sample_id = unique id for each patient
#' time_id = unique id for each observation
#' @return a list for stan
#' @export
#' @importFrom rlang !!
#' @importFrom magrittr %>%
#' @import  dplyr
gen_stan_data <- function(d, formula = as.formula(~1), time = t_dur, status = os_event, sample_id = patient_id , time_id = time_id) {
  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  X <- d %>%
    stats::model.matrix(formula, data = .)
  P <- ncol(X)
  if (P > 1){
    if("(Intercept)" %in% colnames(X))
      X <- array(X[,-1], dim = c(nrow(d), P - 1))
    P <- ncol(X)
  }

  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  sample_id <- dplyr::enquo(sample_id)
  time_id <- dplyr::enquo(time_id)

  patient_id <- factor(unique(d %>% dplyr::select(!!sample_id) %>% unlist),
                levels = unique(d %>% dplyr::select(!!sample_id) %>% unlist ))

  stan_data <- list(
    N = nrow(d),
    ID = length(unique(d %>% dplyr::select(!!sample_id))),
    "T" = dplyr::n_distinct(d %>% dplyr::select(!!time)),
    s = as.integer(patient_id),
    t_dur = d %>% dplyr::select(!!time) %>% unlist,
    status = d %>% dplyr::select(!!status) %>% unlist,
    t_id = d %>% dplyr::select(!!time_id) %>% unlist,
    X = X,
    P = P
  )
  return(stan_data)
}







  # stan_data <- list(
# t_obs <- data %>% select(os_months) %>% unlist %>% unique %>% sort()
# t_dur <- diff(tau)
# X_bg <- longdata %>%
#   model.matrix(formula, data = .)
# M_bg <- ncol(X_bg)
# if (M_bg > 1){
#   if("(Intercept)" %in% colnames(X_bg))
#     X_bg <- array(X_bg[,-1], dim = c(nrow(longdata), M_bg - 1))
#   M_bg <- ncol(X_bg)
# }
  #   N = nrow(longdata),
  #   S = length(unique(longdata$num_id)),
  #   "T" = dplyr::n_distinct(longdata$t_id),
  #   s = as.integer(longdata$num_id),
  #   log_t_dur = array(log(t_dur), dim = length(t_dur)),
  #   M = M_bg,
  #   status = as.integer(longdata$deceased),
  #   t = longdata$t_id,
  #   x = X_bg
  # )



#
#
#
# #---Set initial values---#
# gen_inits <- function(M, data) {
#   tau <- data %>% select(os_months) %>% unlist %>% unique %>% sort()
#   if(tau[1] != 0){
#     tau <- c(0, tau)
#   }
#   t_dur <- diff(tau)
#   function()
#     list(
#       log_baseline_raw = rnorm(n= length(t_dur)),
#       sigma_baseline = runif(1, 0.01, 100),
#       log_baseline_mu = rnorm(1),
#       beta_bg_raw = rnorm(M),
#       tau_s_bg_raw = 0.1*abs(rnorm(1)),
#       tau_bg_raw = abs(rnorm(M))
#     )
# }
