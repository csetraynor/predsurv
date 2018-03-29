#' Fit MAP estimates for the Poisson Model PEM.   Important aspects of the modelling #'approach: Dependent variable is the binary event indicator. Adding an intercept is #'recommendable to provide a decent and interpretable reference category for factor #'variables. Use the logarithm of the interval length as an offset.
#'
#'
#' @param
#' d a dataset \cr
#' formula  \cr
#' a log baseline hazard
#' t_dur duration time for individtual ith in the hth interval
#' lambda hazard rate
#' os_event censoring indicator (actual event) for the ith subject in the hth interval
#' @return a MAP fit
#' @export
#' @importFrom rethinking map
#' @importFrom rlang !!
map_pem <- function(l , formula = tobacco_smoking_history_indicator){

  m1 <- rethinking::map(
    alist(
      status ~ dpois(lambda),
      lambda <- a[t_id] + offset(log(t_dur)),
      a ~ dgamma(2, .1)
    ), data=l,  start=list(a=0, lambda = 0.1) )

}
