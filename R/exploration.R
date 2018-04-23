#' Time to Event distribution
#'
#' This function plots the TTE distribution
#'
#' @param
#' d a dataset \cr
#' ... vars to be evaluated \cr
#' n optional parameter to assure number of rows deleted, default NA
#' @return d a clean dataset
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
plot_tte_dist <- function(d, time = os_months, status = os_status, alpha = 0.5){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  theme_set(theme_bw())
  p <- d %>% dplyr::mutate(status = as.factor( !!status ),
                           time = !!time) %>%
    ggplot(aes(x = time,
                      group = status,
                      colour = status,
                      fill =  status)) +
    geom_density( alpha = alpha)+
    scale_color_discrete(name = "Overall Survival", label=c("Deceased", "Living")) +
    scale_fill_discrete(name = "Overall Survival", label=c("Deceased", "Living"))
  p + labs(x = "Time (months)", title= "Time to event density")
}

#' Kaplan and Meier estimates
#'
#' This function plots the TKaplan and Meier estimates
#'
#' @param
#' d a dataset \cr
#' x customise plot title \cr
#' event_type codification of event type , default DECEASED
#' @return d a clean dataset
#' @export
#' @importFrom rlang !!
#' @import ggfortify
#' @importFrom magrittr %>%
plot_km <- function(d, time = time, status = status, event_type = 1, x= "survival"){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  mle.surv <- survival::survfit(survival::Surv( time2event , os_event) ~ 1,
                      data = d %>%
                        dplyr::mutate(os_event = ( !!status == event_type),
                                      time2event = !!time) )

  ggplot2::autoplot(mle.surv, conf.int = F) +
    labs(x = "Time (months)", y = "Survival probability", title= "K-M")
}


#' Fit Exponential distribution
#'
#' This function fits the Weibull distribution
#'
#' @param
#' d a dataset \cr
#' x customise plot title \cr
#' event_type codification of event type , default DECEASED
#' @return d a clean dataset
#' @export
#' @importFrom rlang !!
#' @import ggfortify
#' @importFrom magrittr %>%
fit_weibull <- function(d){
  expfit <-   survival::survreg(survival::Surv(os_months, os_deceased) ~ age, data = d %>% mutate(os_deceased = (os_status == "DECEASED")) %>% filter(os_months > 0), dist = "exp")
  scale = exp(expfit$coefficient[1])
  beta = weifit$coefficients[2]
}



#' Fit Weibull distribution
#'
#' This function fits the Weibull distribution
#'
#' @param
#' d a dataset \cr
#' x customise plot title \cr
#' event_type codification of event type , default DECEASED
#' @return d a clean dataset
#' @export
#' @importFrom rlang !!
#' @import ggfortify
#' @importFrom magrittr %>%
fit_weibull <- function(d){
   weifit <-   survival::survreg(survival::Surv(os_months, os_deceased) ~ age, data = d %>% mutate(os_deceased = (os_status == "DECEASED")) %>% filter(os_months > 0))
  shape  =   1/weifit$scale
  scale = exp(weifit$coefficient[1])
  beta = weifit$coefficients[2]
}

#' The function for plot correlation between explanatory variables
#'
#' @param
#' d a dataset \cr
#' x customise plot title \cr
#' event_type codification of event type , default DECEASED
#' @return d a clean dataset
#' @export
#' @importFrom rlang !!
plot_correlation <- function(d, time = os_months, status = os_status, event_type = "DECEASED", x= "survival"){

}
