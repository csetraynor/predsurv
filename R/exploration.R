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
plot_tte_dist <- function(d, time = os_months, status = os_status, alpha = 0.5){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  theme_set(theme_bw())
  p <- d %>%
    ggplot(aes_string(x = rlang::quo_text(time),
                      group = rlang::quo_text(status),
                      colour = rlang::quo_text(status),
                      fill =  rlang::quo_text(status))) +
    geom_density( alpha = alpha)+
    scale_color_discrete(name = "Overall Survival", label=c("Deceased", "Living")) +
    scale_fill_discrete(name = "Overall Survival", label=c("Deceased", "Living"))
  p + labs(x = "Time (months)", title= "Time to event density distribution")
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
plot_km <- function(d, time = os_months, status = os_status, event_type = "DECEASED", x= "survival"){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  mle.surv <- survival::survfit(survival::Surv( time2event , os_event) ~ 1,
                      data = d %>%
                        dplyr::mutate(os_event = ( !!status == event_type),
                                      time2event = !!time) )

  ggplot2::autoplot(mle.surv, conf.int = F) +
    ggtitle(paste0('KM for', x ,' Cohort'))
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
