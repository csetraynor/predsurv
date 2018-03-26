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
