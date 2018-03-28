#' Fit a Poisson distribution
#'
#' This function plots the TTE distribution,   Important aspects of the modelling approach: Dependent variable is the binary event indicator.  Adding an intercept is recommendable to provide a decent and interpretable reference category for factor variables.
#Use the logarithm of the interval length as an offset.  Add a smoothing-spline to model the baseline hazard.
#'
#' @param
#' d a dataset \cr
#' formula covariate matrix \cr
#' @return a fit objects ( smoth splines).
#' @export
#' @importFrom magrittr %>%
fit_poisson <- function(d , t_start = tstart, t_stop = time2event, t_dur = t_dur, status = os_event, formula = tobacco_smoking_history_indicator){

  t_start = dplyr::enquo(t_start)
  t_stop = dplyr::enquo(t_stop)
  t_dur = dplyr::enquo(t_dur)
  status = dplyr::enquo(status)
  formula = dplyr::enquo(formula)
  ######Poisson trick

  #Important aspects of the modelling approach: * Dependent variable is the binary event indicator. * Adding an intercept is recommendable to provide a decent and interpretable reference category for factor variables.
  #Use the logarithm of the interval length as an offset. * Add a smoothing-spline to model the baseline hazard.

  fit <- mgcv::gam(os_event~1+offset(log(t_dur))+s(t_stop)+factor(formula),
             d %>% dplyr::mutate(
               t_dur = !!t_dur,
               t_stop = !!t_stop,
               formula = !!formula,
               os_event = !!status
             ),family='poisson')
  fit2 <- mgcv::gam(os_event~1+offset(log(t_dur))+s(t_start)+factor(tobacco_smoking_history_indicator),
              d %>% dplyr::mutate(
                t_dur = !!t_dur,
                t_start = !!t_start,
                formula = !!formula,
                os_event = !!status),
              family='poisson')
  output <- list(fit, fit2)
  return(output)
}


#' Make predictions for the PEM. Series of functions in the PEM MLE.
#'
#'
#' @param
#' d a dataset \cr
#' @return a dataframe ready to be plotted
#' @export
#' @importFrom magrittr %>%
pred_poisson <- function(d , time = os_months){

  time = dplyr::enquo(time)

  time <- c(sort(unique(d %>% dplyr::select(!!time) %>% unlist)))
  t_dur <- diff(c(0, time))
  p1<-stats::predict(mgcv::gam(os_event~1+offset(log(t_dur))+s(time2event),d,family='poisson'),data.frame(t_dur=t_dur,time2event=time))
  S1<-exp(-cumsum(exp(p1)))
  plotframe<-data.frame(Time=c(0,time),S1=c(1,S1))
  return(plotframe)
}
#' Prediction plot for the PEM
#'
#' This function plots the actual Kaplan and Meier and overlies the PEM prediction, follows the previous predict function
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
plot_pred_poisson <- function(d, time = os_months, status = os_status, event_type = "DECEASED", plotdata = plotframe, predtime = Time, predsurv = S1 , formula = as.formula(~1)){

  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  form <- as.formula(paste("Surv(time2event, os_event)", " ~ ", rlang::get_expr(formula)))

  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  predtime <- dplyr::enquo(predtime)
  predsurv <- dplyr::enquo(predsurv)

  plot(survfit(form,
                    data = d %>%
                       dplyr::mutate(os_event = ( !!status == event_type),
                            time2event = !!time)  ),
       ylab='Survival Probability',xlab='Time',
       main='Baseline Survival Curve')
  lines(plotframe %>% dplyr::select(!!predtime) %>% unlist,
        plotframe %>% dplyr::select(!!predsurv) %>% unlist,
        col='green')
  return(form)
}

