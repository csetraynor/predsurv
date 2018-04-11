#' Calculate Integrated Brier Score
#'
#' This function recieves a training set and a test set and calculates
#' different measures of accuracy and loss functions
#' @param
#' pred_data a dataset prediction made by the training set (or CV fold)  \cr
#' test_data the hold out dataset \cr
#' cutoff an optional parameter to integrate the Brier score in a different time, default is max observed time \cr
#' @return pp_bried the posterior predicted integrated brier score
#' @export
brier_score <- function(pred_data, test_data, cutoff=NA){
  smod <- with(test_data %>% dplyr::mutate(
    os_deceased = os_status == 'DECEASED'),
    survival::Surv(os_months, os_deceased))

  pred_KM <- lapply(pred_data %>%
                      purrr::map(~ dplyr::mutate(.,
                                   os_deceased = os_status == 'DECEASED')),
                    function(x){
                      survival::survfit(Surv(os_months, os_deceased) ~ 1, data = x)
                    })

  # integrated Brier score up to max time
  pp_brier <- purrr::map(.x = pred_KM,
                         .f = ~ ipred::sbrier(obj = smod,
                                             pred = .x))
  # ibs[[i]] <- pp_brier %>% unlist()
  return(pp_brier)
}

#' Calculate Time dependent ROC Curve
#' There has been a lot of discussing about the usefullness of ROC curves\cr
#' The time dependent ROC Curves can be used in Survival analysis but only\cr
#' make sense if all the models to be evaluated are based on a similar link\cr
#' function. For example all models are Cox Models with different priors,\cr
#' penalised likelihoods or parameters. However they don't make sense to compare\cr
#' a tree based method with a Cox Model, or a Weibull against a Cox model \cr
#' because they structure.
#' @param
#' pred_newdata a dataset prediction made by the training set (or CV fold)  \cr
#' tet_data the hold out dataset \cr
#' t an optional parameter to integrate the Brier score in a different time, default is max observed time \cr
#' @return pp_bried the posterior predicted integrated brier score
#' @export
roc_curve <- function(newdata, time = os_months, status=os_status,
                      npi = npi, cutoff= 12, plot = F){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  nobs <- NROW(newdata)
  npi_fit <- survivalROC::survivalROC(Stime = newdata %>% dplyr::select(!!time) %>% unlist,
                         status = nedata %>% dplyr::select(!!status) %>% unlist,
                         marker = nedata %>% dplyr::select(!!npi) %>% unlist,
                         predict.time = cutoff, span = 0.25*nobs^(-0.2))
  if(plot == F){
    return(npi_fit)
  }else{
    plot(npi_fit$FP, npi_fit$TP, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab=paste( "FP", "\n", "AUC = ",round(npi_fit$AUC,3)),
         ylab="TP",main=paste0("NPI, Method = NNE \n Months = ", cutoff))
    abline(0,1)
  }
}
