


#' Extract particular features from cv output
#'
#' This function should not be called directly
#'
#' @param
#' data : cv_obj should be a list \cr
#' feature: brier_pred, roc_pred...
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @importFrom rlang !!!
#' @import prodlim

extract_pars_cv <- function(cv_obj, feature, model, time = os_months, status = os_deceased,  ...){

  outList <- lapply(seq_along(cv_obj), function(i) attr(cv_obj[[i]], feature))
  if(feature == 'brier_pred'){
  out.data <- plyr::ldply(outList, function(i)
    out.data <- data.frame(
        ibrier = pec::crps(i, models = "matrix")[1],
        model = model
    ))}
  if(feature == 'roc_pred'){
    out.data <- plyr::ldply(outList, function(i)
      out.data <- data.frame(
        ROC = i$AUC[1],
        model = model
      ))}
  if(feature == 'ci_pred'){
    out.data <- plyr::ldply(outList, function(i)
      out.data <- data.frame(
        "c-index" = i$AppCindex$coxph,
        model = model
      ))}
  if(feature == 'dev_pred'){
    out.data <- plyr::ldply(outList, function(i)
      out.data <- data.frame(
        deviance = i[2],
        model = model
      ))}
  return(out.data)
}

