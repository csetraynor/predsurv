#' Relaxed Cox Lasso
#'
#' This function performs the so called relaxed Lasso for a Cox model.
#'
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
#' @import survival
rel_shrink <- function(obj, data, time = os_months, status = os_deceased, subject = patient_id, table_bool = TRUE, ...){

  ###### abstract dependent vars
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  subject <- dplyr::enquo(subject)

  ##### convert to dataset
  train <- as.data.frame(data$data)[-data$in_id,] %>%
     dplyr::select(-!!subject)
  # find lambda for which dev.ratio is max
  selectedBeta <- obj$Predictor
  inits <- obj$Coefficient

  # Create train X take only covariates for which beta is not zero
  X <- train %>% dplyr::select(-!!time, -!!status)
  X   <- as.data.frame(X[,colnames(X) %in% selectedBeta])
  # prepare for coxph model

  fit <-  survival::coxph(survival::Surv(time = train %>%
                                           dplyr::select(!!time) %>%
                                           unlist,
                                         event = train %>%
                                           dplyr::select(!!status) %>%
                                           unlist)~.,
                          init = inits, iter = 20,
                          data = X)

  if(table_bool){
    out <- predsurv::par.table(fit)
  }else{
    out <- fit$coefficients
  }
  return(out)
}

#' Function par table
#'
#' Extract features of a classical fit model
#' @param
#' obj : survival coxph fit. \cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
#' @import survival
#' @author Sahota Tarj - caret

par.table <- function(fit){ ## needs glm object
  d1 <- summary(fit)$coefficients
  d1 <- as.data.frame(d1)
  dc <- as.data.frame(matrix(confint(fit),ncol=2))
  names(dc) <- c("lower","upper")
  d1 <- cbind(data.frame(Parameter=row.names(d1)),d1,dc)
  rownames(d1) <- NULL
  d1$description <- NA

  d1 <- d1 %>%
    dplyr::rename(Estimate = coef,
           "se_Estimate" = "se(coef)",
           HR = "exp(coef)"
           ) %>%
    dplyr::mutate(se = exp(se_Estimate),
                  lower = exp(lower),
                  upper = exp(upper)) %>%
    dplyr::select(Parameter, Estimate, HR, se_Estimate, se, dplyr::everything()) %>%
    dplyr::select(-Estimate,-se_Estimate)

  d1$description <-"Hazard ratio (relative SE)"

  return(d1)
}
