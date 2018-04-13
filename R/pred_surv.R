#' Cox Model Prediction
#'
#' This function takes a survival fit and following the Cox model, estimates the hazard for individual from the model coefficient, while the baseline hazard is estimated with the Nelson-Aalen method.
#'
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
pred_surv <- function(obj, time = os_months, status = os_deceased, beta = beta, test_data = test, train_data = train){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  # take optimal beta from model object
  optimal.beta <- beta(obj)
  # find non zero beta coef
  nonzero.coef <- abs(optimal.beta)>0 & !is.na(optimal.beta)
  selectedBeta <- optimal.beta[nonzero.coef]

  # take only covariates for which beta is not zero
  trainX <- train_data %>% dplyr::select(-!!time, -!!status)
  selectedTrainX   <- trainX[,nonzero.coef]

  traincoxphdata <- cbind(train_data %>% dplyr::select(!!time, !!status), selectedTrainX)

  #create formula
  testBeta <- paste(names(selectedBeta), collapse = " + ")

  survformula <- as.formula(paste("Surv(time, status)", " ~ ", testBeta))

  # create coxph object with pre-defined coefficients
  coxph.model<- coxph(survformula ,init=selectedBeta,iter=0, data = traincoxphdata %>% dplyr::mutate(
    time = !!time,
    status = !!status))

  # take test covariates for which beta is not zero
  testX <- test_data %>% dplyr::select(-!!time, -!!status)
  selectedTestX <- testX[,nonzero.coef]
  ndata <- cbind(test_data %>% dplyr::select(!!time, !!status), selectedTestX)


  pred_coxph <- predict(coxph.model,
                       newdata = ndata %>%
                         dplyr::mutate(
                         time = !!time,
                         status = !!status),
                       "expected")

  surv_pred <- data.frame(time = test %>% select(!!time) %>% unlist,
                          status = round(exp(-pred_coxph), 0))


  return(surv_pred)
}
