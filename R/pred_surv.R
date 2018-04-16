#' Fit Forward Selection
#'
#' Fit the Forward Selection method
#'
#' @param
#' train : train data \cr
#' time : time variable \cr
#' status : status variable \cr
#' method: forward, backward or both \cr
#' fit : fit model: "stepwise", "random forest"
#' @return a coxph forward stepwise fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!

fun_fit <- function(train, time = os_months, status = os_deceased, method = "forward", fit, penalty = "BIC" ){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)


  # take only covariates for which beta is not zero
  trainX <- train %>% dplyr::select(-!!time, -!!status)

  traincoxphdata <- cbind(train %>% dplyr::select(!!time, !!status), trainX)

  #create formula


   form <- as.formula(paste("Surv(time, status)", paste("~", paste(names(trainX), collapse = " + "), sep = " " )))

  if(fit == "stepwise"){
    # create coxph object
    fit1 <- suppressWarnings(survival::coxph(Surv(time, status)~1,
                            data = traincoxphdata %>% dplyr::mutate(
                              time   = !!time,
                              status = !!status) ))

    testBeta <- as.formula(paste("~", paste(names(trainX), collapse = " + "), sep = " " ) )

    if(penalty == "AIC"){
      k = 2
    }else{
      if(penalty == "BIC"){
        k = log(length(names(trainX)))
      }
    }

    #forward stepwise selection
    mod <- suppressWarnings(MASS::stepAIC(fit1,
            scope = list(upper = testBeta, lower = ~1) ,
            trace = FALSE,
            direction = method,
            k = k,
            steps =  nrow(train)-2 ))
  }

  if(fit == "random forest"){
    mod <- party::cforest(form, data = traincoxphdata %>%
                            dplyr::mutate(
                              time   = !!time,
                              status = !!status),
                          control = cforest_unbiased(ntree = 50))
  }

  return(mod)
}
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
#' @import prodlim
mod_pred <- function(obj, time = os_months, status = os_deceased, beta = beta, test_data = test, train_data = train, fit = "stepwise"){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  # take optimal beta from model object
  optimal.beta <- beta(obj)
  # find non zero beta coef
  nonzero.coef <- abs(optimal.beta)>0 & !is.na(optimal.beta)
  selectedBeta <- optimal.beta[nonzero.coef]

  # take only covariates for which beta is not zero
  trainX <- train_data %>% dplyr::select(-!!time, -!!status)
  selectedTrainX   <- trainX[,colnames(trainX) %in% names(nonzero.coef)]

  traincoxphdata <- cbind(train_data %>% dplyr::select(!!time, !!status), selectedTrainX)

  #create formula
  testBeta <- paste(names(selectedBeta), collapse = " + ")

  form <- as.formula(paste("Surv(time, status)", " ~ ", testBeta))

  # create coxph object with pre-defined coefficients
  if(fit == "stepwise"){
    mod <-  coxph(form , data = traincoxphdata %>% dplyr::mutate(
      time = !!time,
      status = !!status))
  }else{
    mod <- coxph(form , data = traincoxphdata %>% dplyr::mutate(
      time = !!time,
      status = !!status), init=selectedBeta, iter=0)
  }

  # take test covariates for which beta is not zero
  testX <- test_data %>% dplyr::select(-!!time, -!!status)
  selectedTestX <- testX[,colnames(testX) %in% names(nonzero.coef)]
  ndata <- cbind(test_data %>% dplyr::select(!!time, !!status), selectedTestX)
  timepoints <-  ndata %>% dplyr::select(!!time) %>%  unlist

  probs = pec::predictSurvProb(mod, newdata = ndata, times = timepoints)
  print(probs)

  return(probs)
}


#' Brier Score
#'
#' Predict integration Brier Score in a new dataset.
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
fun_brier_score <- function(obj = mod.fs.coxph, pred=pred.fs.coxph, time = os_months, status = os_deceased,  test_data = test, beta = beta){

  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  # take optimal beta from model object
  optimal.beta <- beta(obj)
  # find non zero beta coef
  nonzero.coef <- abs(optimal.beta)>0 & !is.na(optimal.beta)

  # take test covariates for which beta is not zero
  testX <- test_data %>% dplyr::select(-!!time, -!!status)
  selectedTestX <- testX[,nonzero.coef]
  ndata <- cbind(test_data %>% dplyr::select(!!time, !!status), selectedTestX)
  timepoints <-  ndata %>% dplyr::select(!!time) %>% unlist

  brier <- pec::pec(pred, Surv(time, status) ~ 1, data = ndata %>% dplyr::mutate(
    time   = !!time,
    status = !!status), maxtime = max(timepoints), times = timepoints)


  return(brier)

}


#' Integrated Brier Score
#'
#' Predict integration Brier Score in a new dataset.
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
fun_ibrier_score <- function(brier){

  ibrier <- pec::crps(brier, models = "matrix")

  return(ibrier)

}


#' Plot Random Forest
#'
#' Plot Fitted Individual Tree
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!

random_forest_plot <- function(bst){

  ### if you can't resist to look at individual trees ...
  party:::prettytree(bst@ensemble[[1]], names(bst@data@get("input")))

}




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
#' @import prodlim
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

  surv.formula <- as.formula(paste("Surv(time, status)", " ~ ", testBeta))

  # create coxph object with pre-defined coefficients
  coxph.model<- coxph(surv.formula ,init=selectedBeta,iter=0, data = traincoxphdata %>% dplyr::mutate(
    time = !!time,
    status = !!status))

  # take test covariates for which beta is not zero
  testX <- test_data %>% dplyr::select(-!!time, -!!status)
  selectedTestX <- testX[,nonzero.coef]
  ndata <- cbind(test_data %>% dplyr::select(!!time, !!status), selectedTestX)

  pred_ibrier = pec::pec(list(coxph.model),
                         data = ndata %>%
                           dplyr::mutate(
                             time = !!time,
                             status = !!status),
                         formula=Surv(time,status)~1)


  #
  #
  # pred_coxph <- predict(coxph.model,
  #                      newdata = ndata %>%
  #                        dplyr::mutate(
  #                        time = !!time,
  #                        status = !!status))
  #
  #   surv_pred <- data.frame(time = test %>% select(!!time) %>% unlist,
  #                           surv_prob = exp(-pred_coxph))


  return(pred_ibrier)
}

