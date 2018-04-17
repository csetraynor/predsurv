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
#' @import partykit
#' @import grid
#' @import libcoin
#' @import mvtnorm
#' @import rpart

fun_fit <- function(train, time = os_months, status = os_deceased, method = "forward", fit, penalty = "BIC" , tree_control = "Univariate"){
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

  if(fit == "tree"){
    ### with weight-dependent log-rank scores
    ### log-rank trafo for observations in this node only (= weights > 0)
    h <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
      if (is.null(weights)) weights <- rep(1, NROW(y))
      s <- coin::logrank_trafo(y[weights > 0,,drop = FALSE])
      r <- rep(0, length(weights))
      r[weights > 0] <- s
      list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
    }

    survtree_control <-
      partykit::ctree_control(testtype = tree_control, maxdepth = 2)
    mod <- partykit::ctree(form,
                      data = traincoxphdata %>%
                      dplyr::mutate(
                          time   = !!time, status = !!status),
                      ytrafo = h,
                      control = survtree_control )
  }



  if(fit == "random forest"){

    # survtree_control <-
    #   partykit::cforest_unbiased()
    mod <- party::cforest(form,
                           data = as.data.frame(traincoxphdata %>%
                             dplyr::mutate(
                              time   = !!time,
                              status = !!status) ),
                        cforest_control( ntree = 500) )
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
  timepoints <-  sort(unique(ndata %>% dplyr::select(!!time) %>% unlist))

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
  timepoints <-  sort(unique(ndata %>% dplyr::select(!!time) %>% unlist))

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

print_random_forest <- function(forest = bst){
  print(partykit::gettree(bst))
}
#' Var importance
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

var_importance <- function(forest = bst, type = "Permutation"){
  randomForest::importance(bst, type = 2)
}





#' Plot Cox Model Prediction
#'
#' This function takes a survival fit and plots, right now only works with factor
#'
#' @param
#' data : survival dataset \cr
#' ... : covariates to evaluate
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @importFrom rlang !!!
#' @import prodlim
#' @examples http://staff.pubhealth.ku.dk/~tag/Teaching/share/R-tutorials/SurvivalAnalysis.html

plot_cox_pred <- function(train_data=train, time = os_months, status = os_deceased,  ...){

  my_var <- dplyr::quos(...)
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  X <- train_data %>% dplyr::select(!!!my_var)

  form <- as.formula(paste("Hist(time, status)", paste("~", paste( colnames(X), collapse = " + "), sep = " " )))

  print(form)

  #library(pec)
  km.grade <- prodlim::prodlim(form,
    data = train_data %>% dplyr::mutate(
    time = !!time,
    status = !!status))

  cox.grade <- survival::coxph(form,
    data = train_data %>% dplyr::mutate(
    time = !!time,
    status = !!status) ,x=TRUE,y=TRUE)

  newdata <- data.frame(feature270 = c(1,0,-1),
                        feature350 = c(1,0,-1),
                        )
  ## first show Kaplan-Meier without confidence limits
  plot(km.grade, lty=1, lwd=3,
       col=c("darkgreen","darkorange","red"), confint=FALSE)
  ## now add survival estimates based on Cox regression
  pec::plotPredictSurvProb(cox.grade, lty=2,
                      col=c("darkgreen","darkorange","red"),
                      add=TRUE, sort(unique(train_data %>% dplyr::select(!!time) %>% unlist)),
                      newdata=newdata)
  mtext("Comparison of univariate  Cox regression and stratified Kaplan-Meier")
  }

