#' Compute C-V iBrier
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
fun_cv <- function(data, method = "forward", time = os_months, status = os_deceased, K = 10){

  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  #Create Folds
  folds <- caret::createFolds(y = data %>% dplyr::select(!!status) %>% unlist , k = K)

  #Initiate ibrier vector
  ibrier <- rep(NA, K)

  #start cross-validation
  for(k in 1:K){

    #create train and test fold
    train = data[-folds[[k]] ,]
    test = data[folds[[k]] ,]

    #apply method
    if(method == "forward"){
      mod.fs.coxph = predsurv::fun_fit_fs(train)
    }


    #make prediction
    if(method == "forward"){
    pred.forward.coxph = predsurv::mod_pred(obj = mod.fs.coxph, beta = coef)
    }

    #calculate Brier score
    if(method == "forward"){
    brier_score = predsurv::fun_brier_score(obj = mod.fs.coxph, pred = pred.forward.coxph, data = test , beta = coef)
    }

    #save ibrier
    ibrier[k] = predsurv::fun_ibrier_score(brier = brier_score)

  } #end for loop end cross validation

  return(ibrier)
}
