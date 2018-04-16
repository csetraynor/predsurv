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
fun_cv <- function(data, model_selection = "stepwise", time = os_months, status = os_deceased, K = 10){

  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  #Create Folds
  folds <- caret::createFolds(y = data %>% dplyr::select(!!status) %>% unlist , k = K)

  #Initiate ibrier vector
  ibrier <- rep(NA, K)

  #start cross-validation
  for(k in 1:K){

    #create train and test fold
    train <-  data[-folds[[k]] ,]
    test <-  data[folds[[k]] ,]

    #apply method
    if(model_selection == "stepwise"){
      mod.fs.coxph <-  predsurv::fun_fit(train, fit = "stepwise")
    }


    #make prediction
    if(model_selection == "stepwise"){
      pred.fs.coxph <- predsurv::mod_pred(obj = mod.fs.coxph, beta = coef, test_data = test, train_data = train)
    }

    #calculate Brier score
    if(model_selection == "stepwise"){
      brier_score <-  predsurv::fun_brier_score(obj = mod.fs.coxph,pred=pred.fs.coxph, test_data = test , beta = coef)
    }

    #save ibrier
    ibrier[k] <-  predsurv::fun_ibrier_score(brier_score)

    if(is.na(ibrier[k][1])){
      print(pred.fs.coxph);
      print(mod.fs.coxph);
      print(head(test[1:6,1:6]));
      print(head(train[1:6,1:6]));
      print(brier_score);
    }

    print(ibrier)

  } #end for loop end cross validation
  ibrier <- do.call(rbind, ibrier)

  return(ibrier)
}
