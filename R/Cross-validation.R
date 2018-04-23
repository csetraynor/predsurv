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
fun_cv <- function(..., data, iter = FALSE, time = os_months, status = os_deceased, KMC = 20){
  M <- list(...)

  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  #Create folds
  set.seed(9)
  mc_samp <- rsample::mc_cv(data, strata = "os_deceased", times = KMC)
  train_cv <- purrr::map(mc_samp$splits,
             function(x) {
               as.data.frame(x)})

  #start cross-validation
  out <- rep(NA,length(M))
  out <- for(m in seq_along(M)){ lapply(seq_along(train_cv), function(k){
    print(k)
    #create train and test fold
    train <- train_cv[[k]];

    test <-  data[setdiff(data$subject, train$subject ),];

    trained_model <- predsurv::fun_train(train = train, fit = M[[m]], iterative = iter)
    test_model <- predsurv::fun_test(test_data = test, train_data = train, obj = trained_model, all = TRUE)
  })}
  return(out)
}
