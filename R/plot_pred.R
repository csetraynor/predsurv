#' Plot Brier Score
#'
#' Plot Brier Score
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
plot_brier <- function(brier_score){
  Brier_Score = brier_score$AppErr[[2]]
  Time = brier_score$time
  plot(Brier_Score ~ Time, type= "l")

}
