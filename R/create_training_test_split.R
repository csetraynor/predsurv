#' Create training test split
#'
#' This function plots the TTE distribution
#'
#' @param
#' d a dataset \cr
#' seed  \cr
#' percent the percentage of data that goes to training
#' @return a list with train and test set
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
create_training_test_set <- function(d, status = os_deceased,  seed = 111, percent = 0.8){
  set.seed(seed)

  status <- dplyr::enquo(status)

  tmp <- rsample::mc_cv(d, strata = d %>% dplyr::select(!!status)%>%colnames , times = 1, prop = 3/4)
  train_data <- purrr::map(tmp$splits,
                        function(x) {
                        as.data.frame(x)})
  train <- train_data[[1]];
  test <-  d[!(d$subject %in% train$subject ),];

  # train <- train %>% dplyr::select(- !!y)
  # test <- test %>% dplyr::select(- !!y)


  return(list(train = train, test = test))
}
