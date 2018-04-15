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
create_training_test_set <- function(d, status = status,  seed = 111, percent = 0.8){
  set.seed(seed)

  y <- dplyr::enquo(status)

  tmp <- caret::createDataPartition(y = d %>% dplyr::select(!!y) %>% unlist,
                                    p = percent, list =  FALSE)
  train <-  d[tmp,]
  test <- d[-tmp,]

  # train <- train %>% dplyr::select(- !!y)
  # test <- test %>% dplyr::select(- !!y)


  return(list(train = train, test = test))
}
