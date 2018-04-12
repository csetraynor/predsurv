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
create_training_test_set <- function(d, y = censoring.status,  seed = 111, percent = 0.8){
  set.seed(seed)

  status <- dplyr::enquo(y)

  tmp <- caret::createDataPartition(y = d %>% select(!!status) %>% unlist,
                                    p = percent, list =  FALSE)
  train <-  d[tmp,] ; test <- d[-tmp,]
  return(list(train = train, test = test))
}
