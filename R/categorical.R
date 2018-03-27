#' Coerce categorical variables in index
#'
#' This function allows to coerce categorical variables into index, \cr
#' @param
#' d a dataset \cr
#' ... vars to be coerced \cr
#' @return d a dataset with coerced covariates
#' @export
categories_to_index <- function(d, ...){
  fit_var <- dplyr::quos(...)
  d <- d %>%
    mutate_at(vars(!!!fit_var), rethinking::coerce_index)
  return(d)
}
