#' Convert all blank to na in dataset
#'
#' This function converts all blank values in the dataset to NA,
#' this is very useful in real world datasets that contain masked NA values
#' @param  d a dataset from cgdsr
#' @return dataset complete
#' @export
blank_to_na <- function(d){
  `%>%` <- magrittr::`%>%`
  convert_blank_to_na <- function(x) {
    if(!purrr::is_character(x)){
      return(x)
      print("Error not character")
    } else {
      ifelse(x == "" | x == "[Not Available]" , NA, x)
    }
  }
  d <- d %>% dplyr::mutate_all(dplyr::funs(convert_blank_to_na))
  return(d)
}
#' Plot all na
#' This function gives an explanatory plot of NA values using VIM package
#' @param  d a dataset from cgdsr
#' @return dataset complete
#' @export
plot_na <- function(d){
    d %>%
    VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE,
                sortVars = TRUE, sortCombs = TRUE, plot = TRUE, only.miss = FALSE)
}
#' Clean dataset
#'
#' This function is only to be used if there is the will to delete missing data.
#' All rows with missing data from the vars passed will be deleted
#' @param
#' d a dataset \cr
#' ... vars to be evaluated \cr
#' n optional parameter to assure number of rows deleted, default NA
#' @return d a clean dataset
#' @export
drop_na <- function(d, ..., n = NA){
  fit_var <- dplyr::quos(...)
  d_old <- d
  d <- d[complete.cases(d %>% dplyr::select(!!!fit_var)),]
  # d <- d %>%
  #   dplyr::filter(!is.na(!!fit_var))
  if(!is.na(n)){
    #Check n fewer obsrvations than original
    assertthat::assert_that(nrow(d) == nrow(d_old) - n)
  }
  return(d)
}


