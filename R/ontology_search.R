#' Examine Ontology
#'
#' Check list of Ontology to find genes relevant assessed by cancer biology.
#' @param
#' hallmark_list: list of hallmarks from ontology database \cr
#' predictor: feature list \cr
#' coeff : TRUE return hazard coefficient
#'
#' @return predicted linear predictor for each individual in the test set
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!

ontology_search <- function(hallmark_list, predictor, coeff = TRUE){

  p <- as.character(predictor$Predictor)
  hall <- lapply(hallmark_list, function(x){
    sapply(p, function(p)
      ifelse(p %in% x, p, NA)
    )
  }
  )
  hall <- lapply(hall, function(x)  x[!is.na(x)])
  hall <- hall[sapply(hall, length)>0]

  if(coeff){
    hall <- lapply(hall, function(h){
      sapply(h, function(i){
        predictor %>% filter(Predictor == i) %>% dplyr::select(median_coeff)
      })
    })
    print(hall)
    hall <- lapply(hall, function(h){
      data.frame( features =  gsub("\\..*","",names(h) ),
                  hazard_coeff = unlist(h), row.names = NULL)
    })
  }else{
    hall <- lapply(hall, function(h) unname(h))
  }
  return(hall)
}
