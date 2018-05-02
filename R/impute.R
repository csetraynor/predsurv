#' Data Impute
#'
#' This function imputes a gene matrix, currently with knn
#' @param
#' genedata: genomic dataset \cr
#' @return a full dataset
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
impute_gene_matrix <- function(genedata){
  genedata <- as.matrix(genedata)
  genes_impute <- impute::impute.knn(genedata, k = 10)$data

  return(genes_impute)
}
