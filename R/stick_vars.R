#' Stick genes
#'
#' This function makes easy to stick the gene vars to the clinical vars
#' @param
#' clidata : clinical dataset \cr
#' genedata: genomic dataset \cr
#' @return a full dataset
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
stick_vars <- function(clidata, genedata){

  #get gene expression matrix obtain id
  gene_expression <- genedata %>% dplyr::select(- Hugo_Symbol, -Entrez_Gene_Id)
  Hugo_Symbol <- genedata %>% dplyr::select(Hugo_Symbol) %>% unlist
  id_gene_expression <- colnames(gene_expression)
  #from clinical data keep observations with gene expression
  #obtain sample with gene expression measurements and match
  Y <- clidata[clidata$patient_id %in% id_gene_expression,]
  Y <- arrange(Y,  match(Y$patient_id, id_gene_expression))
  assertthat::assert_that(all(id_gene_expression == Y$patient_id))

  ##Transpose gene expression matrix
  gene_matrix <- as.matrix(gene_expression, ncol = ncol(gene_expression))
  t_gene_matrix <- t(gene_matrix)

  ##Get dataframe back
  t_gene_expression <- as.data.frame(t_gene_matrix)

  colnames(t_gene_matrix) <- Hugo_Symbol

  ### cbind clinical and genomic data
  out <- cbind(Y, t_gene_matrix)

  return(out)
}

