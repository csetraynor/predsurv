#' Clinical data from the TCGA dataset
#'
#' A dataset containing the survival outcome in lung cancer TCGA trial
#'
#' @format A data frame with230  rows and  19 variables:
#' \describe{
#'   \item{patient_id}{Randomly completely Anonamysed patient ID unique for each row}
#'   \item{cancer_type_detailed}{A longer explanation of cancer subtype}
#'   ...
#' }
#' @source \url{https://cancergenome.nih.gov/}
"clinicaldata"

#' Data Simulated for a Hypothetic Microarray experiment
#'
#' A dataframe containing the survival outcome and a gene expression matrix refere to data_sim_script to check the algorithm.
#'
#' @format Gene expression matrix with 1000 rows and 40 columns, survival outcome data frame with 40 rows and 2 variables (time, censoring).
#' \describe{
#'   \item{surv_months}{Actual survival months observed}
#'   \item{censoring.status}{censoring indicator 1 observed, 0 censored}
#'   \item{featurenames}{Hypothetic feature names in the gene expression matrix}
#'   ...
#' }
"survdata"
