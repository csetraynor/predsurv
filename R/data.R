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

#' Microarray Dataset
#'
#' Lung adenocarcinomas were profiled by Beer et al. (2002) using Affymetrix hu6800 microarrays. The data here were normalized from raw .CEL files by RMAExpress (v0.3). The expression matrix contains expression data for 86 patients with 7,129 probe sets.
#'
#' @format Gene expression matrix and survival outcome data frame (time, cens).
#' \describe{
#'   \item{time}{Actual survival months observed}
#'   \item{cens}{censoring indicator 1 observed, 0 censored}
#'   \item{featurenames}{feature names in the gene expression matrix}
#'   ...
#' }
"survdata"
