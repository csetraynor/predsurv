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
"lungdata"

#' GBM data
#'
#'This is a subset of the glioblastoma dataset from the cancer genome atlas (TCGA) GBM study (2009) used in Shen et al. (2012).
#'
#' @format A list object containing three data matrices: copy number, methylation and mRNA expression in 55 samples.
#' \describe{
#'   \item{time}{Actual survival months observed}
#'   \item{cens}{censoring indicator 1 observed, 0 censored}
#'   \item{featurenames}{feature names in the gene expression matrix}
#'   ...
#' }
#' @references Ronglai Shen, Qianxing Mo, Nikolaus Schultz, Venkatraman E. Seshan, Adam B. Olshen, Jason Huse, Marc Ladanyi, Chris Sander. (2012). Integrative Subtype Discovery in Glioblastoma Using iCluster. PLoS ONE 7, e35236
"gbm"


