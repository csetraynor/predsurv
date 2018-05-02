######
# #download data
library(cgdsr)
library(dplyr)
library(readr)

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)
# Get list of cancer studies at server
list_cancer = getCancerStudies(mycgds)
# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[25,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
# Get clinical data for the case list
clinicaldata = getClinicalData(mycgds,mycaselist)

# Get an easier code to read
names(clinicaldata) <- tolower(names(clinicaldata))
clinicaldata <- tibble::rownames_to_column(clinicaldata, var = "patient_id")
clinicaldata$patient_id <- gsub( "\\.","-", clinicaldata$patient_id)
glimpse(clinicaldata)
clinicaldata <- clinicaldata %>%
  select(patient_id, intclust, npi, age_at_diagnosis, os_months, os_status)
metabric_clinical_data <- clinicaldata
devtools::use_data(metabric_clinical_data)
# Get genomic data
##### Can be downloaded from http://www.cbioportal.org/study?id=brca_metabric#summary to download, follow the link-> Download data
genedata <- readr::read_tsv("C:/RFactory/Downloads/brca_metabric/data_expression.txt", col_names = TRUE)
object.size(gendata)
###Which is too heavy to upload

rm(list=ls())

