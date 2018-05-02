###iClust2 Analyses

library(predsurv)
data("metabric_clinical_data")

##### Data pre-process #######
#### Clean-data    #######
metabric_clinical_data <- blank_to_na(metabric_clinical_data)
plot_na(metabric_clinical_data)
metabric_clinical_data <- drop_na(metabric_clinical_data, os_status)
plot_na(metabric_clinical_data)

# Get genomic data
##### Can be downloaded from http://www.cbioportal.org/study?id=brca_metabric#summary to download, follow the link-> Download data
genedata <- readr::read_tsv("C:/RFactory/Downloads/brca_metabric/data_expression.txt", col_names = TRUE)
object.size(genedata)
###Which is too heavy to upload

##Impute data
sum(is.na(genedata[-2])) ## 1 and 2 are Hugo_Symbol and Entrez_Gene_Id
genedata[,c(-1,-2)] <- impute_gene_matrix(genedata[,c(-1,-2)])
sum(is.na(genedata[-2])) ## 1 and 2 are Hugo_Symbol and Entrez_Gene_Id

###stick vars
brca <- stick_vars(clidata = metabric_clinical_data, genedata = genedata)
###Remove gene matrix but save gene names
Hugo_Symbol <- unlist(genedata[,"Hugo_Symbol"])
Entrez_Gene_Id <- unlist(genedata[,"Entrez_Gene_Id"])
rm(list = c("genedata", "metabric_clinical_data"))

#Standardise predictors
brca[,c("age_at_diagnosis","npi",Hugo_Symbol)] <- std_dat(brca[,c("age_at_diagnosis","npi",Hugo_Symbol)])

#Create os_deceased and remove os_status
unique(brca$os_status)
brca$os_deceased <- brca$os_status == "DECEASED"

#Get brca iclust2
iclust2 <- brca[brca$intclust == 2,]

#create resamples first
set.seed(9666)
mc_samp <- rsample::mc_cv(iclust2, strata = "os_deceased", times = 100, prop = 1/4)

####Train models with iclust2
mc_samp$mod_lasso <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Lasso", data = iclust2)
saveRDS(mc_samp$mod_lasso, "mod_lasso2.RDS")
mc_samp$mod_lasso <- NULL

mc_samp$mod_ridge <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Ridge regression", data = iclust2)
saveRDS(mc_samp$mod_ridge, "mod_ridge2.RDS")
mc_samp$mod_ridge <- NULL

mc_samp$mod_enet <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Elastic net", data = iclust2)
saveRDS(mc_samp$mod_enet, "mod_enet2.RDS")
mc_samp$mod_enet <- NULL

mc_samp$mod_srf <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Random forest", data = iclust2)
saveRDS(mc_samp$mod_srf, "mod_srf2.RDS")
mc_samp$mod_srf <- NULL

####Train models with pooled data
mc_samp$mod_lasso <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Lasso", data = brca)
saveRDS(mc_samp$mod_lasso, "mod_lasso_pool.RDS")
mc_samp$mod_lasso <- NULL

mc_samp$mod_ridge <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Ridge regression", data = brca)
saveRDS(mc_samp$mod_ridge, "mod_ridge_pool.RDS")
mc_samp$mod_ridge <- NULL

mc_samp$mod_enet <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Elastic net", data = brca)
saveRDS(mc_samp$mod_enet, "mod_enet_pool.RDS")
mc_samp$mod_enet <- NULL

mc_samp$mod_srf <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Random forest", data = brca)
saveRDS(mc_samp$mod_srf, "mod_srf_pool.RDS")
mc_samp$mod_srf <- NULL

### Perform adap

