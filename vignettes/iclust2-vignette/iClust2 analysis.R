###iClust2 Analyses

library(predsurv)
#devtools::document()
data("metabric_clinical_data")

##### Data pre-process #######
#### Clean-data    #######
metabric_clinical_data <- blank_to_na(metabric_clinical_data)
plot_na(metabric_clinical_data)
metabric_clinical_data <- drop_na(metabric_clinical_data, os_status)
plot_na(metabric_clinical_data)
metabric_clinical_data$os_months <- ifelse(metabric_clinical_data$os_months == 0, 1e-10, metabric_clinical_data$os_months)

# Get genomic data
##### Can be downloaded from http://www.cbioportal.org/study?id=brca_metabric#summary to download, follow the link-> Download data
genedata <- readr::read_tsv("C:/RFactory/Downloads/brca_metabric/data_expression.txt", col_names = TRUE)
genedata <- readr::read_tsv("/home/mtr/rfactory/brca_metabric/data_expression.txt", col_names = TRUE)
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
brca[,Hugo_Symbol] <- std_dat(brca[,Hugo_Symbol])
brca$age_std <- (brca$age_at_diagnosis - mean(brca$age_at_diagnosis)) / sd(brca$age_at_diagnosis)
brca$age_at_diagnosis <- NULL

#Create os_deceased and remove os_status
unique(brca$os_status)
brca$os_deceased <- brca$os_status == "DECEASED"
brca$os_status <- NULL

brca <- readRDS("/media/mtr/SeagateExpansionDrive/R-Factory/brca_data.RDS")

#Get brca iclust2
iclust2 <- brca[brca$intclust == 2,]
iclust2$intclust <- NULL
brca$intclust <- NULL


#create resamples first
set.seed(9666)
mc_samp <- rsample::mc_cv(iclust2, strata = "os_deceased", times = 100, prop = 1/4)

mc_samples <- lapply(mc_samp$splits, function(s) s$in_id)

saveRDS(mc_samples, "mc_samples.RDS")

memory.limit(10e10)

####Train models with pooled data
mc_samp$mod_enet_pool <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Elastic net", data = brca, lambda = 0.0001)
saveRDS(mc_samp$mod_enet_pool, "mod_enet_pool.RDS")
mc_samp$mod_enet_pool <- NULL

####Train models with iclust2
mc_samp$mod_enet_2 <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Elastic net", data = iclust2, lambda = 0.001)
saveRDS(mc_samp$mod_enet_2, "mod_enet_2.RDS")
mc_samp$mod_enet_2 <- NULL



memory.limit(10e10)
####Train models with iclust2
mc_samp$mod_lasso2 <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Lasso", data = iclust2, lambda = 0.001)
saveRDS(mc_samp$mod_lasso2, "mod_lasso2.RDS")
mc_samp$mod_lasso2 <- NULL

####Train models with pooled data
mod_lasso_pool <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Lasso", data = brca, lambda = 0.001)
saveRDS(mod_lasso_pool, "mod_lasso_pool.RDS")
mod_lasso_pool <- NULL

memory.limit(10e10)
####Train models with iclust2
mc_samp$mod_lasso2 <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Lasso", data = iclust2, lambda = 0.01)
saveRDS(mc_samp$mod_lasso2, "mod_lasso2.RDS")
mc_samp$mod_lasso2 <- NULL

####Train models with pooled data
mod_lasso_pool <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Lasso", data = brca, lambda = 0.001)
saveRDS(mod_lasso_pool, "mod_lasso_pool.RDS")
mod_lasso_pool <- NULL

mc_samp$mod_lasso2 <- readRDS("mod_lasso2.RDS")

####Train models with iclust2
mc_samp$mod_ridge <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Ridge regression", data = iclust2, lambda = 0.01)
saveRDS(mc_samp$mod_ridge, "mod_ridge2.RDS")
mc_samp$mod_ridge <- NULL

mc_samp$mod_enet <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Elastic net", data = iclust2, lambda = 0.01)
saveRDS(mc_samp$mod_enet, "mod_enet2.RDS")
mc_samp$mod_enet <- NULL

mc_samp$mod_srf <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Random forest", data = iclust2)
saveRDS(mc_samp$mod_srf, "mod_srf2.RDS")
mc_samp$mod_srf <- NULL

####Train models with pooled data
mc_samp$mod_lasso <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Lasso", data = brca, lambda = 0.001)
saveRDS(mc_samp$mod_lasso, "mod_lasso_pool.RDS")
mc_samp$mod_lasso <- NULL

mc_samp$mod_ridge <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Ridge regression", data = brca, lambda = 0.001)
saveRDS(mc_samp$mod_ridge, "mod_ridge_pool.RDS")
mc_samp$mod_ridge <- NULL

mc_samp$mod_enet <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Elastic net", data = brca, lambda = 0.001)
saveRDS(mc_samp$mod_enet, "mod_enet_pool.RDS")
mc_samp$mod_enet <- NULL

mc_samp$mod_srf <- purrr::map(mc_samp$splits, predsurv::fun_train2, fit = "Random forest", data = brca)
saveRDS(mc_samp$mod_srf, "mod_srf_pool.RDS")
mc_samp$mod_srf <- NULL

### Perform adap

mc_samp$mod_ridge <- readRDS("mod_ridge2.RDS")
mc_samp$mod_lasso2 <- readRDS("mod_lasso2.RDS")

library(doMC)
?mclapply()




########## Performance measure

mc_cv <-  readRDS("mc_cv.RDS")
saveRDS(mc_cv, "mc_cv.RDS")

mc_cv <- readRDS("C:/RFactory/mc_cv.RDS")


mc_cv$brier_enet_iclust2 <- purrr::pmap_dbl(list(mc_cv$splits,
                                                   mc_cv$enet_iclust2),
                                       function(splits, mod){predsurv::fun_mp(
                                         obj = mod,
                                         test_data = splits,
                                         fit = "Elastic net",
                                         pred = "Brier",
                                         subject = patient_id
                                       )
                                       })


saveRDS(mc_cv, "mc_cv.RDS")

mc_cv$brier_enet_lasso <- purrr::pmap_dbl(list(mc_cv$splits, mc_cv$lasso_iclust2),
                                           function(splits, mod){predsurv::fun_mp(
                                             obj = mod,
                                             test_data = splits,
                                             fit = "Lasso",
                                             pred = "Brier",
                                             subject = patient_id
                                           )
                                           })


mc_cv$brier_reference <- purrr::pmap_dbl(list(mc_cv$splits, mc_cv$lasso_iclust2),
                                          function(splits, mod){predsurv::fun_mp(
                                            obj = mod,
                                            test_data = splits,
                                            fit = "Lasso",
                                            pred = "Brier",
                                            subject = patient_id,
                                            reference = TRUE
                                          )
                                          })

saveRDS(mc_cv, "mc_cv.RDS")

library(tidyr)
library(ggplot2)
##Plot Brier

brier_dens <- mc_cv %>%
  dplyr::select(contains("brier")) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "bottom")
brier_dens <- brier_dens +
  labs(x = "Integrated Brier Score",
       title = "Density of iBrier") +
  geom_vline(xintercept =  0.25 , linetype = "dotted" )
brier_dens

library(tidyposterior)
mc_samp_brier <- tidyposterior::perf_mod(mc_samp %>%
                                           dplyr::select(-splits, - mod_lasso),
                                         seed = 6507, iter = 5000, transform = logit_trans,
                                         hetero_var = FALSE)

mbri_tab <- summary(tidy(mc_samp_brier))
mbri_tab <- as.data.frame(mbri_tab)

star = stargazer(mbri_tab, type = "latex", summary = FALSE, digits.extra = 3,digits = 3, digit.separator = ".",
                 title = "Bayesian analysis of resampling AUC")


stargazer(mbri_tab, type = "latex", summary = FALSE, digits.extra = 3,
          digits = 3, digit.separator = ".",
          title = "Bayesian analysis of resampling AUC")


posterior_brier <- ggplot(tidy(mc_samp_brier)) +
  theme_bw()+
  labs(
    title = "Posterior probability for integrated Brier Score")
posterior_brier <- posterior_brier +   labs(
  title = "Posterior probability of iBrier")
posterior_brier

comparisons_brier <- contrast_models(
  mc_samp_brier,
  list_1 = rep("brier_reference", 1),
  list_2 = "brier_lasso",
  seed = 4654
)

compare_brier <- ggplot(comparisons_brier, size = 0.05) +
  theme_bw()+
  labs(
    title = "Posterior probability of iBrier.",
    subtitle ="Benchmark: Reference")

compare_brier <- compare_brier +   labs(
  title = "Posterior probability for iBrier",
  subtitle ="Benchmark: Reference")
compare_brier


####### Variable Importance


#Calculate Brier Skill Score

mc_cv$bss_enet <- purrr::map2_dbl(mc_cv$brier_enet_iclust2, mc_cv$brier_reference, bss)
mc_cv$bss_lasso <- purrr::map2_dbl(mc_cv$brier_lasso_iclust2, mc_cv$brier_reference, bss)



library(tidyr)
library(ggplot2)
##Plot Brier

bss_dens <- mc_cv %>%
  dplyr::select(contains("bss")) %>%
  gather() %>%
  ggplot(aes(x = value, col = key, fill = key)) +
  geom_density(alpha = 0.1) +
  theme_bw() +
  theme(legend.position = "bottom")
bss_dens <- bss_dens +
  labs(x = "Brier Skill Score",
       title = "Density Brier Skill score") +
  geom_vline(xintercept =  0 , linetype = "dotted" )
bss_dens


############# Model averaging
mc_cv$weight_enet <- bss(mc_cv$brier_enet_iclust2, mc_cv$brier_reference, averaging = TRUE)

mc_cv %>%
  dplyr::select(contains("weight")) %>%
  gather() %>%
  ggplot(aes(x = value, col = key, fill = key)) +
  geom_density(alpha = 0.1) +
  theme_bw()
############ Select features

par(mfrow=c(1,1))


summary_features <- do.call(rbind, mc_cv$enet_iclust2)
summary_features$split <- as.numeric(gsub("\\..*","",rownames(summary_features)))

summary_features$weights <- mc_cv$weight_enet[summary_features$split]*100

summary_weights <- summary_features %>%
  group_by(Predictor) %>%
  dplyr::summarize(
    count = n(),
    weight = sum(weights),
    min_coeff = min(Coefficient),
    max_coeff = max(Coefficient),
    median_coeff = median(Coefficient),
    weight_coeff = sum(weights*Coefficient)/100) %>%
  arrange(desc(abs(weight_coeff)))

    ,
    coeff_weight = sum(weights*coefficients()))

library(data.table)
summary_features <- rbindlist(lapply(mc_cv$enet_iclust2, FUN = function(x) x))



count_features <- dplyr::count(summary_features, Predictor) %>%
  arrange(desc(n))

summary_features %>%
  filter(Predictor == "NGF")

summary_features %>%
  filter(Predictor == "MAP1B")





#####GSEA

require(clusterProfiler)

