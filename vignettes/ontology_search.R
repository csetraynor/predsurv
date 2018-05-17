# !diagnostics off
library(qusage)

hallmark_all <- qusage::read.gmt("/home/mtr/rfactory/predsurv/vignettes/Gene sets/h.all.v6.1.symbols.gmt")

hallmark_all <- qusage::read.gmt("C:/RFactory/predsurv/vignettes/Gene sets/c6.all.v6.1.symbols.gmt")

results <- read.csv("/home/mtr/rfactory/predsurv/vignettes/iclust2-vignette/summary_iclust2.csv")

results <- read.csv("C:/RFactory/predsurv/performance_results/summary_iclust2.csv")

onco_signature <- ontology_search(hallmark_list = hallmark_all, predictor = results,
                             coef = median_coeff)

onco_signature_weight <- ontology_search(hallmark_list = hallmark_all, predictor = results,
                                  coef = weight)

capture.output(print(onco_signature), file = "onco_signature.txt")


library(glmnet)
enet_model <- readRDS("C:/RFactory/predsurv/performance_results/enet_model_iclust2_filter.RDS")

data("geneList")
genedata <- readr::read_tsv("C:/RFactory/Downloads/brca_metabric/data_expression.txt", col_names = TRUE)

genedata$Entrez_Gene_Id[match(names(geneList), genedata$Hugo_Symbol)]

geneList <- coef(clinico_genomic_fit)[1:10]



names(geneList) <- genedata$Entrez_Gene_Id[match(names(geneList), genedata$Hugo_Symbol)]

sort(geneList, decreasing = TRUE)


gseDO(geneList = sort(geneList, decreasing = TRUE),
      nPerm = 100,
      minGSSize = 1,
      pvalueCutoff = 0.2,
      pAdjustMethod = "BH",
      seed = 10)
