library(qusage)

hallmark_all <- read.gmt("/home/mtr/rfactory/predsurv/vignettes/Gene sets/h.all.v6.1.symbols.gmt")


results <- read.csv("/home/mtr/rfactory/predsurv/vignettes/iclust2-vignette/summary_iclust2.csv")


ontology_search(hallmark_list = hallmark_all, predictor = results)
