library(qusage)

hallmark_all <- read.gmt("/home/mtr/rfactory/predsurv/vignettes/Gene sets/h.all.v6.1.symbols.gmt")


results <- read.csv("/home/mtr/rfactory/predsurv/vignettes/iclust2-vignette/summary_iclust2.csv")






hall <- lapply(hallmark_all, function(x){
  sapply(results$Predictor, function(p)
    ifelse(p %in% x, p, NA)
  )
}
  )
hall <- lapply(hall, function(x)  x[!is.na(x)])
hall <- hall[sapply(hall, length)>0]

hall2 <- lapply(hall, function(h){
  sapply(h, function(i){
    results %>% filter(Predictor == i) %>% dplyr::select(median_coeff)
  })
})


hall2 <- lapply(hall2, function(h){
  data.frame( features =  gsub("\\..*","",names(h) ),
              hazard_coeff = unlist(h), row.names = NULL)
})



#Filter(Negate(is.null), hall)


