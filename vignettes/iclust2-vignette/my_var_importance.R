mc_samp <- readRDS("C:/RFactory/data_test.RDS")
mc_samp <- readRDS("/home/mtr/rfactory/data_test.RDS")

survdata <- readRDS("my_musings/mc_samp.RDS")
survdata$splits <- mc_samp$splits


#Calculate Brier Skill Score

survdata$bss <- purrr::map2(survdata$brier_ic2, survdata$brier_reference, bss)

library(tidyr)
par(mfrow=c(1,1))


features <- purrr::map(survdata$mod_lasso, function(x){
  output <- as.data.frame(x);
  data.frame( features = rownames(output),
              coefficient = output %>% unlist,
              number = 1)
}
  )

summary_features <- do.call(rbind, features)

count_features <- dplyr::count(summary_features, features) %>%
  arrange(desc(n))






