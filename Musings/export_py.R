
####### Export to Python

hold_out <- as.data.frame(mc_samp$splits$`1`)
library(dplyr)
x <- brca[!( (brca %>%
                dplyr::select(patient_id) %>%
                unlist) %in% (hold_out %>%
                                dplyr::select(patient_id) %>%
                                unlist) ),];
x <- x %>% select(-patient_id)
y <- x %>% select(os_months, os_deceased)
X <- x %>% select(-os_months, -os_deceased)

write.csv(X, "X_split_1.csv", row.names = FALSE)
write.csv(y, "y_split_1.csv", row.names = FALSE)

##### create penalty.factor vector
p.fac = rep(1, ncol(X))
p.fac[match(c("npi","age_std"), colnames(X))] = 0

write.csv(as.matrix(p.fac), "p.fac_split_1.csv", row.names = FALSE)