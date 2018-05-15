#mc-cross-validation.R

#Load libraries
library(dplyr)
require(doMC)

#Load data
brca_data <- readRDS("brca_data.RDS")
mc_samples <- readRDS("mc_samples.RDS")

#Get brca iclust2
iclust2 <- brca_data[brca_data$intclust == 2,]
iclust2$intclust <- NULL
brca_data$intclust <- NULL
hold_out <- lapply(mc_samples, function(x)
  iclust2$patient_id[x])

#Train : Elastic net model
rocky <- function(holdout){
  #### obtain training datset
  if('patient_id' %in% colnames(survdata)){
    train <- survdata[!( (survdata %>%
                        dplyr::select(patient_id) %>%
                        unlist) %in% 
                          (holdout) ),];
    #unselect subject
    train <- train %>% dplyr::select(- patient_id)
  }
  
  # create predictor matrix
  x <- train %>% dplyr::select(-os_months, -os_deceased)
  
  ###### create fold id for CV in glmnet
  set.seed(9)
  foldid <-  caret::createFolds(train %>% select(os_deceased) %>% unlist,
                                k = 10, list = FALSE)  
  ###### Fit models
  ##################################################
  #grouped enet
  p.fac = rep(1, ncol(x))
  p.fac[match(c("npi","age_std"), colnames(x))] = 0
  #prepare 
  x <- as.matrix(x)
  y <- as.matrix(train %>%
                   dplyr::select(time = os_months,
                                 status = os_deceased), ncol = 2)
  require(doMC)
  registerDoMC(cores=nw)
  ### Apply elastic net with a=0.8 closer to Lasso
  mod <-  glmnet::cv.glmnet(x, y, family = "cox",
                            grouped = TRUE,
                            lambda.min.ratio = lambda, alpha = 0.8,
                            foldid = foldid, 
                            parallel = TRUE, penalty.factor = p.fac)
  # find optimised lambda
  optimal <- as.matrix(coef(mod, s = "lambda.min"))
  optimal <- as.data.frame(optimal)
  colnames(optimal) <- "mod"
  optimal$Hugo_symbol <- rownames(optimal)
  optimal <-  optimal %>% filter(mod != 0)
  return(optimal)
}

# This assumes that length(datlist) is much less than ncores
ncores <-  parallel::detectCores()
m <- 100
nw <- ncores %/% m

# survdata = iclust2
# lambda = 0.0001
# system.time(iclust2_fit <- mclapply(hold_out, rocky, mc.cores=m))
# 
# saveRDS(iclust2_fit, "iclust2_fit.RDS")

survdata = brca_data
lambda = 0.001
system.time(pooled_fit <- mclapply(hold_out, rocky, mc.cores=m ) )

saveRDS(pooled_fit, "pooled_fit.RDS")

print("Done!")

