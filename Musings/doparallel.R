# mod1 <- foreach::foreach(independentvar = colnames(traincoxphdata %>%
#                                  dplyr::select(- !!time, - !!status)),
#                          .combine = list,
#                          .multicombine = TRUE,
#                          .packages=c('survival','dplyr','rlang'),
#                          .verbose = TRUE) %dopar%
#   coxph(as.formula(paste("Surv(time, status)"," ~ ", independentvar)),
#                           data = traincoxphdata %>%
#                             dplyr::mutate(time   = !!time,
#                                           status = !!status))
# doParallel::stopImplicitCluster()
mod1 <- lapply(colnames(trainX), FUN = reg, dep_var = "Surv(time, status)", data_source = traincoxphdata %>% dplyr::mutate(
  time   = !!time,
  status = !!status) )
print("2")
#Bonferroni correction
doParallel::registerDoParallel(cores = parallel::detectCores()-1)
mod2 <- foreach::foreach(model = mod1,
                         .combine = list,
                         .verbose=TRUE,
                         .multicombine = TRUE) %dopar%
  model[model$logtest["pvalue"] < 0.0005]$coefficients

doParallel::stopImplicitCluster()

doParallel::registerDoParallel(cl)
mod2[foreach::foreach(model = mod2,
                      .combine = c,
                      .verbose = TRUE) %dopar% is.null(model)] <- NULL
doParallel::stopImplicitCluster()

doParallel::registerDoParallel(cl)
features <- foreach::foreach(model = mod2,
                             .combine = c,
                             verbose = TRUE) %dopar%
  rownames(model)
doParallel::stopImplicitCluster()
print("4")
doParallel::registerDoParallel(cl)
coefficients <- foreach::foreach(model = mod2,
                                 .combine = c) %dopar%
  model[,1]
doParallel::stopImplicitCluster()



# reg <- function(indep_var,dep_var,data_source) {
#   formula <- as.formula(paste(dep_var," ~ ", indep_var))
#   res     <- survival::coxph(formula, data = data_source)
#   summary(res)
# }
doParallel::registerDoParallel(parallel::detectCores() - 1)
mod1 <- foreach::foreach(i = 1:ncol(traincoxphdata %>%
                                      dplyr::select(- !!time, - !!status)),
                         .combine = list,
                         .multicombine = TRUE,
                         .packages=c('survival','dplyr','rlang'),
                         .verbose = TRUE) %dopar%{
                           vars <- names((traincoxphdata %>%
                                            dplyr::select(- !!time, - !!status))[i])
                           coxmod <-  survival::coxph(as.formula(paste("Surv(time, status)"," ~ ", vars)),
                                                      data = traincoxphdata %>%
                                                        dplyr::mutate(time   = !!time,
                                                                      status = !!status))
                           coxmod
                         }
doParallel::stopImplicitCluster()
