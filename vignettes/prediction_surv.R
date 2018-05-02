library(predsurv)
theme_set(theme_bw())


####Simulation study
set.seed(12318)
survdata <- surv_sim_data(N = 1000, features = 100, CenRate = 1/10)

#Exploratory Analysis
tte <- plot_tte_dist(survdata, time = os_months, status = os_deceased)
km <- plot_km(survdata, time = os_months, status = os_deceased)
km
ggpubr::ggarrange(tte, km,
          labels = c("A", "B"),
          ncol = 2)

##Create train-test split
set.seed(83742)
fold = create_training_test_set(survdata, p = 0.5, status = os_deceased)
train = fold[["train"]]
test = fold[["test"]]
head(train[,1:6]); rm(fold)
##Train Models##
mod_uni <- predsurv::fun_train(train = train, fit = "Univariate")
mod_lasso <- predsurv::fun_train(train = train, fit = "Lasso")
mod_bst <- predsurv::fun_train(train = train, fit = "Random forest")
mod_ridge <- predsurv::fun_train(train = train, fit = "Ridge regression")
mod_enet <- predsurv::fun_train(train = train, fit = "Elastic net")
mod_iter_enet <- predsurv::fun_train(train = train, fit = "Elastic net", iterative = TRUE)
##ROC
roc_lasso <- predsurv::fun_test(obj = mod_lasso, train_data = train,  pred = "ROC", integrated = FALSE, noboot = 10, test_data = test)
roc_ridge <- predsurv::fun_test(obj = mod_ridge, train_data = train, pred = "ROC", integrated = FALSE, noboot = 10 , test_data = test)
roc_uni <- predsurv::fun_test(obj = mod_uni, train_data = train,  pred = "ROC", integrated = FALSE, noboot = 10 , test_data = test)
roc_enet <- predsurv::fun_test(obj = mod_enet, train_data = train,  pred = "ROC", integrated = FALSE, noboot = 10 , test_data = test)
roc_iter_enet <-  predsurv::fun_test(obj = mod_iter_enet, train_data = train,  pred = "ROC", integrated = FALSE, noboot = 10 , test_data = test)
##brier
brier_lasso <- predsurv::fun_test(obj = mod_lasso, train_data = train, test_data = test, pred = "Brier", integrated = FALSE)
brier_ridge <- predsurv::fun_test(obj = mod_ridge, train_data = train,  test_data = test ,pred = "Brier", integrated = FALSE)
brier_uni <- predsurv::fun_test(obj = mod_uni, train_data = train, test_data = test ,pred = "Brier", integrated = FALSE)
brier_enet <- predsurv::fun_test(obj = mod_enet, train_data = train, test_data = test, pred = "Brier", integrated = FALSE)
brier_bst<- predsurv::fun_test(obj = mod_bst, train_data = train, test_data = test, pred = "Brier", integrated = FALSE)
brier_iter_enet <- predsurv::fun_test(obj = mod_iter_enet, train_data = train, test_data = test, pred = "Brier", integrated = FALSE)

### Gives attr for plot
attr(roc_uni, 'prediction.of.model') <- "Uni"
attr(roc_lasso, 'prediction.of.model') <- "Lasso"
attr(roc_enet, 'prediction.of.model') <- paste0("ENet (a = ", attr(mod_enet, 'chosen.alpha'), ")")
attr(roc_iter_enet, 'prediction.of.model') <- paste0("Iter (a = ", attr(mod_iter_enet, 'chosen.alpha'), ")")
attr(roc_ridge, 'prediction.of.model') <- "Ridge"

### Gives attr for plot
attr(brier_uni, 'prediction.of.model') <- "Uni"
attr(brier_lasso, 'prediction.of.model') <- "Lasso"
attr(brier_enet, 'prediction.of.model') <- paste0("ENet (a = ", attr(mod_enet, 'chosen.alpha'), ")")
attr(brier_ridge, 'prediction.of.model') <- "Ridge"
attr(brier_bst, 'prediction.of.model') <- "Random forest"
attr(brier_iter_enet, 'prediction.of.model') <- paste0("Iter (a = ", attr(mod_iter_enet, 'chosen.alpha'), ")")

library(dplyr)
rocplot <- roc.plot2(roc_uni, roc_lasso, roc_enet, roc_ridge, roc_iter_enet) +
  labs(subtitle = paste0("Time = ", round(quantile(test %>%
                                                     dplyr::filter(os_deceased) %>%
                                                     dplyr::select(os_months) %>%
                                                     unlist, .73),0 ) , " months"))
rocplot


brierplot <- plot_brier2(brier_uni, brier_lasso, brier_enet, brier_bst, brier_iter_enet, brier_ridge) +
  labs(title = "Brier score", subtitle =
         paste0("Max time : " , round(max(test %>%
                                            dplyr::filter(os_deceased) %>%
                                            dplyr::select(os_months) %>%
                                          unlist )) , " months") )
brierplot

pdf('sim_results.pdf')
ggpubr::ggarrange(tte, km , brierplot, rocplot,
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
dev.off()


methods <- c("Univariate", "Lasso", "Ridge", "Elastic net", "Uni-elastice net", "Random forest")
brierList <- list(brier_uni, brier_lasso, brier_ridge, brier_enet, brier_iter_enet,  brier_bst)
ibrier <- lapply(brierList, function(i) pec::crps(i, models = "matrix"))
iBrier = sapply(seq_along(ibrier), function(i) round(ibrier[[i]][1], 2))

rocList <- list(roc_uni, roc_lasso, roc_ridge, roc_enet, roc_iter_enet)
AUC <-  unlist(c(sapply(seq_along(rocList), function(i) round(rocList[[i]]$AUC[1], 2)), NA))
       # "(", round(rocList[[i]]$AUC[3], 2),
       #   "-", round(rocList[[i]]$AUC[4], 2), ")" ) ), NA)

colAUC <- paste0("AUC", "(" ,"lower", "-" ,"upper" ,")")
tab_sim <- data.frame(Methods = methods,
           iBrier = iBrier,
           AUC = AUC)

library(stargazer)
stargazer(tab_sim, type = "latex", summary = FALSE,
          title = "Simulation study", digits = 2, digit.separate = 2)


# ##PCR
# mod.pcr = fun_train(train, fit = "PCR")
# pcr.brier_score = pred_error(obj = mod.pcr, fit = "PCR")
# plot_brier(pcr.brier_score)
#
# pcr.roc = pred_error(obj = mod.pcr, fit = "PCR", pred = "ROC")
# tdROC::plot.tdROC(pcr.roc)

##Real Data study
data("lungdata")
tte <- plot_tte_dist(lungdata, time = os_months, status = os_deceased)
km <- plot_km(lungdata, time = os_months, status = os_deceased)
lungdata <- tibble::rownames_to_column(lungdata, var = "subject")

#Standardise gene matrix
lungdata[,4:7132] <- std_dat(lungdata[,4:7132])


pairs(survdata[,3:10])

correlations <- cor(na.omit(lungdata))
# correlations
row_indic <- apply(correlations, 1, function(x) sum(x > 0.3 | x < -0.3) > 1)
correlations<- correlations[row_indic ,row_indic ]
pdf('corrplot.pdf')
corrplot::corrplot(correlations, method="square")
dev.off()


##Create train-test split
set.seed(83742)
fold = create_training_test_set(lungdata, status = os_deceased)
train = fold[["train"]]
test = fold[["test"]]
head(train[,1:6]); rm(fold)
##Train Models##
mod_uni <- predsurv::fun_train(train = train, fit = "Univariate")
mod_lasso <- predsurv::fun_train(train = train, fit = "Lasso")
mod_bst <- predsurv::fun_train(train = train, fit = "Random forest")
mod_ridge <- predsurv::fun_train(train = train, fit = "Ridge regression")
mod_enet <- predsurv::fun_train(train = train, fit = "Elastic net")
mod_iter_enet <- predsurv::fun_train(train = train, fit = "Elastic net",
                                     iterative = TRUE)

##ROC
roc_lasso <- predsurv::fun_test(obj = mod_lasso, train_data = train, data = lungdata, pred = "ROC", integrated = FALSE, noboot = 10)
roc_ridge <- predsurv::fun_test(obj = mod_ridge, train_data = train, data = lungdata, pred = "ROC", integrated = FALSE, noboot = 10)
roc_uni <- predsurv::fun_test(obj = mod_uni, train_data = train, data = lungdata, pred = "ROC", integrated = FALSE, noboot = 10)
roc_enet <- predsurv::fun_test(obj = mod_enet, train_data = train, data = lungdata, pred = "ROC", integrated = FALSE, noboot = 10)
roc_iter_enet <- predsurv::fun_test(obj = mod_iter_enet, train_data = train, data = lungdata, pred = "ROC", integrated = FALSE, noboot = 10)

##brier
brier_lasso <- predsurv::fun_test(obj = mod_lasso, train_data = train, data = lungdata, pred = "Brier", integrated = FALSE)
brier_ridge <- predsurv::fun_test(obj = mod_ridge, train_data = train, data = lungdata, pred = "Brier", integrated = FALSE)
brier_uni <- predsurv::fun_test(obj = mod_uni, train_data = train, data = lungdata, pred = "Brier", integrated = FALSE)
brier_enet <- predsurv::fun_test(obj = mod_enet, train_data = train, data = lungdata, pred = "Brier", integrated = FALSE)
brier_iter_enet <- predsurv::fun_test(obj = mod_iter_enet, train_data = train, data = lungdata, pred = "Brier", integrated = FALSE)
brier_bst<- predsurv::fun_test(obj = mod_bst, train_data = train, data = lungdata, pred = "Brier", integrated = FALSE)

### Gives attr for plot
attr(roc_uni, 'prediction.of.model') <- "Uni"
attr(roc_lasso, 'prediction.of.model') <- "Lasso"
attr(roc_enet, 'prediction.of.model') <- paste0("ENet (a = ", attr(mod_enet, 'chosen.alpha'), ")")
attr(roc_iter_enet, 'prediction.of.model') <- paste0("Iter (a = ", attr(mod_iter_enet, 'chosen.alpha'), ")")
attr(roc_ridge, 'prediction.of.model') <- "Ridge"

### Gives attr for plot
attr(brier_uni, 'prediction.of.model') <- "Uni"
attr(brier_lasso, 'prediction.of.model') <- "Lasso"
attr(brier_enet, 'prediction.of.model') <- paste0("ENet (a = ", attr(mod_enet, 'chosen.alpha'), ")")
attr(brier_ridge, 'prediction.of.model') <- "Ridge"
attr(brier_bst, 'prediction.of.model') <- "Random forest"
attr(brier_iter_enet, 'prediction.of.model') <- paste0("Iter (a = ", attr(mod_iter_enet, 'chosen.alpha'), ")")

library(dplyr)
rocplot <- roc.plot2(roc_uni, roc_lasso, roc_enet, roc_ridge, roc_iter_enet) +
  labs(subtitle = paste0("Time = ", round(quantile(test %>%
                                                     dplyr::filter(os_deceased == 1) %>%
                                                     dplyr::select(os_months) %>%
                                                     unlist, .73),0 ) , " months"))
rocplot


brierplot <- plot_brier2(brier_uni, brier_lasso, brier_enet, brier_bst, brier_iter_enet, brier_ridge) +
  labs(title = "Brier score", subtitle =
         paste0("Max time : " , round(max(test %>%
                                            dplyr::filter(os_deceased == 1) %>%
                                            dplyr::select(os_months) %>%
                                            unlist )) , " months") )
brierplot
# plot(my_trained_models[[3]])

pdf('lung_results.pdf')
ggpubr::ggarrange(tte, km , brierplot, rocplot,
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
dev.off()



methods <- c("Univariate", "Lasso", "Ridge", "Elastic net", "Uni-elastice net", "Random forest")
brierList <- list(brier_uni, brier_lasso, brier_ridge, brier_enet, brier_iter_enet,  brier_bst)
ibrier <- lapply(brierList, function(i) pec::crps(i, models = "matrix"))
iBrier = sapply(seq_along(ibrier), function(i) round(ibrier[[i]][1], 2))

rocList <- list(roc_uni, roc_lasso, roc_ridge, roc_enet, roc_iter_enet)
AUC <-  unlist(c(sapply(seq_along(rocList), function(i) round(rocList[[i]]$AUC[1], 2)), NA))
# "(", round(rocList[[i]]$AUC[3], 2),
#   "-", round(rocList[[i]]$AUC[4], 2), ")" ) ), NA)

colAUC <- paste0("AUC", "(" ,"lower", "-" ,"upper" ,")")
tab_lung <- data.frame(Methods = methods,
                      iBrier = iBrier,
                      AUC = AUC)

library(stargazer)
stargazer(tab_lung, type = "latex", summary = FALSE,
          title = "Simulation study", digits = 2, digit.separate = 2)
