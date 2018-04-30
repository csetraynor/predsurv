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
fold = create_training_test_set(survdata, p = 0.8, status = os_deceased)
train = fold[["train"]]
test = fold[["test"]]
head(train[,1:6]); rm(fold)
##Train Models##
mod_uni <- predsurv::fun_train(train = train, fit = "Univariate")
mod_lasso <- predsurv::fun_train(train = train, fit = "Lasso")
mod_bst <- predsurv::fun_train(train = train, fit = "Random forest")
mod_ridge <- predsurv::fun_train(train = train, fit = "Ridge regression")
mod_enet <- predsurv::fun_train(train = train, fit = "Elastic net")

##ROC
roc_lasso <- predsurv::fun_test(obj = mod_lasso, train_data = train,  pred = "ROC", integrated = FALSE, noboot = 10, test_data = test)
roc_ridge <- predsurv::fun_test(obj = mod_ridge, train_data = train, pred = "ROC", integrated = FALSE, noboot = 10 , test_data = test)
roc_uni <- predsurv::fun_test(obj = mod_uni, train_data = train,  pred = "ROC", integrated = FALSE, noboot = 10 , test_data = test)
roc_enet <- predsurv::fun_test(obj = mod_enet, train_data = train,  pred = "ROC", integrated = FALSE, noboot = 10 , test_data = test)

##brier
brier_lasso <- predsurv::fun_test(obj = mod_lasso, train_data = train, test_data = test, pred = "Brier", integrated = FALSE)
brier_ridge <- predsurv::fun_test(obj = mod_ridge, train_data = train, data = lungdata, pred = "Brier", integrated = FALSE)
brier_uni <- predsurv::fun_test(obj = mod_uni, train_data = train, data = lungdata, pred = "Brier", integrated = FALSE)
brier_enet <- predsurv::fun_test(obj = mod_enet, train_data = train, data = lungdata, pred = "Brier", integrated = FALSE)
brier_bst<- predsurv::fun_test(obj = mod_bst, train_data = train, data = lungdata, pred = "Brier", integrated = FALSE)

### Gives attr for plot
attr(roc_uni, 'prediction.of.model') <- "Uni"
attr(roc_lasso, 'prediction.of.model') <- "Lasso"
attr(roc_enet, 'prediction.of.model') <- paste0("ENet (a = ", attr(mod_enet, 'chosen.alpha'), ")")
attr(roc_ridge, 'prediction.of.model') <- "ridge"

### Gives attr for plot
attr(brier_uni, 'prediction.of.model') <- "Uni"
attr(brier_lasso, 'prediction.of.model') <- "Lasso"
attr(brier_enet, 'prediction.of.model') <- paste0("ENet (a = ", attr(mod_enet, 'chosen.alpha'), ")")
attr(brier_ridge, 'prediction.of.model') <- "ridge"
attr(brier_bst, 'prediction.of.model') <- "Random forest"


rocplot <- roc.plot2(roc_uni, roc_lasso, roc_enet, roc_iter_enet, roc_ridge) +
  labs(subtitle = paste0("Time = ", round(quantile(test$os_months, .67),0 ) , " months")) + guides(fill = guide_legend(title = "Models", title.position = "top", keywidth = 10, keyheight = 1, col = guide_legend(ncol = 2)))

brierplot <- plot_brier2(brier_uni, brier_lasso, brier_enet, brier_bst, brier_iter_enet, brier_ridge)+  guides(fill = guide_legend(title = "Models", title.position = "top", keywidth = 10, keyheight = 1, col = guide_legend(ncol = 2)))

pdf('sim_results.pdf')
ggpubr::ggarrange(tte, km , brierplot, rocplot,
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
dev.off()


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

correlations <- cor(na.omit(lungdata[,3:23]))
# correlations
row_indic <- apply(correlations, 1, function(x) sum(x > 0.3 | x < -0.3) > 1)
correlations<- correlations[row_indic ,row_indic ]
pdf('corrplot.pdf')
corrplot::corrplot(correlations, method="square")
dev.off()


##Create train-test split
set.seed(83742)
fold = create_training_test_set(lungdata, p = 0.8, status = os_deceased)
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
attr(roc_iter_enet, 'prediction.of.model') <- paste0("Iter-ENet (a = ", attr(mod_iter_enet, 'chosen.alpha'), ")")
attr(roc_ridge, 'prediction.of.model') <- "ridge"

### Gives attr for plot
attr(brier_uni, 'prediction.of.model') <- "Uni"
attr(brier_lasso, 'prediction.of.model') <- "Lasso"
attr(brier_enet, 'prediction.of.model') <- paste0("ENet (a = ", attr(mod_enet, 'chosen.alpha'), ")")
attr(brier_iter_enet, 'prediction.of.model') <- paste0("Iter-ENet (a = ", attr(mod_iter_enet, 'chosen.alpha'), ")")
attr(brier_ridge, 'prediction.of.model') <- "ridge"
attr(brier_bst, 'prediction.of.model') <- "Random forest"


rocplot <- roc.plot2(roc_uni, roc_lasso, roc_enet, roc_iter_enet, roc_ridge) +
  labs(subtitle = paste0("Time = ", round(quantile(test$os_months, .67),0 ) , " months")) + guides(fill = guide_legend(title = "Models", title.position = "top", keywidth = 10, keyheight = 1, col = guide_legend(ncol = 2)))

brierplot <- plot_brier2(brier_uni, brier_lasso, brier_enet, brier_bst, brier_iter_enet, brier_ridge)+  guides(fill = guide_legend(title = "Models", title.position = "top", keywidth = 10, keyheight = 1, col = guide_legend(ncol = 2)))

# plot(my_trained_models[[3]])

pdf('lung_results.pdf')
ggpubr::ggarrange(tte, km , brierplot, rocplot,
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
dev.off()
