library(predsurv)
theme_set(theme_bw())


####Simulation study
set.seed(12318)
survdata <- surv_sim_data(N = 500, features = 500, CenRate = 0.25)

#Exploratory Analysis
tte <- plot_tte_dist(survdata, time = os_months, status = os_deceased)
km <- plot_km(survdata, time = os_months, status = os_deceased)
ggpubr::ggarrange(tte, km,
          labels = c("A", "B"),
          ncol = 2)

##Create train-test split
fold = create_training_test_set(survdata, p = 0.8, status = os_deceased)
train = fold[["train"]]
test = fold[["test"]]
head(train[,1:6]); rm(fold)

##Train Models##
my_model_list <- list("Univariate", "Ridge regression", "Lasso", "Elastic net", "Random forest", "PCR")
my_trained_models <- lapply(my_model_list, function(m) fun_train(train = train, fit = m))

plot(my_trained_models[[3]])


##Test Models
my_tested_models <- lapply(my_trained_models[c(1,2,3,4,5)], function(m) fun_test(test_data = test,obj = m, all = TRUE))

uni.roc <- attr(my_tested_models[[1]], 'roc_pred')
ridge.roc <- attr(my_tested_models[[2]], 'roc_pred')
lasso.roc <- attr(my_tested_models[[3]], 'roc_pred')
enet.roc <- attr(my_tested_models[[4]], 'roc_pred')
attr(uni.roc, 'prediction.of.model') <- "Univariate"
attr(ridge.roc, 'prediction.of.model') <- "Ridge regression"
attr(lasso.roc, 'prediction.of.model') <- "Lasso"
attr(enet.roc, 'prediction.of.model') <- "Elastic net"

# stree = fun_train(train, fit = "tree")
# plot(stree)
#Random Forest

#var_imp(bst)
rocplot <- roc.plot2(uni.roc, ridge.roc, lasso.roc, enet.roc)


uni.brier <- attr(my_tested_models[[1]], 'brier_pred')
ridge.brier <- attr(my_tested_models[[2]], 'brier_pred')
lasso.brier <- attr(my_tested_models[[3]], 'brier_pred')
enet.brier <- attr(my_tested_models[[4]], 'brier_pred')
bs.brier <- my_tested_models[[5]]
attr(uni.brier, 'prediction.of.model') <- "Univariate"
attr(ridge.brier, 'prediction.of.model') <- "Ridge regression"
attr(lasso.brier, 'prediction.of.model') <- "Lasso"
attr(enet.brier, 'prediction.of.model') <- "Elastic net"
attr(bs.brier, 'prediction.of.model') <- "Random forest"
brierplot <- plot_brier2(uni.brier, ridge.brier, lasso.brier, enet.brier, bs.brier)
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
