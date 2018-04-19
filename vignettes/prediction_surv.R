library(predsurv)
theme_set(theme_bw())

# data("survdata")
# survdata <- surv_sim_data(N = 500, features = 500, CenRate = 0.25)

data("survdata")
library(dplyr)

#Exploratory Analysis
plot_tte_dist(survdata, time = os_months, status = os_deceased)
plot_km(survdata, time = os_months, status = os_deceased)

##Create train-test split
fold = create_training_test_set(survdata, p = 0.8, status = os_deceased)
train = fold[["train"]]
test = fold[["test"]]
head(train[,1:6]); rm(fold)

##Train Models##
my_model_list <- list("Univariate", "Ridge regression", "Lasso", "Elastic net", "Random forest", "PCR")
my_trained_models <- lapply(my_model_list, function(m) fun_train(train = train, fit = m))

my_trained_models[c(4)]


##Test Models
my_tested_models <- lapply(my_trained_models[c(1,4,5)], function(m) fun_test(test_data = test,obj = m, all = TRUE))

uni.roc <- attr(my_tested_models[[1]], 'roc_pred')
enet.roc <- attr(my_tested_models[[2]], 'roc_pred')
attr(uni.roc, 'predction.of.model') <- "Univariate"
attr(enet.roc, 'predction.of.model') <- "Elastic net"

# stree = fun_train(train, fit = "tree")
# plot(stree)
#Random Forest

#var_imp(bst)
roc.plot2(uni.roc, enet.roc)


uni.brier <- attr(my_tested_models[[1]], 'brier_pred')
enet.brier <- attr(my_tested_models[[2]], 'brier_pred')
bst.brier <- my_tested_models[[3]]
attr(uni.brier, 'predction.of.model') <- "Univariate"
attr(enet.brier, 'predction.of.model') <- "Elastic net"
attr(bst.brier, 'predction.of.model') <- "Random forest"
plot_brier2(uni.brier, enet.brier, bst.brier)


##PCR
mod.pcr = fun_train(train, fit = "PCR")
pcr.brier_score = pred_error(obj = mod.pcr, fit = "PCR")
plot_brier(pcr.brier_score)

pcr.roc = pred_error(obj = mod.pcr, fit = "PCR", pred = "ROC")
tdROC::plot.tdROC(pcr.roc)

