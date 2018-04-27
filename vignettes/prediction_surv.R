library(predsurv)
theme_set(theme_bw())

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

pairs(survdata[,3:10])

correlations <- cor(na.omit(lungdata[,3:23]))
# correlations
row_indic <- apply(correlations, 1, function(x) sum(x > 0.3 | x < -0.3) > 1)
correlations<- correlations[row_indic ,row_indic ]
pdf('corrplot.pdf')
corrplot::corrplot(correlations, method="square")
dev.off()


##Create train-test split
fold = create_training_test_set(lungdata, p = 0.8, status = os_deceased)
train = fold[["train"]]
test = fold[["test"]]
head(train[,1:6]); rm(fold)
##Train Models##
my_model_list <- list("Univariate",  "Ridge regression", "Lasso", "Elastic net", "Random forest")

my_trained_models <- lapply(my_model_list, function(m) fun_train(train = train, fit = m))
#Iterative Method
my_trained_models[[6]] <- fun_train(train = train, fit = "Univariate", iterative = TRUE)

plot(my_trained_models[[4]])

##Test Models
my_tested_models <- lapply(my_trained_models[c(1,2,4,5,6)], function(m) fun_test(test_data = test,obj = m, all = TRUE, train_data = train))


pdf('lung_corr.pdf')
par(mfrow = c(2, 2))
plot(my_trained_models[[3]])
title("Lasso",line=2.5)
plot(my_trained_models[[4]])
title("Elastic net",line=2.5)
plot(my_trained_models[[6]])
title("Iterative elastic net",line=2.5)
corrplot::corrplot(correlations, method="square")
dev.off()


uni.roc <- attr(my_tested_models[[1]], 'roc_pred')
lasso.roc <- attr(my_tested_models[[2]], 'roc_pred')
enet.roc <- attr(my_tested_models[[3]], 'roc_pred')
iter.roc <- attr(my_tested_models[[5]], 'roc_pred')
attr(uni.roc, 'prediction.of.model') <- "Uni"
attr(lasso.roc, 'prediction.of.model') <- "Lasso"
attr(enet.roc, 'prediction.of.model') <- paste0("ENet (α = ", attr(my_trained_models[[4]], 'chosen.alpha'), ")")
attr(iter.roc, 'prediction.of.model') <- paste0("Iter-ENet (α = ", attr(my_trained_models[[6]], 'chosen.alpha'), ")")

rocplot <- roc.plot2(uni.roc, lasso.roc, enet.roc, iter.roc) +
  labs(subtitle = paste0("Time = ", round(quantile(test$os_months, .67),0 ) , " months")) + guides(fill = guide_legend(title = "Models", title.position = "top", keywidth = 10, keyheight = 1, col = guide_legend(ncol = 2)))

uni.brier <- attr(my_tested_models[[1]], 'brier_pred')
lasso.brier <- attr(my_tested_models[[2]], 'brier_pred')
enet.brier <- attr(my_tested_models[[3]], 'brier_pred')
bs.brier <- my_tested_models[[4]]
iter.brier <- attr(my_tested_models[[5]], 'brier_pred')
attr(uni.brier, 'prediction.of.model') <- "Uni"
attr(lasso.brier, 'prediction.of.model') <- "Lasso"
attr(enet.brier, 'prediction.of.model') <- paste0("Enet (α=", attr(my_trained_models[[4]], 'chosen.alpha'), ")")
attr(iter.brier, 'prediction.of.model') <- paste0("Iter-Enet (α=", attr(my_trained_models[[6]], 'chosen.alpha'), ")")
attr(bs.brier, 'prediction.of.model') <- "Random Forest"
brierplot <- plot_brier2(uni.brier, lasso.brier, enet.brier, bs.brier, iter.brier)+  guides(fill = guide_legend(title = "Models", title.position = "top", keywidth = 10, keyheight = 1, col = guide_legend(ncol = 2)))

plot(my_trained_models[[3]])

pdf('lung_results.pdf')
ggpubr::ggarrange(tte, km , brierplot, rocplot,
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
dev.off()

##Cross validation
data("lungdata")
lungdata <- tibble::rownames_to_column(lungdata, var = "subject")

# cv.uni <- fun_cv(data = lungdata, fit = "Univariate", KMC = 10)
# cv.enet <- fun_cv(data = lungdata, fit = "Elastic net", KMC = 10)
# cv.iter <- fun_cv(data = lungdata, iter= TRUE, fit = "Elastic net", KMC = 10)
# cv.bst <- fun_cv(data = lungdata, fit = "Random forest", KMC = 10)
#
# brier_est <- rbind(extract_pars_cv(cv.uni, 'brier_pred', "Univariate"),
#                extract_pars_cv(cv.enet, 'brier_pred', "Elastic net"))
#
# extract_pars_cv(cv.uni, 'roc_pred', "Univariate")
# extract_pars_cv(cv.uni, 'ci_pred', "Univariate")
# extract_pars_cv(cv.uni, 'dev_pred', "Univariate")
#
# pdf('brier_cv.pdf')
# brier_est %>%
#   ggplot(aes(x = ibrier, col = model)) +
#   geom_line(stat = "density") +
#   theme_bw() +
#   theme(legend.position = "top")
# dev.off()
#Create folds
# set.seed(9)
# folds_l <- rsample::vfold_cv(lungdata, strata = "os_deceased", v = 10)
set.seed(9666)
mc_samp <- rsample::mc_cv(lungdata, strata = "os_deceased", times = 100, prop = 3/4)

mc_samp$brier_uni <- purrr::map_dbl(mc_samp$splits, predsurv::fun_score, fit = "Univariate", data = lungdata, prediction = "Brier")
mc_samp$brier_lasso <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Lasso")
mc_samp$brier_enet <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Elastic net")
mc_samp$brier_iter <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Univariate", iterative = TRUE)
mc_samp$brier_bst  <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Random forest")

# purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Random forest")
mc_samp$roc_uni <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_uni2),
                                function(splits, mod){predsurv::fun_test(
                                  obj = mod,
                                  train_data = splits,
                                  data = lungdata,
                                  pred = "ROC"
                                )
                                })
mc_samp$roc_enet <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_enet),
                                 function(splits, mod){predsurv::fun_test(
                                   obj = mod,
                                   train_data = splits,
                                   data = lungdata,
                                   pred = "ROC"
                                 )
                                 })
mc_samp$roc_iter <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_iter),
                                 function(splits, mod){predsurv::fun_test(
                                   obj = mod,
                                   train_data = splits,
                                   data = lungdata,
                                   pred = "ROC"
                                 )
                                 })
mc_samp$roc_bst <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_bst),
                                function(splits, mod){predsurv::fun_test(
                                  obj = mod,
                                  train_data = splits,
                                  data = lungdata,
                                  pred = "ROC"
                                )
                                })
load("C:/RFactory\mc_cv.RData")
