##Cross validation
#If Windonws
memory.limit(5e10)

#Load data
library(predsurv)
data("lungdata")
lungdata <- tibble::rownames_to_column(lungdata, var = "subject")

#Standardise gene matrix
lungdata[,4:7132] <- std_dat(lungdata[,4:7132])

#create resamples
set.seed(9666)
mc_samp <- rsample::mc_cv(lungdata, strata = "os_deceased", times = 100, prop = 3/4)

##Ridge regression
mc_samp$mod_ridge <- readRDS("performance_results/mod_ridge.RDS")
mc_samp$brier_ridge <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_ridge),
                                       function(splits, mod){predsurv::fun_test(
                                         obj = mod,
                                         train_data = splits,
                                         data = lungdata,
                                         pred = "Brier"
                                       )
                                       })
saveRDS(mc_samp$brier_ridge, "brier_ridge.RDS")


###Lasso

mc_samp$mod_lasso <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Lasso")
mc_samp$roc_lasso <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_lasso),
                                    function(splits, mod){predsurv::fun_test(
                                      obj = mod,
                                      train_data = splits,
                                      data = lungdata,
                                      pred = "ROC"
                                    )
                                    })
mc_samp$brier_lasso <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_lasso),
                                    function(splits, mod){predsurv::fun_test(
                                      obj = mod,
                                      train_data = splits,
                                      data = lungdata,
                                      pred = "Brier"
                                    )
                                    })
mc_samp$deviance_lasso <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_lasso),
                                    function(splits, mod){predsurv::fun_test(
                                      obj = mod,
                                      train_data = splits,
                                      data = lungdata,
                                      pred = "Deviance"
                                    )
                                    })
mc_samp$ci_lasso <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_lasso),
                                    function(splits, mod){predsurv::fun_test(
                                      obj = mod,
                                      train_data = splits,
                                      data = lungdata,
                                      pred = "c_index"
                                    )
                                    })
saveRDS(mc_samp$mod_lasso, "mod_lasso.RDS")
saveRDS(mc_samp$roc_lasso, "roc_lasso.RDS")
saveRDS(mc_samp$brier_lasso, "brier_lasso.RDS")
mc_samp$mod_lasso <- NULL
mc_samp$roc_lasso <- NULL
mc_samp$brier_lasso <- NULL


##Univariate
mc_samp$mod_uni <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Univariate")
mc_samp$roc_uni <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_uni),
                                     function(splits, mod){predsurv::fun_test(
                                       obj = mod,
                                       train_data = splits,
                                       data = lungdata,
                                       pred = "ROC"
                                     )
                                     })
mc_samp$brier_uni <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_uni),
                                       function(splits, mod){predsurv::fun_test(
                                         obj = mod,
                                         train_data = splits,
                                         data = lungdata,
                                         pred = "Brier"
                                       )
                                       })
saveRDS(mc_samp$mod_uni, "mod_uni.RDS")
mc_samp$mod_uni <- NULL
saveRDS(mc_samp$roc_uni, "roc_uni.RDS")
saveRDS(mc_samp$brier_uni, "brier_uni.RDS")
mc_samp$roc_uni <- NULL
mc_samp$brier_uni <- NULL


##Elastic net
mc_samp$mod_enet <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Elastic net")
mc_samp$roc_enet <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_enet),
                                   function(splits, mod){predsurv::fun_test(
                                     obj = mod,
                                     train_data = splits,
                                     data = lungdata,
                                     pred = "ROC"
                                   )
                                   })
mc_samp$brier_enet <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_enet),
                                     function(splits, mod){predsurv::fun_test(
                                       obj = mod,
                                       train_data = splits,
                                       data = lungdata,
                                       pred = "Brier"
                                     )
                                     })
saveRDS(mc_samp$mod_enet, "mod_enet.RDS")
mc_samp$mod_enet <- NULL
saveRDS(mc_samp$roc_enet, "roc_enet.RDS")
saveRDS(mc_samp$brier_enet, "brier_enet.RDS")
mc_samp$roc_enet <- NULL
mc_samp$brier_enet <- NULL


##Iter-Elastic net
mc_samp$mod_iter_enet <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Elastic net",
                                    iterative = TRUE)
mc_samp$roc_iter_enet <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_iter_enet),
                                    function(splits, mod){predsurv::fun_test(
                                      obj = mod,
                                      train_data = splits,
                                      data = lungdata,
                                      pred = "ROC"
                                    )
                                    })
mc_samp$brier_iter_enet <- purrr::pmap_dbl(list(mc_samp$splits, mc_samp$mod_iter_enet),
                                      function(splits, mod){predsurv::fun_test(
                                        obj = mod,
                                        train_data = splits,
                                        data = lungdata,
                                        pred = "Brier"
                                      )
                                      })
saveRDS(mc_samp$mod_iter_enet, "mod_iter_enet.RDS")
mc_samp$mod_iter_enet <- NULL
saveRDS(mc_samp$roc_iter_enet, "roc_iter_enet.RDS")
saveRDS(mc_samp$brier_iter_enet, "brier_iter_enet.RDS")
mc_samp$roc_iter_enet <- NULL
mc_samp$brier_iter_enet <- NULL


##Random forest

mc_samp$brier_bst <- purrr::map(mc_samp$splits, predsurv::fun_score, fit = "Random forest",
                              data = lungdata, prediction = "Brier")


saveRDS(mc_samp$brier_bst, "brier_bst.RDS")
mc_samp$brier_bst <- NULL



# mc_samp$brier_uni <- purrr::map_dbl(mc_samp$splits, predsurv::fun_score, fit = "Univariate", data = lungdata, prediction = "Brier")
# mc_samp$brier_lasso <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Lasso")
# mc_samp$brier_lasso <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Elastic net")
# mc_samp$brier_iter <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Univariate", iterative = TRUE)
# mc_samp$brier_bst  <- purrr::map(mc_samp$splits, predsurv::fun_train, fit = "Random forest")

library(dplyr)
library(tidyr)
library(ggplot2)
##Plot Brier
mc_samp$uni <- readRDS("performance_results/brier_uni.RDS")
mc_samp$lasso <- readRDS("performance_results/brier_lasso.RDS")
mc_samp$enet <- readRDS("performance_results/brier_enet.RDS")
mc_samp$ridge <- readRDS("performance_results/brier_ridge.RDS")
mc_samp$iter<- readRDS("performance_results/brier_iter_enet.RDS")
mc_samp$bst <- unlist(readRDS("performance_results/brier_bst.RDS"))


brier_dens <- mc_samp %>%
  dplyr::select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "bottom")
brier_dens <- brier_dens +
  labs(x = "Integrated Brier Score",
       title = "Density of iBrier") +
  geom_vline(xintercept =  0.25 , linetype = "dotted" )

library(tidyposterior)
mc_samp_brier <- tidyposterior::perf_mod(mc_samp, seed = 6507, iter = 5000)

mbri_tab <- summary(tidy(mc_samp_brier))
mbri_tab <- as.data.frame(mbri_tab)

star = stargazer(mbri_tab, type = "latex", summary = FALSE, digits.extra = 3,digits = 3, digit.separator = ".",
                 title = "Bayesian analysis of resampling AUC")


stargazer(mbri_tab, type = "latex", summary = FALSE, digits.extra = 3,
          digits = 3, digit.separator = ".",
          title = "Bayesian analysis of resampling AUC")


posterior_brier <- ggplot(tidy(mc_samp_brier)) +
  theme_bw()+
  labs(
       title = "Posterior probability for integrated Brier Score")
posterior_brier <- posterior_brier +   labs(
  title = "Posterior probability of iBrier")

comparisons_brier <- contrast_models(
  mc_samp_brier,
  list_1 = rep("ridge", 5),
  list_2 = c("uni", "bst", "enet", "lasso", "iter"),
  seed = 4654
)

compare_brier <- ggplot(comparisons_brier, size = 0.05) +
  theme_bw()+
  labs(
    title = "Posterior probability of iBrier.",
subtitle ="Benchmark: ridge regression")

compare_brier <- compare_brier +   labs(
  title = "Posterior probability for iBrier",
  subtitle ="Benchmark: ridge regression")

#### ROC

mc_samp$uni <- readRDS("performance_results/roc_uni.RDS")
mc_samp$lasso <- readRDS("performance_results/roc_lasso.RDS")
mc_samp$enet <- readRDS("performance_results/roc_enet.RDS")
mc_samp$iter <- readRDS("performance_results/roc_iter_enet.RDS")
mc_samp$ridge <- readRDS("performance_results/roc_ridge.RDS")
mc_samp$bst <- NULL
roc_dens <- mc_samp %>%
  dplyr::select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "AUC-ROC",
      title = "Density of AUC-ROC") +
  geom_vline(xintercept =  0.5, linetype = "dotted" )+
  ylim(0,5)

mc_samp_roc <- tidyposterior::perf_mod(mc_samp, seed = 6507, iter = 5000, hetero_var = TRUE)

mcroc_tab <- summary(tidy(mc_samp_roc))
mcroc_tab <- as.data.frame(mcroc_tab)
stargazer(mcroc_tab, type = "latex", summary = FALSE,
          title = "Bayesian analysis of resampling AUC",
          digits = 3, digits.extra = 3)


posterior_roc <- ggplot(tidy(mc_samp_roc)) +
  theme_bw()+
  labs(
    title = "Posterior probability of AUC-ROC")+
  geom_hline(yintercept = 0.5, linetype = "dotted" )


comparisons_roc <- contrast_models(
  mc_samp_roc,
  list_1 = rep("ridge", 4),
  list_2 = c("uni", "enet", "lasso", "iter"),
  seed = 4654
)

compare_roc <- ggplot(comparisons_roc, size = 0.05) +
  theme_bw()+
  labs(
    title = "Posterior probability of AUC-ROC.",
subtitle = "Benchmark: ridge")



pdf('lung_mcmc.pdf')
ggpubr::ggarrange(brier_dens,  posterior_brier,
                  roc_dens, posterior_roc,
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
dev.off()

pdf('lung_mcmc_compare.pdf')
ggpubr::ggarrange(compare_brier, compare_roc,
                  labels = c("A", "B"),
                  ncol = 2, nrow = 1)
dev.off()

library(stargazer)

tab_roc <- summary(comparisons_roc, size = 0.05) %>%
  select(contrast, starts_with("pract"))
tab_roc <- as.data.frame(tab_roc)


stargazer(tab_roc, type = "latex", summary = FALSE,
          title = "Posterior distribution of AUC compared to ridge regression", digits = 3, digit.separate = 3)


tab_brier <- summary(comparisons_brier, size = 0.05) %>%
  select(contrast, starts_with("pract"))

tab_brier <- as.data.frame(tab_brier)

stargazer(tab_brier, type = "latex", summary = FALSE,
          title = "Posterior distribution of iBrier compared to ridge regression", digits = 3, digit.separate = 3)
