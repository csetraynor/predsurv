##Cross validation
memory.limit(5e10)
data("lungdata")
lungdata <- tibble::rownames_to_column(lungdata, var = "subject")

#Standardise gene matrix
lungdata[,4:7132] <- std_dat(lungdata[,4:7132])

#create resamples
set.seed(9666)
mc_samp <- rsample::mc_cv(lungdata, strata = "os_deceased", times = 100, prop = 3/4)


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

mc_samp$uni <- readRDS("brier_uni.RDS")
mc_samp$lasso <- readRDS("brier_lasso.RDS")
mc_samp$enet <- readRDS("brier_enet.RDS")
mc_samp$iter_enet <- readRDS("brier_iter_enet.RDS")
mc_samp$bst <- readRDS("brier_bst.RDS")


mc_samp$uni <- readRDS("roc_uni.RDS")
mc_samp$lasso <- readRDS("roc_lasso.RDS")
mc_samp$enet <- readRDS("roc_enet.RDS")
mc_samp$iter_enet <- readRDS("roc_iter_enet.RDS")
mc_samp$ridge <- readRDS("roc_ridge.RDS")

mc_samp %>%
  select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "top")

library(tidyposterior)
mc_samp <- perf_mod(mc_samp, seed = 6507, iter = 5000)

ggplot(tidy(mc_samp)) +
  theme_bw()


comparisons <- contrast_models(
  mc_samp,
  list_1 = rep("full", 3),
  list_2 = c("ph.ecog", "age", "sex"),
  seed = 4654
)
