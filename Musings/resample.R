library(tidyposterior)
library(survival)
str(lung)
lung_mod <- survreg(Surv(time, status) ~ ph.ecog + age + strata(sex), data = lung)
summary(lung_mod)
library(rsample)
set.seed(9666)
mc_samp <- mc_cv(lung, strata = "status", times = 100)

library(purrr)
cens_rate <- function(x) mean(analysis(x)$status == 1)
summary(map_dbl(mc_samp$splits, cens_rate))

three_fact <- as.formula(Surv(time, status) ~ ph.ecog + age + strata(sex))
rm_ph.ecog <- as.formula(Surv(time, status) ~           age + strata(sex))
rm_age     <- as.formula(Surv(time, status) ~ ph.ecog +       strata(sex))
rm_sex     <- as.formula(Surv(time, status) ~ ph.ecog + age              )

mod_fit <- function(x, form, ...)
  survreg(form, data = analysis(x), ...)

get_concord <- function(split, mod, ...) {
  pred_dat <- assessment(split)
  pred_dat$pred <- predict(mod, newdata = pred_dat)
  survConcordance(Surv(time, status) ~ pred, pred_dat, ...)$concordance
}


mc_samp$mod_full    <- map(mc_samp$splits, mod_fit, form = three_fact)
mc_samp$mod_ph.ecog <- map(mc_samp$splits, mod_fit, form = rm_ph.ecog)
mc_samp$mod_age     <- map(mc_samp$splits, mod_fit, form = rm_age)
mc_samp$mod_sex     <- map(mc_samp$splits, mod_fit, form = rm_sex)

mc_samp$full    <- map2_dbl(mc_samp$splits, mc_samp$mod_full, get_concord)
mc_samp$ph.ecog <- map2_dbl(mc_samp$splits, mc_samp$mod_ph.ecog, get_concord)
mc_samp$age     <- map2_dbl(mc_samp$splits, mc_samp$mod_age, get_concord)
mc_samp$sex     <- map2_dbl(mc_samp$splits, mc_samp$mod_sex, get_concord)

library(dplyr)
concord_est <- mc_samp %>%
  select(-matches("^mod"))

library(tidyr)
library(ggplot2)
concord_est %>%
  select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "top")

library(tidyposterior)
concord_est <- perf_mod(concord_est, seed = 6507, iter = 5000)

ggplot(tidy(concord_est)) +
  theme_bw()

comparisons <- contrast_models(
  concord_est,
  list_1 = rep("full", 3),
  list_2 = c("ph.ecog", "age", "sex"),
  seed = 4654
)

ggplot(comparisons, size = 0.05) +
  theme_bw()




##Baseline hazard
cox.fit.null <- coxph(Surv(time, status) ~ 1,data = as.data.frame(mc_samp$splits$`1`)  )
cox.fit1 <- coxph(three_fact, data =as.data.frame(mc_samp$splits$`1`) )
cox.fit2 <- coxph(rm_ph.ecog, data = as.data.frame(mc_samp$splits$`1`))
bhnull=basehaz(cox.fit.null, centered = FALSE)
bh1=basehaz(cox.fit1, centered = FALSE)[,1:2]
bh2=basehaz(cox.fit2, centered = FALSE)[,1:2]
plot(bhnull)
plot(bh1)
plot(bh2)


plot(bh1[,2],bh1[,1],main="Cumulative hazard function",xlab="Time",ylab="H0(t)")
