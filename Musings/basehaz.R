mod <-  survival::coxph(Surv(os_months, os_deceased)~1,
                        data = lungdata )
bh <- basehaz(mod)
coef <- mc_samp$mod_ridge$`1`

exp(-bh[,1])^(exp(coef))
