### Possible new method ###
## Not run:
n <- 200
p <- 100
beta <- c(rep(1,10),rep(0,p-10))
x <- matrix(rnorm(n*p),n,p)
real.time <- -(log(runif(n)))/(10*exp(drop(x %*% beta)))
cens.time <- rexp(n,rate=1/10)
status <- ifelse(real.time <= cens.time,1,0)
time <- ifelse(real.time <= cens.time,real.time,cens.time)
test = data.frame(time = time, status = status)

######## Method now###
set.seed(332)
x<-matrix(rnorm(1000*80),ncol=80)
y<-10+svd(x)$v[,1]+ .1*rnorm(80)
censoring.status<- sample(c(rep(1,60),rep(0,20)))
featurenames <- paste("feature",as.character(1:1000),sep="")
x <- t(x)
colnames(x) <- featurenames
survdata<- data.frame(surv_months=y, censoring.status=censoring.status)
survdata <- cbind(survdata, x)

surv_sim_data <- function(alpha = 0.8, mu = -3, n = 40, features = 1000) {

  x<-matrix(rnorm(features*n),ncol= n)

  survdata <- data.frame(surv_months = rweibull(n = n, 0.8, exp(-(0.8 + svd(x[1:60,] )$v[,1] + + .1*rnorm(n))/alpha)),
                     censor_months = rexp(n = n, rate = 1/100),
                     stringsAsFactors = F
  ) %>%
    dplyr::mutate(os_status = ifelse(surv_months < censor_months,
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months,
                       surv_months, censor_months
    )
    )
  x <- t(x)
  colnames(x) <- featurenames

  survdata <- cbind(survdata, x)

  return(survdata)
}
