library(survival) 
library(glmnet) 
load(system.file("doc","VignetteExample.rdata",package="glmnet")) 
attach(patient.data) 

# leave the first patient for testing 
# and train glmnet on all other patients 
trainX      <-x[-1,] 
trainTime   <-time[-1] 
trainStatus <- status[-1] 

# fit Coxnet 
fit <- glmnet(trainX,Surv(trainTime,trainStatus),family="cox",alpha=0.5,maxit=10000) 

# find lambda for which dev.ratio is max 
max.dev.index     <- which.max(fit$dev.ratio) 
optimal.lambda <- fit$lambda[max.dev.index] 

# take beta for optimal lambda 
optimal.beta  <- fit$beta[,max.dev.index] 

# find non zero beta coef 
nonzero.coef <- abs(optimal.beta)>0 
selectedBeta <- optimal.beta[nonzero.coef] 

# take only covariates for which beta is not zero 
selectedTrainX   <- trainX[,nonzero.coef] 

# create coxph object with pre-defined coefficients 
coxph.model<- coxph(Surv(trainTime,trainStatus) ~selectedTrainX,init=selectedBeta,iter=0) 

# take test covariates for which beta is not zero 
selectedTestX <- x[1,nonzero.coef] 

# find survival curve for test subject 
sfit<- survfit(coxph.model,newdata=selectedTestX) 
cat("\ntime ") 
cat(sfit$time) 
cat("\nsurvival ") 
cat(sfit$surv) 


http://r.789695.n4.nabble.com/estimating-survival-times-with-glmnet-and-coxph-td4614225.html