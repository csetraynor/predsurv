set.seed(332)
x<-matrix(rnorm(1000*80),ncol=80)
y<-10+svd(x[1:60,])$v[,1]+ .1*rnorm(80)
censoring.status<- sample(c(rep(1,30),rep(0,10)))
featurenames <- paste("feature",as.character(1:1000),sep="")
x <- t(x)
colnames(x) <- featurenames
survdata<- data.frame(surv_months=y, censoring.status=censoring.status)
survdata <- cbind(survdata, x)
