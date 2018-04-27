library(dplyr)
alphaList <-  (1:19) * 0.05
x <- as.matrix(lungdata %>% dplyr::select(-os_months, -os_deceased, -subject), ncol = 7129
               )

y <- as.matrix(lungdata %>%
                 dplyr::select(time = os_months, status = os_deceased), ncol = 2)
#register do Parallel
doMC::registerDoMC(cores=4)
elasticnet <-  lapply(alphaList, function(a){
glmnet::cv.glmnet(x, y, family = "cox", grouped = TRUE , alpha = a, lambda.min.ratio = 0.001,  parallel = TRUE)});

cvm <- sapply(seq_along(alphaList), function(i) min(elasticnet[[i]]$cvm ) );

a <- alphaList[match(min(cvm), cvm)];
mod <- elasticnet[[match(min(cvm), cvm)]]
# find lambda for which dev.ratio is max
optimal.coef <- as.matrix(coef(mod, s = "lambda.min") )
mod <-  optimal.coef[optimal.coef != 0,]
mod <- t(as.data.frame(mod))
names(mod) <- colnames(mod)
mod

# find lambda for which dev.ratio is max
selectedBeta <- names(mod)


selectedTrainX   <- x[,colnames(x) %in% selectedBeta]
colnames(selectedTrainX) <- names(selectedBeta)
traincoxphdata <- cbind(train_data %>% dplyr::select(!!time, !!status), selectedTrainX)

lp <- selectedTrainX %*% as.vector(mod)
