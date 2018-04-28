#' Fit Model to Training data
#'
#' Train data to specific feature selection method
#'
#' @param
#' train : train data \cr
#' time : time variable \cr
#' status : status variable \cr
#' method: forward, backward or both \cr
#' fit : fit model: "Stepwise", "Random forest"
#' @return a coxph forward Stepwise fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import partykit
#' @import grid
#' @import libcoin
#' @import mvtnorm
#' @import rpart
#' @import doMC

fun_train <- function(train, time = os_months, status = os_deceased, fit, penalty = "BIC" , iterative = FALSE){

  ###### abstract dependent vars
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  ##### create training dataset
  train <- as.data.frame(train)
  if('subject' %in% colnames(train)){
    train <- train %>% dplyr::select(-subject)
  }
  # take only covariates for which beta is not zero
  trainX <- train %>% dplyr::select(-!!time, -!!status)
  traincoxphdata <- cbind(train %>% dplyr::select(!!time, !!status), trainX)

  ##### create formula
   form <- as.formula(paste("Surv(time, status)", paste("~", paste(names(trainX), collapse = " + "), sep = " " )))

   ###### create fold id for CV in glmnet
   set.seed(9)
   foldid <-  caret::createFolds(train %>% select(!!status) %>% unlist,
                                 k = 10, list = FALSE)
   ###### Fit models
   #### Univariate
  if(fit == "Univariate" | iterative){
    #create formula for apply loop
    reg <- function(indep_var,dep_var,data_source) {
      formula <- as.formula(paste(dep_var," ~ ", indep_var))
      res     <- survival::coxph(formula, data = data_source)
      summary(res)
    }

    #### fit univariate cox models with each covariate
    mod1 <- lapply(colnames(trainX), FUN = reg, dep_var = "Surv(time, status)", data_source = traincoxphdata %>% dplyr::mutate(
      time   = !!time,
      status = !!status) )
    if(iterative){
      #soft correction
      mod2 <- lapply(mod1, function(p) p[p$logtest["pvalue"] < 0.001]$coefficients)
    }else{
      #Bonferroni correction
      mod2 <- lapply(mod1, function(p) p[p$logtest["pvalue"] < 0.0005]$coefficients)
    }
    ### drop covariates that are not "statistically significant"
    mod2[sapply(mod2, is.null)] <- NULL
    features <- sapply(mod2, function(p) rownames(p))
    coefficients <- sapply(mod2, function(p) p[,1])
    mod <- data.frame(coef = coefficients)
    rownames(mod) <- features
    # print(mod)

    if(iterative & length(mod) > 1){
      # change to elastic net
      fit <-  "Elastic net"
    }
  }
   #Stepwise fit, currently not recommended computationally too expensive
  if(fit == "Stepwise"){
    # create coxph object
    fit1 <- suppressWarnings(survival::coxph(Surv(time, status)~1,
                            data = traincoxphdata %>% dplyr::mutate(
                              time   = !!time,
                              status = !!status) ))

    testBeta <- as.formula(paste("~", paste(names(trainX), collapse = " + "), sep = " " ) )
    if(penalty == "AIC"){
      k = 2
    }else{
      if(penalty == "BIC"){
        k = log(nrow(trainX))
      }
    }
    #forward Stepwise selection
    mod <- suppressWarnings(MASS::stepAIC(fit1,
            scope = list(upper = testBeta, lower = ~1) ,
            trace = FALSE,
            direction = "forward",
            k = k,
            steps =  nrow(train)-1 ))
  }
   #Lasso
  if(fit == "Lasso" | fit == "Adaptive Lasso"){
    #prepare data for glmlasso
    x <- as.matrix(trainX)
    y <- as.matrix(train %>%
                     dplyr::select(time = !!time, status = !!status), ncol = 2)
    #register for parallelisation cv
    doMC::registerDoMC(cores=parallel::detectCores() - 1)
    #fit glmlasso
    mod <-  glmnet::cv.glmnet(x, y, family = "cox", grouped = TRUE, lambda.min.ratio = 0.001, foldid = foldid, parallel = TRUE)

    # find optimised lambda
    optimal.coef <- as.matrix(coef(mod, s = "lambda.min"))
    optimal.coef <- as.data.frame(optimal.coef)
    colnames(optimal.coef) <- "mod"
    optimal.coef <- tibble::rownames_to_column(optimal.coef, var = "gene")
    optimal.coef <-  optimal.coef[optimal.coef$mod != 0,]
    mod <- data.frame(coef = optimal.coef$mod)
    rownames(mod) <- optimal.coef$gene
  }

  if(fit == "Ridge regression"){
    #prepare for ridge regression
     x <- as.matrix(trainX)
     y <- as.matrix(train %>%
                      dplyr::select(time = !!time, status = !!status), ncol = 2)
     #register do Parallel
     doMC::registerDoMC(cores=parallel::detectCores() - 1)
     mod <-  glmnet::cv.glmnet(x, y, family = "cox", alpha = 0 ,grouped = TRUE, lambda.min.ratio = 0.001, foldid = foldid, parallel = TRUE) #alpha 0 ridge penalty
     # find optimised lambda
     optimal.coef <- as.matrix(coef(mod, s = "lambda.min"))
     optimal.coef <- as.data.frame(optimal.coef)
     colnames(optimal.coef) <- "mod"
     optimal.coef <- tibble::rownames_to_column(optimal.coef, var = "gene")

     optimal.coef <-  optimal.coef[optimal.coef$mod != 0,]
     mod <- data.frame(coef = optimal.coef$mod)
     rownames(mod) <- optimal.coef$gene
   }
   if(fit == "Elastic net"){
     #create a sequence of alpha
     alphaList <-  (1:19) * 0.05
     if(iterative){
       # print(mod)
       # obtain the subset of covariates from univariate
       x <- as.matrix(trainX %>% select(rownames(mod)))
     }else{
       #prepare for glmnet
       x <- as.matrix(trainX)
     }
    y <- as.matrix(train %>%
                        dplyr::select(time = !!time, status = !!status), ncol = 2)
    #register for parallelisation
    doMC::registerDoMC(cores=parallel::detectCores()-1)
    #do crossvalidation to find optimal alpha and lambda
       elasticnet <-  lapply(alphaList, function(a){
         glmnet::cv.glmnet(x, y, family = "cox", grouped = TRUE , alpha = a, lambda.min.ratio = 0.001, foldid = foldid, parallel = TRUE)});
       #extract optimal lambda
       cvm <- sapply(seq_along(alphaList), function(i) min(elasticnet[[i]]$cvm ) );
       doParallel::stopImplicitCluster()
       a <- alphaList[match(min(cvm), cvm)];
       mod <- elasticnet[[match(min(cvm), cvm)]]

       # find lambda for which dev.ratio is max
       optimal.coef <- as.matrix(coef(mod, s = "lambda.min"))
       optimal.coef <- as.data.frame(optimal.coef)
       colnames(optimal.coef) <- "mod"
       optimal.coef <- tibble::rownames_to_column(optimal.coef, var = "gene")

       optimal.coef <-  optimal.coef[optimal.coef$mod != 0,]
       mod <- data.frame(coef = optimal.coef$mod)
       rownames(mod) <- optimal.coef$gene
   }
  if(fit == "Tree"){
    ### with weight-dependent log-rank scores
    ### log-rank trafo for observations in this node only (= weights > 0)
    survtree_control <-
      party::ctree_control(testtype = "Univariate", maxdepth = 2)
    mod <- party::ctree(form,
                      data = traincoxphdata %>%
                      dplyr::mutate(
                          time   = !!time, status = !!status),
                      control = survtree_control )
  }
  if(fit == "Random forest"){
    train_tree <- as.data.frame(traincoxphdata %>%
                    dplyr::mutate(
                      time   = !!time,
                      status = !!status) )
    forest_control <- party::cforest_unbiased()

    mod <- pec::pecCforest(form,
                           data = train_tree,
                          controls = forest_control)
  }
   #currently not recommended either PCR or PLS because predictions are "tricky"
  if(fit == "PLS"){
    X_train <- apply((as.matrix(trainX)), FUN="as.numeric",MARGIN = 2)
    Y_train <- train %>% dplyr::select(time = !!time) %>% unlist %>% as.numeric
    C_train <- train %>% dplyr::select(status = !!status) %>% unlist %>% as.integer

    mod <- plsRcox::plsRcoxmodel.default(X_train, time = Y_train,event = C_train, nt=5)
  }
   if(fit == "PCR"){
     X_train <- apply((as.matrix(trainX)),FUN="as.numeric",MARGIN=2)
     X_train <- t(X_train)
     Y_train <- train %>% dplyr::select(time = !!time) %>% unlist
     C_train <- train %>% dplyr::select(status = !!status) %>% unlist
     featurenames <- colnames(trainX)
     data<-list(x=X_train,y=Y_train, censoring.status=C_train, featurenames=featurenames)
     mod <- superpc::superpc.train(data, type="survival")
   }
   attr(mod, 'fit.model') <- fit
  return(mod)
}
#' Test model
#'
#' This function takes a survival fit and following the Cox model (or random forest), estimates the hazard for individual from the model coefficient, while the baseline hazard is estimated with the Nelson-Aalen method. Gives prediction error measures Brier, c index, loglik, auc
#'
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' pred : prediction error Brier, ROC or C-Index
#' adapted : in Lasso allows to swap to "adapated Lasso", otherwise will take usual Lasso \cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
fun_test <- function(obj, train_data = train, data = lungdata, subject = subject, time = os_months, status = os_deceased, event_type = 1,  pred = "Brier", adapted = "default", all = FALSE, integrated = TRUE, ...){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  subject <- dplyr::enquo(subject)

  # print(train_data$id)
  #Preprocess data
  fit <- attr(obj, "fit.model")
  train_data <- as.data.frame(train_data)
  test_data <- data[!( (data %>%
                        dplyr::select(!!subject) %>%
                        unlist) %in% (train_data %>%
                        dplyr::select(!!subject) %>%
                        unlist) ),];
  train_data <- train_data %>% dplyr::select(- !!subject)
  test_data <- test_data %>% dplyr::select(- !!subject)
  # create coxph object with pre-defined coefficients
    if(fit == "Stepwise"){
      # take optimal beta from model object
      optimal.beta <- coef(obj)
      # find non zero beta coef
      nonzero.coef <- abs(optimal.beta)>0 & !is.na(optimal.beta)
      selectedBeta <- optimal.beta[nonzero.coef]
    }

    if(fit == "Univariate" | fit == "Lasso" | fit == "Adaptive Lasso" | fit == "Ridge regression" | fit == "Elastic net" ){
      # find lambda for which dev.ratio is max
      selectedBeta <- rownames(obj)
    }
    if(fit == "PLS" | fit == "PCR"){
      mod <- obj;
      if(fit == "PCR"){
        selectedBeta <- mod$feature.scores
      }
      if(fit == "PLS"){
        selectedBeta <- mod$Coeffs
        names(selectedBeta) <- rownames(selectedBeta)
      }
    }

  if(fit != "Random forest"){
    ##### No variables selected complete "shrinkage" tipycal in Lasso
    if(length(selectedBeta) == 0 ){
      mod <-  survival::coxph(Surv(time, status)~1, iter = 0,
                              data = train_data %>%
                                dplyr::mutate(time = !!time, status = !!status))
      # Create Test vars
      testX <- test_data %>% dplyr::select(-!!time, -!!status)
      ndata <- test_data %>% dplyr::mutate(time = !!time, status =  !!status)
      # if(pred == "Brier"){
      #   out <- 0.25
      # }
      # if(pred == "ROC"){
      #   out <- 0.5
      # }
      # if(pred == "c-index"){
      #   out <- 0.5
      # }
      # if(pred == "Deviance"){
      #   out <- -2*(mod$loglik)
      # }
    }
    #Start normal prediction
    if(length(selectedBeta) > 0 ){
      # Create train X take only covariates for which beta is not zero
      trainX <- train_data %>% dplyr::select(-!!time, -!!status)
      selectedTrainX   <- as.data.frame(trainX[,colnames(trainX) %in% selectedBeta])
      # print(selectedBeta)
      colnames(selectedTrainX) <- selectedBeta
      traincoxphdata <- cbind(train_data %>% dplyr::select(!!time, !!status), selectedTrainX)
      #Create Cox formula
      testBeta <- paste(selectedBeta, collapse = " + ")

      form <- as.formula(paste("Surv(time, status)", " ~ ", testBeta))
      #Create Cox Model
      if(fit == "Univariate" | fit == "Stepwise" | fit == "Adaptive Lasso"){
        #Fit model
        mod <-  survival::coxph(form , data = traincoxphdata %>% dplyr::mutate(
          time = !!time,
          status = !!status))
      }
      if(fit == "Lasso" |  fit == "Elastic net" | fit == "Ridge regression"){
        #Fit model
        mod <-  survival::coxph(form , data = traincoxphdata %>% dplyr::mutate(
          time = !!time,
          status = !!status), init = as.vector(unlist(obj)), iter = 0,
          control = coxph.control(iter.max = 0))
        # print(mod)
      }
      # Create Test vars
      testX <- test_data %>% dplyr::select(-!!time, -!!status)
      selectedTestX <- as.data.frame(testX[,colnames(testX) %in% selectedBeta])
      colnames(selectedTestX) <- selectedBeta
      ndata <- cbind(test_data %>% dplyr::select(!!time, !!status) , selectedTestX)
    }
    #Create grid of equidistant time points for testing
    timepoints <-  seq(0, max(test_data %>%
                                dplyr::filter(!!status == event_type) %>%
                                dplyr::select(!!time) %>% unlist
    ) , length.out = 100L)
    ####Prediction Error
    if(pred == "Brier" | all){
      if(length(selectedBeta) == 0){
        probs <- matrix(runif(nrow(test_data)*100), ncol = 100)
        brier <- pec::pec(probs, Surv(time, status) ~ 1,
                          data = test_data %>%
                            dplyr::mutate( time   = !!time, status = !!status),
                          maxtime = max(timepoints),
                          exact = F,
                          exactness = 99L)
        out <- brier
        if(integrated){
          out <- pec::crps(brier, models = "Reference")[1]
        }
      }else{
        #Calculate probs
        probs <- pec::predictSurvProb(mod,
                                      newdata = ndata, times = timepoints)
        #Calculate brier score
        brier <- pec::pec(probs, Surv(time, status) ~ 1,
                          data = ndata %>%
                            dplyr::mutate( time = !!time, status = !!status),
                          maxtime = max(timepoints),
                          exact = F,
                          exactness = 99L)
        out <- brier
        if(integrated){
          out <- predsurv::fun_ibrier_score(out)
        }
      }
    }
    if(pred == "ROC" | all ){
      if(length(selectedBeta) == 0){
        probs <- rep(0, nrow(testX))
        names(probs) <- rownames(testX)
      }else{
        if(fit == "Ridge regression"){
          probs <- as.matrix(testX, ncol = length(selectedBeta)) %*% as.vector(unlist(obj))
        }else{
          probs <- predict(mod, newdata = testX, type = "lp")
        }
      }
      roc <- tdROC::tdROC(X = probs,
                          Y = test_data %>% dplyr::select(time = !!time)%>% unlist,
                          delta = test_data %>% dplyr::select(status = !!status)%>% unlist,
                          tau = quantile(test_data %>% dplyr::select(time = !!time)%>% unlist, .73), nboot = 0, alpha = 0.05, n.grid = 1000,  type = "uniform"
      )
      out <- roc
      if(integrated){
        out <- out$AUC[1] %>% unlist
      }
    }
    if(pred  == "Deviance" | all  ){
      logl <- -2*(mod$loglik)
      out <- logl
      if(integrated){
        out <- out[2]
      }
    }
    if(pred == "c_index" | all ){
      if(length(selectedBeta) == 0){
        out <- 0.5
      }else{
        ###Create your survival estimates
        ci <-  pec::cindex(mod, formula = Surv(time, status) ~ 1,
                           data = ndata %>% dplyr::mutate(
                             time   = !!time, status = !!status))
        out <- ci

        if(integrated){
          out <-  out$AppCindex$coxph
        }
      }
    }
  }else{
    if(fit == "Random forest"){
      # Create Test vars
      testX <- test_data %>% dplyr::select(-!!time, -!!status)
      ndata <- cbind(test_data %>% dplyr::select(!!time, !!status), testX)
      #Create grid of equidistant time points for testing
      timepoints <-  seq(0,
                         max(test_data %>%
                               dplyr::filter(!!status == event_type) %>%
                               dplyr::select(!!time) %>% unlist), length.out = 100L)
      #Calculate probs
      mod <- obj;
      probs <- pec::predictSurvProb(mod,
                                    newdata = ndata, times = timepoints)
      #Calculate brier score
      out <- pec::pec(probs, Surv(time, status) ~ 1,
                      data = ndata %>% dplyr::mutate( time   = !!time, status = !!status),
                      maxtime = max(timepoints),
                      exact = F,
                      exactness = 99L)
      if(integrated){
        out <- predsurv::fun_ibrier_score(out)
      }
    }else{
      print("Not accepted fit")
    }
  }

  if(all & fit != "Random forest"){
    attr(out, 'brier_pred') <- brier
    attr(out, 'ci_pred') <- ci
    attr(out, 'roc_pred') <- roc
    attr(out, 'dev_pred') <- logl
  }
  # attr(out, 'number.of.individuals.cohort') <- nrow(train_data) + nrow(test_data)
  # attr(mod, 'predction.of.model') <- fit
  return(out)
}


#' Integrated Brier Score
#'
#' Predict integration Brier Score in a new dataset.
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
fun_ibrier_score <- function(brier){

  ibrier <- pec::crps(brier, models = "matrix")[1]

  return(ibrier)

}


#' Fit and measure measures of prognostic performance
#' Prediction of various scores
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
fun_score <- function(train, fit, data, prediction){

  obj <- predsurv::fun_train( fit = fit, train = train)
  out <- predsurv::fun_test(obj = obj, train_data = train, data = lungdata, pred = prediction)

  return(out)

}




#' Var Importance
#'
#' Experimental var importance for survival trees
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!

var_imp <- function(obj){

    plot(party::varimp(bst, conditional = TRUE)) #obtain variable importance
}
#' Var importance
#'
#' Plot Fitted Individual Tree
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!

var_importance <- function(forest = bst, type = "Permutation"){
  randomForest::importance(bst, type = 2)
}





#' Plot Cox Model Prediction
#'
#' This function takes a survival fit and plots, right now only works with factor
#'
#' @param
#' data : survival dataset \cr
#' ... : covariates to evaluate
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @importFrom rlang !!!
#' @import prodlim
#' @examples http://staff.pubhealth.ku.dk/~tag/Teaching/share/R-tutorials/SurvivalAnalysis.html

plot_cox_pred <- function(train_data=train, time = os_months, status = os_deceased,  ...){

  my_var <- dplyr::quos(...)
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  X <- train_data %>% dplyr::select(!!!my_var)

  form <- as.formula(paste("Hist(time, status)", paste("~", paste( colnames(X), collapse = " + "), sep = " " )))

  print(form)

  #library(pec)
  km.grade <- prodlim::prodlim(form,
    data = train_data %>% dplyr::mutate(
    time = !!time,
    status = !!status))

  cox.grade <- survival::coxph(form,
    data = train_data %>% dplyr::mutate(
    time = !!time,
    status = !!status) ,x=TRUE,y=TRUE)

  newdata <- data.frame(feature270 = c(1,0,-1),
                        feature350 = c(1,0,-1)
                        )
  ## first show Kaplan-Meier without confidence limits
  plot(km.grade, lty=1, lwd=3,
       col=c("darkgreen","darkorange","red"), confint=FALSE)
  ## now add survival estimates based on Cox regression
  pec::plotPredictSurvProb(cox.grade, lty=2,
                      col=c("darkgreen","darkorange","red"),
                      add=TRUE, sort(unique(train_data %>% dplyr::select(!!time) %>% unlist)),
                      newdata=newdata)
  mtext("Comparison of Univariate  Cox regression and stratified Kaplan-Meier")
  }

