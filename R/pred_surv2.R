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

<<<<<<< HEAD
fun_train2 <- function(hold_out, data, time = os_months, status = os_deceased, fit, subject = patient_id, lambda = 0.001, preCut = TRUE, ...){
=======
fun_train2 <- function(hold_out, data, time = os_months, status = os_deceased, fit, subject = patient_id, lambda = 0.001, ...){
>>>>>>> 80e0aa1757dcfc9309be6ff879774478dca18126

  ###### abstract dependent vars
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  subject <- dplyr::enquo(subject)
  # drop_var <- dplyr::quos(...)
  print(unlist(hold_out$id))
  ##### convert to dataset
  hold_out <- as.data.frame(hold_out)

  #### obtain training datset
  if('patient_id' %in% colnames(data)){
  train <- data[!( (data %>%
                           dplyr::select(!!subject) %>%
                           unlist) %in% (hold_out %>%
                                           dplyr::select(!!subject) %>%
                                           unlist) ),];
  #unselect subject
  train <- train %>% dplyr::select(- !!subject)
  }

  # create predictor matrix
  trainX <- train %>% dplyr::select(-!!time, -!!status)

  ###### create fold id for CV in glmnet
  set.seed(9)
  foldid <-  caret::createFolds(train %>% select(!!status) %>% unlist,
                                k = 10, list = FALSE)

<<<<<<< HEAD
=======
  ##### create penalty.factor vector
  p.fac = rep(1, ncol(trainX))
  p.fac[match(c("npi","age_std"), colnames(trainX))] = 0
>>>>>>> 80e0aa1757dcfc9309be6ff879774478dca18126
  ###### Fit models

  #Lasso
  if(fit == "Lasso" ){
    #prepare data for glmlasso
    x <- as.matrix(trainX)
    y <- as.matrix(train %>%
                     dplyr::select(time = !!time, status = !!status), ncol = 2)
    #register for parallelisation cv
<<<<<<< HEAD

    # doParallel::registerDoParallel(cores=parallel::detectCores() - 1)
    # #fit glmlasso
    # doMC::registerDoMC(cores=4)
    mod <-  glmnet::cv.glmnet(x, y, family = "cox", grouped = TRUE, lambda.min.ratio = lambda, foldid = foldid, parallel = FALSE, penalty.factor = p.fac)
    print("Done!")
=======
    # doParallel::registerDoParallel(cores=parallel::detectCores() - 1)
    #fit glmlasso

    mod <-  glmnet::cv.glmnet(x, y, family = "cox", grouped = TRUE, lambda.min.ratio = lambda, foldid = foldid, parallel = FALSE, penalty.factor = p.fac)

    #Stop Parallel
    # doParallel::stopImplicitCluster()
>>>>>>> 80e0aa1757dcfc9309be6ff879774478dca18126

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

    mod <-  glmnet::cv.glmnet(x, y, family = "cox", alpha = 0 ,grouped = TRUE, lambda.min.ratio = 0.001, foldid = foldid, parallel = FALSE, penalty.factor = p.fac)
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
<<<<<<< HEAD
    ##### For computational burden apply preUni
      #create formula for apply loop
      reg <- function(indep_var,dep_var,data_source) {
        formula <- as.formula(paste(dep_var," ~ ", indep_var))
        res     <- survival::coxph(formula, data = data_source)
        summary(res)
      }
=======
    #create a sequence of alpha
    alphaList <-  (1:19) * 0.05
    #prepare for glmnet
      x <- as.matrix(trainX)
    y <- as.matrix(train %>%
                     dplyr::select(time = !!time, status = !!status), ncol = 2)


    #do crossvalidation to find optimal alpha and lambda
    elasticnet <-  lapply(alphaList, function(a){
      glmnet::cv.glmnet(x, y, family = "cox", grouped = TRUE , alpha = a, lambda.min.ratio = 0.001, foldid = foldid, parallel = FALSE, penalty.factor = p.fac)});
    #extract optimal lambda
    cvm <- sapply(seq_along(alphaList), function(i) min(elasticnet[[i]]$cvm ) );

    a <- alphaList[match(min(cvm), cvm)];
    mod <- elasticnet[[match(min(cvm), cvm)]]
>>>>>>> 80e0aa1757dcfc9309be6ff879774478dca18126

      #### fit univariate cox models with each covariate
      mod1 <- lapply(colnames(trainX %>%
                                dplyr::select(-npi, -age_std)),
                     FUN = reg, dep_var = "Surv(time, status)",
                     data_source = train %>% dplyr::mutate(
        time   = !!time,
        status = !!status) )
      #soft correction
      mod2 <- lapply(mod1,
                     function(p)
                       p[p$logtest["pvalue"] < 0.01]$coefficients)
      ### drop covariates that are not "statistically significant"
      mod2[sapply(mod2, is.null)] <- NULL
      features <- sapply(mod2, function(p) rownames(p))
      coefficients <- sapply(mod2, function(p) p[,1])
      mod <- data.frame(coef = coefficients)
      rownames(mod) <- features
      #prepare for enet
      x <- cbind(trainX %>% dplyr::select(rownames(mod)),
                           trainX %>% dplyr::select(npi, age_std))
      ##### create penalty.factor vector
      p.fac = rep(1, ncol(x))
      p.fac[match(c("npi","age_std"), colnames(x))] = 0
      #### prepare matrix
      x <- as.matrix(x)
      y <- as.matrix(train %>%
                       dplyr::select(time = !!time,
                                     status = !!status), ncol = 2)
      ##################################################
      #Elastic net

      ### Apply elastic net with a=0.8 closer to Lasso
      mod <-  glmnet::cv.glmnet(x, y, family = "cox",
        grouped = TRUE, lambda.min.ratio = lambda, alpha = 0.8,
        foldid = foldid, parallel = FALSE, penalty.factor = p.fac)
      print("Done!")

      # find optimised lambda
      optimal.coef <- as.matrix(coef(mod, s = "lambda.min"))
      optimal.coef <- as.data.frame(optimal.coef)
      colnames(optimal.coef) <- "mod"
      optimal.coef <- tibble::rownames_to_column(optimal.coef, var = "gene")
      optimal.coef <-  optimal.coef[optimal.coef$mod != 0,]
      mod <- data.frame(coef = optimal.coef$mod)
      rownames(mod) <- optimal.coef$gene
    }
   if(fit == "Random forest"){
    train_tree <- as.data.frame(train %>%
                                  dplyr::rename(
                                    time   = !!time,
                                    status = !!status) )
    forest_control <- party::cforest_unbiased()

    mod <- pec::pecCforest(form,
                           data = train_tree,
                           controls = forest_control)
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
fun_test2 <- function(obj, train_data = train, test_data = NA,  data = NA, subject = subject, time = os_months, status = os_deceased, event_type = 1,  pred = "Brier", adapted = "default", all = FALSE, integrated = TRUE, noboot = 0, mc = FALSE, ...){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)
  subject <- dplyr::enquo(subject)

  # print(train_data$id)
  #Preprocess data
  #Get fit name
  fit <- attr(obj, "fit.model")
  #transform splits to train data
  train_data <- as.data.frame(train_data)
  #get test data from data source
  if('subject' %in% colnames(train)){
    test_data <- data[!( (data %>%
                            dplyr::select(!!subject) %>%
                            unlist) %in% (train_data %>%
                                            dplyr::select(!!subject) %>%
                                            unlist) ),];
    #unselect subject
    train_data <- train_data %>% dplyr::select(- !!subject)
    test_data <- test_data %>% dplyr::select(- !!subject)
  }
  # create coxph object with pre-defined coefficients for different models
  #currently not supported stepwise, pls, pcr
  if(fit == "Stepwise"){
    # take optimal beta from model object
    optimal.beta <- coef(obj)
    # find non zero beta coef
    nonzero.coef <- abs(optimal.beta)>0 & !is.na(optimal.beta)
    selectedBeta <- optimal.beta[nonzero.coef]
  }
  #this is the most important bit, gets the fitted model and extracts the selected features
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
  #Now prepares for getting performance measures, random forest at the bottom
  if(fit != "Random forest"){
    ##### No variables selected complete "shrinkage" tipycal in Lasso
    if(length(selectedBeta) == 0 ){
      mod <-  survival::coxph(Surv(time, status)~1, iter = 0,
                              data = train_data %>%
                                dplyr::mutate(time = !!time, status = !!status))
      # Create Test vars
      testX <- test_data %>% dplyr::select(-!!time, -!!status)
      ndata <- test_data %>% dplyr::mutate(time = !!time, status =  !!status)
    }
    #Create training and test set and formula
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
        ##### Fit multicovariate Cox model
        mod <-  survival::coxph(form , data = traincoxphdata %>% dplyr::mutate(
          time = !!time,
          status = !!status))
      }
      if(fit == "Lasso" |  fit == "Elastic net"){

        #### Fit Cox model with inits and iter 0
        mod <-  survival::coxph(form , data = traincoxphdata %>% dplyr::mutate(
          time = !!time,
          status = !!status), init = as.vector(unlist(obj)), iter = 0,
          control = coxph.control(iter.max = 0, eps= 10e9),
          singular.ok = TRUE)
        # print(mod)
      }
      #stack overflow for paste
      if(fit == "Ridge regression"){
        mod <-  survival::coxph(Surv(time = train_data %>%
                                       dplyr::select(!!time) %>%
                                       unlist,
                                     event = train_data %>%
                                       dplyr::select(!!status) %>%
                                       unlist)~.,
                                init = as.vector(unlist(obj)), iter = 0,
                                data = train_data %>% dplyr::select(- !!time, - !!status))
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
    ######################################################
    #### Prediction Error
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
                          tau = quantile(test_data %>% dplyr::select(time = !!time)%>% unlist, .73), nboot = noboot, alpha = 0.05, n.grid = 1000,  type = "uniform"
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

