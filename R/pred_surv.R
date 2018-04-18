#' Fit Forward Selection
#'
#' Fit the Forward Selection method
#'
#' @param
#' train : train data \cr
#' time : time variable \cr
#' status : status variable \cr
#' method: forward, backward or both \cr
#' fit : fit model: "stepwise", "random forest"
#' @return a coxph forward stepwise fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import partykit
#' @import grid
#' @import libcoin
#' @import mvtnorm
#' @import rpart

fun_fit <- function(train, time = os_months, status = os_deceased, method = "forward", fit, penalty = "BIC" ){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)


  # take only covariates for which beta is not zero
  trainX <- train %>% dplyr::select(-!!time, -!!status)

  traincoxphdata <- cbind(train %>% dplyr::select(!!time, !!status), trainX)

  #create formula


   form <- as.formula(paste("Surv(time, status)", paste("~", paste(names(trainX), collapse = " + "), sep = " " )))

  if(fit == "Univariate"){
    reg <- function(indep_var,dep_var,data_source) {
      formula <- as.formula(paste(dep_var," ~ ", indep_var))
      res     <- survival::coxph(formula, data = data_source)
      summary(res)
    }
  mod1 <- lapply(colnames(trainX), FUN = reg, dep_var = "Surv(time, status)", data_source = traincoxphdata %>% dplyr::mutate(
    time   = !!time,
    status = !!status) )
  mod2 <- lapply(mod1, function(p) p[p$logtest["pvalue"] < 0.01]$coefficients)
  mod2[sapply(mod2, is.null)] <- NULL
  features <- sapply(mod2, function(p) rownames(p))
  coefficients <- sapply(mod2, function(p) p[,1])
  mod <- t(data.frame(coef = coefficients))
  names(mod) <- features
  }

  if(fit == "stepwise"){
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

    #forward stepwise selection
    mod <- suppressWarnings(MASS::stepAIC(fit1,
            scope = list(upper = testBeta, lower = ~1) ,
            trace = FALSE,
            direction = method,
            k = k,
            steps =  nrow(train)-1 ))
  }

  if(fit == "lasso"){
    x <- as.matrix(trainX)
    y <- as.matrix(train %>%
                     dplyr::select(time = !!time, status = !!status), ncol = 2)
    mod <-  glmnet::cv.glmnet(x, y, family = "cox")
  }
   if(fit == "ridge"){
     x <- as.matrix(trainX)
     y <- as.matrix(train %>%
                      dplyr::select(time = !!time, status = !!status), ncol = 2)
     mod <-  glmnet::cv.glmnet(x, y, family = "cox", alpha = 0) #alpha 0 ridge penalty
   }

  if(fit == "tree"){
    ### with weight-dependent log-rank scores
    ### log-rank trafo for observations in this node only (= weights > 0)
    h <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
      if (is.null(weights)) weights <- rep(1, NROW(y))
      s <- coin::logrank_trafo(y[weights > 0,,drop = FALSE])
      r <- rep(0, length(weights))
      r[weights > 0] <- s
      list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
    }

    survtree_control <-
      party::ctree_control(testtype = "Univariate", maxdepth = 2)
    mod <- party::ctree(form,
                      data = traincoxphdata %>%
                      dplyr::mutate(
                          time   = !!time, status = !!status),
                      ytrafo = h,
                      control = survtree_control )
  }
  if(fit == "random forest"){
    train_tree <- as.data.frame(traincoxphdata %>%
                    dplyr::mutate(
                      time   = !!time,
                      status = !!status) )
    forest_control <- party::cforest_unbiased(n = 500)

    mod <- pec::pecCforest(form,
                           data = train_tree,
                          controls = forest_control)
  }

  if(fit == "PLS"){
    X_train <- apply((as.matrix(trainX)),FUN="as.numeric",MARGIN=2)
    Y_train <- train %>% dplyr::select(time = !!time) %>% unlist
    C_train <- train %>% dplyr::select(status = !!status) %>% unlist

    mod <- plsRcox::coxpls(X_train,Y_train,C_train,nt=6)
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

  return(mod)
}
#' Cox Model Prediction
#'
#' This function takes a survival fit and following the Cox model, estimates the hazard for individual from the model coefficient, while the baseline hazard is estimated with the Nelson-Aalen method.
#'
#' @param
#' obj : survival model object for example a glmnet fit. \cr
#' coef: coefficient beta is the default. \cr
#' new_data: a test holdout dataframe default test\cr
#' train_data : train dataframe default train\cr
#' pred : prediction error Brier, ROC or C-Index
#' adapted : in Lasso allows to swap to "adapated Lasso", otherwise will take usual lasso \cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
pred_error <- function(obj, train_data = train, test_data = test, time = os_months, status = os_deceased,  pred = "Brier", fit = "stepwise", adapted = "default"){
  time <- dplyr::enquo(time)
  status <- dplyr::enquo(status)

  # create coxph object with pre-defined coefficients
  if(fit == "Univariate" | fit == "stepwise" | fit == "lasso" | fit == "random forest" | fit == "PLS" | fit == "PCR"){

    if(fit == "Univariate"){
      selectedBeta <- obj
    }

    if(fit == "stepwise"){
      # take optimal beta from model object
      optimal.beta <- coef(obj)
      # find non zero beta coef
      nonzero.coef <- abs(optimal.beta)>0 & !is.na(optimal.beta)
      selectedBeta <- optimal.beta[nonzero.coef]
    }

    if(fit == "lasso"){
      # find lambda for which dev.ratio is max
      optimal.coef <- as.matrix(coef(obj, s = "lambda.min") )
      selectedBeta <-  optimal.coef[optimal.coef != 0,]
    }
    if(fit == "random forest" | fit == "PLS" | fit == "PCR"){
      mod <- obj;
      if(fit == "random forest"){
        selectedBeta <- train %>% dplyr::select(-!!time, -!!status)
      }
      if(fit == "PCR"){
        selectedBeta <- mod$feature.scores
      }
    }
    if(length(selectedBeta) == 0){
      selectedBeta <- 1;
      names(selectedBeta) <- "Null"
    }
    # take only covariates for which beta is not zero
    trainX <- train_data %>% dplyr::select(-!!time, -!!status)
    selectedTrainX   <- trainX[,colnames(trainX) %in% names(selectedBeta)]

    traincoxphdata <- cbind(train_data %>% dplyr::select(!!time, !!status), selectedTrainX)

    #create formula
    testBeta <- paste(names(selectedBeta), collapse = " + ")
    form <- as.formula(paste("Surv(time, status)", " ~ ", testBeta))

    if(fit == "Univariate" | fit == "stepwise" | adapted == "adapted lasso"){
      #Fit model
      mod <-  rms::cph(form , data = traincoxphdata %>% dplyr::mutate(
        time = !!time,
        status = !!status), surv = TRUE, x=TRUE, y=TRUE)
    }

    if(fit == "lasso" | fit == "PCR"){
      #Fit model
      mod <-  rms::cph(form , data = traincoxphdata %>% dplyr::mutate(
        time = !!time,
        status = !!status), init = selectedBeta, iter = 0, surv = TRUE)
    }

    # Create test vars
    testX <- test_data %>% dplyr::select(-!!time, -!!status)
    selectedTestX <- testX[,colnames(testX) %in% names(selectedBeta)]
    ndata <- cbind(test_data %>% dplyr::select(!!time, !!status), selectedTestX)

    #grid of equidistant time points
    timepoints <-  seq(0, max(train_data %>% dplyr::select(!!time) %>% unlist), length.out = 100L)

    if(pred == "Brier"){
      #Calculate probs
      probs <- pec::predictSurvProb(mod,
                      newdata = ndata, times = timepoints)
      #Calculate brier score
      out <- pec::pec(probs, Surv(time, status) ~ 1,
            data = ndata %>% dplyr::mutate( time   = !!time, status = !!status),
            maxtime = max(timepoints),
            exact = F,
            exactness = 99L)
    }
    if(pred == "ROC"){
      probs <- predict(mod, newdata = testX, type = "lp")
      out <- tdROC::tdROC(X = probs,
        Y = test_data %>% dplyr::select(time = !!time)%>% unlist,
         delta = test_data %>% dplyr::select(status = !!status)%>% unlist,
  tau = quantile(test_data %>% dplyr::select(time = !!time)%>% unlist, .67),
  span = 0.1, nboot = 10, alpha = 0.05, n.grid = 1000, cut.off = 5:9
      )
      attr(out, 'fit.model.roc') <- fit
    }
    if(pred == "c_index"){
      surv.obj=with(ndata %>%
          dplyr::mutate( time   = !!time, status = !!status),
          Surv(time,status))
      ###Create your survival estimates
      estimates= rms::survest(mod,newdata=testX,
                              times = max(timepoints))$surv
      ###Determine concordance
      out <- Hmisc::rcorr.cens(x=estimates,S=surv.obj)
    }

  }else{
    print("Not allowed fit argument")
  }

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

  ibrier <- pec::crps(brier, models = "matrix")

  return(ibrier)

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
                        feature350 = c(1,0,-1),
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

