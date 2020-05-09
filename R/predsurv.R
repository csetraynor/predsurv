#' Fit univariate survival
#' Get the cost and coefficients for univariate survival analysis.
#' @param dat dataframe
#' @param x_vars vars to analysis
#' @param ncores for parallel.
#' @param time_var defults to time
#' @param status_var defaults to status
#' @param cost cost measure
#' @param save_file should results be saved on file
#' @importFrom parallel mclapply
#' @export
univariate_survival <- function(d_train,
                                x_vars,
                                d_test = NULL,
                                weights = NULL,
                                ncores = 1L,
                                time_var = "time",
                                status_var = "status",
                                cost = "ll",
                                mod = c("cox", "reg"),
                                save_file) {

  .check_args(d_train, x_vars, time_var, status_var)
  x_vars <- .check_x_vars(d_train, x_vars, time_var, status_var)
  base_form <- get_surv_form_lefthand(time_var, status_var)
  allowed_cost <- c("ll","concordance", "gamma", "harrelc","coef")
  if(cost %!in% allowed_cost) {
    stop("cost function ", cost, " not recognised. ")
  }

  mod_fit_specific <- paste0("mod_fit_", mod)
  mod_fit_specific <- get(mod_fit_specific, mode = "function")

  get_cost_specific <- paste0("get_cost_", cost)
  get_cost_specific <- get(get_cost_specific, mode = "function")

  if (is.null(d_test)) {
    d_test <- d_train
    d_type <- 'train'
  } else {
    d_test <- .check_data(d_test)
    d_type <- 'test'
  }

  if(!is.null(weights)) {
    w <- weights
  } else {
    w <- rep(1, nrow(d_test))
    w <- w/sum(w)
  }

  if(ncores > 1) {
    res <- unlist( parallel::mclapply(x_vars, function(var){
      survform <- as.formula(paste0(base_form, var ))
      mod <- mod_fit_specific( d_train, survform, time_var, status_var)
      return( get_cost_specific(x = mod, dat = d_test, time = time_var, status = status_var, weights = w )) }, mc.cores = ncores) )
  } else {
    res <- sapply(x_vars, function(var){
      survform <- as.formula(paste0(base_form, var ))
      mod <- mod_fit_specific( d_train, survform, time_var, status_var)
      return( get_cost_specific(x = mod, dat = d_test, time = time_var, status = status_var, weights = w ))    })
  }
  names(res) <- x_vars
  if(missing(save_file)) {
    return(res)
  } else {
    saveRDS(res, save_file)
    return(0)
  }
}

#' Fit and select univariate selection survival model
#' Get the cost and coefficients for univariate survival analysis.
#' @param d_train dataframe
#' @param g additional data frame for bootstrapping
#' @param x_vars vars to analysis
#' @param ncores for parallel.
#' @param time_var defults to time
#' @param status_var defaults to status
#' @param save_file should results be saved on file
#' @param threshold p value
#' @param correction correction to apply.
#' @importFrom parallel mclapply
#' @export
univariate_sel <- function(d_train, g, x_vars, time_var = "time", status_var = "status", threshold = 0.01, correction = "bonferroni", ncores =1L, save_file){
  if(!missing(g)){
    d_train <- flirt_rows(d_train, g)
  }
  .check_args(d_train, x_vars, time_var, status_var)
  x_vars <- .check_x_vars(d_train, x_vars, time_var, status_var)
  base_form <- get_surv_form_lefthand(time_var, status_var)
  if(ncores > 1) {
    wald_test <- unlist( parallel::mclapply(x_vars, function(var){
      survform <- as.formula(paste0(base_form, var ))
      mod <- mod_fit_cox( d_train, survform, time_var, status_var)
      return( get_waldtest(x = mod)) }, mc.cores = ncores) )
  } else {
    wald_test <- sapply(x_vars, function(var){
      survform <- as.formula(paste0(base_form, var ))
      mod <- mod_fit_cox( d_train, survform, time_var, status_var)
      return( get_waldtest(x = mod ))    })
  }
  wald_test <- p.adjust(wald_test, correction)
  sel_vars <- !is.na(wald_test) & wald_test < threshold
  ## fit multivariate model
  x_sel <- x_vars[sel_vars]
  surv_form <-  get_surv_form(time_var, status_var, x_sel)
  fit <- coxph(surv_form, data = d_train)
  res <- padd_coefs(x_vars, get_cost_coef(fit))

  if(missing(save_file)) {
    return(res)
  } else {
    saveRDS(res, save_file)
    return(0)
  }
}

#' Fit forward stepwise selection survival model
#' Get the cost and coefficients for forward search survival analysis.
#' @param d_train dataframe
#' @param g additional data frame for bootstrapping
#' @param x_vars vars to analysis
#' @param ncores for parallel.
#' @param time_var defults to time
#' @param status_var defaults to status
#' @param save_file should results be saved on file
#' @param penalty information criteria to use: AIC or BIC
#' @importFrom parallel mclapply
#' @export
forward_sel <- function(d_train, g, x_vars, time_var = "time", status_var = "status", penalty = "AIC", save_file){
  if(!missing(g)){
    d_train <- flirt_rows(d_train, g)
  }
  .check_args(d_train, x_vars, time_var, status_var)
  x_vars <- .check_x_vars(d_train, x_vars, time_var, status_var)
  if("time" %!in% colnames(d_train)){
    d_train$time   <- d_train[[time_var]]
  }
  if("status" %!in% colnames(d_train)){
    d_train$status <- d_train[[status_var]]
  }
  fit0 <- coxph(Surv(time, status) ~ 1, data = d_train,
                x = TRUE, y = TRUE)
  k <- ifelse(penalty == "AIC", 2,  log(nrow(d_train)))
  testBeta <- reformulate(x_vars)
  fit1 <- step(fit0, direction = "forward", scope = list(upper = testBeta, lower = ~1), k = k, trace = F)
  res <- padd_coefs(x_vars, get_cost_coef(fit))

    if(missing(save_file)) {
    return(res)
  } else {
    saveRDS(res, save_file)
    return(0)
  }
}

#' Fit elastic net selection survival model
#' Get the cost and coefficients for elastic net survival analysis.
#' @param d_train dataframe
#' @param g additional data frame for bootstrapping
#' @param x_vars vars to analysis
#' @param ncores for parallel.
#' @param time_var defults to time
#' @param status_var defaults to status
#' @param alpha alpha for elastic net
#' @param save_file should results be saved on file
#' @param relaxed apply relaxed lasso set to TRUE
#' @param correction correction to apply.
#' @importFrom doMC registerDoMC
#' @importFrom glmnet cv.glmnet
#' @export
lasso_sel <- function(d_train, g, x_vars, time_var = "time", status_var = "status", alpha = 1, relaxed = F, ncores = 1L, save_file, ...) {
  if(!missing(g)){
    d_train <- flirt_rows(d_train, g)
  }
  .check_args(d_train, x_vars, time_var, status_var)
  x_vars <- .check_x_vars(d_train, x_vars, time_var, status_var)
  if(relaxed) {
    alpha <- 1
  }
  # fit the classical Lass model with glmenet
  X <- as.matrix(d_train[, x_vars, drop = F])
  Y <- as.matrix(d_train[, c(time_var, status_var), drop = F])
  colnames(Y) <- c("time", "status")

  if(ncores > 1 ){
    if(!require(doMC)) {
      stop("Package doMC required for parallel computation in cv.glmnet")
    }
    doMC::registerDoMC(cores=ncores)
    enetfit <-  glmnet::cv.glmnet(X, Y, family = "cox", grouped = TRUE, alpha = alpha, parallel = TRUE, ...)
  } else {
    enetfit <-  glmnet::cv.glmnet(X, Y, family = "cox", grouped = TRUE, alpha = alpha, parallel = FALSE, ...)
  }
  optimal_coef <- t(as.matrix(coef(enetfit, s = "lambda.min")))[1, ]

  if(!relaxed){
    res <- optimal_coef
  } else {
    sel_vars <- names(optimal_coef[abs(optimal_coef) > 0.0])
    enet_form <- get_surv_form(time_var, status_var, sel_vars)
    modenet <- survival::coxph(enet_form, data = d_train,
        init = optimal_coef[abs(optimal_coef) > 0.0],
        x = T, y = T, iter.max = 5)
    res <- padd_coefs(x_vars, get_cost_coef(modenet))
  }
  if(missing(save_file)) {
    return(res)
  } else {
    saveRDS(res, save_file)
    return(0)
  }
}

#' Fit projective prediction survival model
#' Get the cost and coefficients for projective prediction survival analysis.
#' @param d_train dataframe
#' @param g additional data frame for bootstrapping
#' @param x_vars vars to analysis
#' @param ncores for parallel.
#' @param time_var defults to time
#' @param status_var defaults to status
#' @param alpha alpha for elastic net
#' @param save_file should results be saved on file
#' @param relaxed apply relaxed lasso set to TRUE
#' @param p0  prior guess for the number of relevant variables
#' @importFrom rstanarm stan_surv
#' @importFrom projpred cv_varsel suggest_size project
#' @export
project_sel <- function(d_train, g, x_vars, time_var = "time", status_var = "status", method = "L1", p0 = 5, nv_max = 40, save_file, precond_cutoff, cores = 1L, ...) {
  if(!missing(g)){
    d_train <- flirt_rows(d_train, g)
  }
  .check_args(d_train, x_vars, time_var, status_var)
  x_vars <- .check_x_vars(d_train, x_vars, time_var, status_var)
  x_vars_back <- x_vars
  if(!missing(precond_cutoff)) {
    precond_vars <- univariate_survival(d_train,
                        x_vars,
                        d_test = NULL,
                        weights =rep(1, nrow(d_train)),
                        ncores = cores,
                        time_var = time_var,
                        status_var = status_var,
                        cost = "ll",
                        mod = "reg")
    x_vars <- .ll_precond(x_vars, precond_vars, precond_cutoff)
  }
  D <- length(x_vars)
  n <- nrow(d_train)
  sigma <- 2 # approximate plug-in value for observation information
  tau0 <- p0/(D-p0) * sigma/sqrt(n)
  prior_coeff <- rstanarm::hs(global_scale = tau0, slab_scale = 1)
  # fit the full model
  surv_formula <- get_surv_form(time_var, status_var , x_vars)
  mod_surv <- rstanarm::stan_surv(formula = surv_formula,
                                  data    = d_train,
                                  prior   = prior_coeff,
                                  cores   = cores,
                                  ...)
  cvs <-  projpred::cv_varsel(mod_surv)
  size <- projpred::suggest_size(cvs)
  vind <- cvs$vind[seq_len(size)]
  res <- projpred::project(mod_surv, vind = vind)
  res <- padd_coefs(x_vars_back, apply(as.matrix(res), 2, mean)[-1])
  if(missing(save_file)) {
    return(res)
  } else {
    saveRDS(res, save_file)
    return(0)
  }
}



#' Fit stan survival model
#' Fit stan survival model
#' @param d_train dataframe
#' @param g additional data frame for bootstrapping
#' @param x_vars vars to analysis
#' @param ncores for parallel.
#' @param time_var defults to time
#' @param status_var defaults to status
#' @param alpha alpha for elastic net
#' @param save_file should results be saved on file
#' @param relaxed apply relaxed lasso set to TRUE
#' @param p0  prior guess for the number of relevant variables
#' @importFrom rstanarm stan_surv
#' @importFrom projpred cv_varsel suggest_size project
#' @export
fit_stan_surv <- function(d_train, g, x_vars, time_var = "time", status_var = "status", method = "L1", p0 = 5, nv_max = 40, save_file, ...) {
  if(!missing(g)){
    d_train <- flirt_rows(d_train, g)
  }
  .check_args(d_train, x_vars, time_var, status_var)
  x_vars <- .check_x_vars(d_train, x_vars, time_var, status_var)
  D <- length(x_vars)
  n <- nrow(d_train)
  tau0 <- p0/(D-p0) * 1/sqrt(n) # scale for tau
  prior_coeff <- rstanarm::hs(global_scale = tau0, slab_scale = 1)
  # fit the full model
  surv_formula <- get_surv_form(time_var, status_var ,x_vars)
  mod_surv <- rstanarm::stan_surv(formula = surv_formula,
                                  data    = d_train,
                                  prior   = prior_coeff,
                                  ...)
  cvs <- projpred::cv_varsel(mod_surv, method = method, nv_max = nv_max)
  if(missing(save_file)) {
    return(list(survfit = mod_surv, cvs = cvs))
  } else {
    saveRDS(list(survfit = mod_surv, cvs = cvs), save_file)
    return(0)
  }
}



######## tree
tree_sel <- function(d_train, x_vars, time_var = "time", status_var = "status") {
  .check_args(d_train, x_vars, time_var, status_var)
  x_vars <- .check_x_vars(d_train, x_vars, time_var, status_var)
  surv_form <- get_surv_form(time_var, status_var, x_vars)
  mod <-  party::cforest(surv_form, data = d_train, controls = cforest_unbiased(ntree = 20))
  # t1 <- Sys.time()
  # vi_res <- varimp(mod)
  # t2 <- Sys.time()
  # t2 - t1
}
