get_surv_form_lefthand <- function(time_var, status_var){
  paste0("Surv(", time_var, ", ", status_var, ") ~ ") }

get_surv_form <- function(time_var, status_var, vars){
  lh <- get_surv_form_lefthand(time_var, status_var)
  return(as.formula(paste0(lh, paste(vars, collapse = " + ")))) }

flirt_rows <- function(d, g) {
  g$id <- 1:nrow(g) ## dummy id
  left_join(d, g, by = "id") }

pasteplus <- function(x) {
  paste(x, collapse = " + ") }

pasteformula <- function(x) {
  as.formula( paste0("~", x) ) }

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

mod_fit_reg <- function(x, form, time_var = "time",
                    status_var = "status", ...)
  survival::survreg(form, data = x)

mod_fit_cox <- function(x, form, time_var = "time",
                        status_var = "status", ...)
  survival::coxph(form, data = x, ...)


get_cost_ll <- function(x, ...) {
  cost <- x$loglik[2]
  if(is.null(cost)) {
    warning("ll extraction was not succesfull check get_cost_ll for object of class ", class(x))
  }
  return(cost)
}

get_cost_coef <- function(x, ...) {
  coef <- x$coef
  if(is.null(coef)) {
    warning("coef extraction was not succesfull check get_cost_coef for object of class ", class(x))
  }
  return(coef)
}

get_cost_concordance <- function(x, dat, time, status, weights) {
  risk <- predict(x, newdata = dat)
  cost_data <- .get_cost_data(risk, dat, time, status, weights)
  cindex_concordance(cost_data)
}

get_cost_gamma <- function(x, dat, time, status, weights) {
  risk <- predict(x, newdata = dat)
  cost_data <- .get_cost_data(risk, dat, time, status, weights)
  cindex_gamma(cost_data)
}

get_cost_harrelc <- function(x, dat, time, status, weights) {
  risk <- predict(x, newdata = dat)
  cost_data <- .get_cost_data(risk, dat, time, status, weights)
  cindex_harrell(cost_data)
}


.get_cost_data <- function(risk, dat, time, status, w) {
  cost_data <- data.frame(
    time = dat[[time]],
    status = dat[[status]],
    risk = risk,
    weights = w
  )
  return(cost_data)
}

get_waldtest <- function(x) {
  summary(x)$waldtest[3]
}
padd_coefs <- function(coefs, x_vars) {
  dummy_coefs <- rep(0.0, length(x_vars))
  coefs <- dummy_coefs
  coefs[match( names(coefs), x_vars)] <- coefs
  names(coefs) <- x_vars
  coefs
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
  if(ncores > 1) {
    res <- unlist( parallel::mclapply(x_vars, function(var){
      survform <- as.formula(paste0(base_form, var ))
      mod <- mod_fit_cox( d_train, survform, time_var, status_var)
      return( get_waldtest(x = mod)) }, mc.cores = ncores) )
  } else {
    res <- sapply(x_vars, function(var){
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
  coefs <- padd_coefs(x_vars, get_cost_coef(fit))

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
  if("time" %!in% colnames(dat)){
    dat$time <- dat[[time_var]]
  }
  if("status" %!in% colnames(dat)){
    dat$status <- dat[[status_var]]
  }
  fit0 <- coxph(Surv(time, status) ~ 1, data = d_train,
                x = TRUE, y = TRUE)
  k <- ifelse(penalty == "AIC", 2,  log(nrow(d_train)))
  testBeta <- reformulate(x_vars)
  fit1 <- step(fit0, direction = "forward", scope = list(upper = testBeta, lower = ~1), k = k, trace = F)
  coefs <- padd_coefs(x_vars, fit$coefficients)

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
lasso_sel <- function(d_train, g, x_vars, time_var = "time", status_var = "status", alpha = 1, relaxed = F, ncores, save_file) {
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
    coefs <- optimal_coef
  } else {
    sel_vars <- names(optimal_coef[abs(optimal_coef) > 0.0])
    enet_form <- get_surv_form(time_var, status_var, sel_vars)
    modenet <- survival::coxph(enet_form, data = d_train,
        init = optimal_coef[abs(optimal_coef) > 0.0],
        x = T, y = T, iter.max = 5)
    coefs <- padd_coefs(x_vars, get_cost_coef(modenet))
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
project_sel <- function(dat, g, x_vars, time_var = "time", status_var = "status", method = "L1", p0 = 5, nv_max = 40, save_file, ...) {
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
  cvs <-  projpred::cv_varsel(mod_surv, method = method, nv_max = nv_max)
  size <- projpred::suggest_size(cvs)
  vind <- cvs$vind[1:size]
  res <- projpred::project(mod_surv, vind = vind)

  coefs <- padd_coefs(x_vars, apply(as.matrix(res), 2, mean)[-1])

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
fit_stan_surv <- function(dat, g, x_vars, time_var = "time", status_var = "status", method = "L1", p0 = 5, nv_max = 40, save_file, ...) {
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
tree_sel <- function(dat, x_vars, time_var = "time", status_var = "status") {
  surv_form <- get_surv_form(time_var, status_var, x_vars)
  mod <-  party::cforest(surv_form, data = dat, controls = cforest_unbiased(ntree = 20))
  # t1 <- Sys.time()
  # vi_res <- varimp(mod)
  # t2 <- Sys.time()
  # t2 - t1
}
