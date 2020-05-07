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
padd_coefs <- function(x_vars, coefs) {
  dummy_coefs <- rep(0.0, length(x_vars))
  coefs <- dummy_coefs
  coefs[match( names(coefs), x_vars)] <- coefs
  names(coefs) <- x_vars
  coefs
}
