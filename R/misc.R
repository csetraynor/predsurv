'%!in%' <- function(x,y)!('%in%'(x,y))


.check_args <- function(dat, x_vars, time_var, status_var) {
  if(!missing(x_vars)){
    if(any(x_vars %!in% colnames(dat))) {
      stop("Not all features in column names passed dataframe.")
    }
  }
  if(time_var %!in% colnames(dat)) {
    stop("time_var not in column names passed dataframe.")
  }
  if(status_var %!in% colnames(dat)) {
    stop("status_var not in column names passed dataframe.")
  }
}

.check_x_vars <- function(dat, x_vars, time_var, status_var) {
  if(missing(x_vars)) {
    warning("x_vars not passed using all vars as default.")
    time_var <- match(time_var, colnames(dat))
    status_var <- match(status_var, colnames(dat))
    x_vars <- colnames(dat)[-c(time_var, status_var)]
  }
  return(x_vars)
}


.check_data <- function(obj) {
  warning("Check data test needs to check if all tested levels of factor variables are in the test dataset.")
  as.data.frame(obj)
}
