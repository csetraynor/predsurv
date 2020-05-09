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

.ll_precond <- function(x_vars, precond_vars, precond_cutoff) {
  x_vars[x_vars %in% names(head( sort(-2*precond_vars) , precond_cutoff))]
}


#' @export
model_matrix <- function(X, xvars, ncores = 1){
  if(missing(xvars)){
    xvars <- colnames(X)
  }
  div <- length(xvars) %% 1000
  res <- round( length(xvars) / 1000, 0)
  if(ncores > 1) {
    resX <- parallel::mclapply(seq_len(res), function(j){
      ind <- seq( 1+(j-1)*1000, j*1000, by = 1)
      Xchunk <- X[,ind,drop=F]
      Xchunk <- as.data.frame(Xchunk) %>% dplyr::mutate_all(as.factor)
      return( model.matrix( ~ . , Xchunk )) },
      mc.cores = ncores)
    resX <- do.call(cbind, resX)
  } else {
    resX <- sapply(seq_len(res), function(j){
      ind <- seq( 1+(j-1)*1000, j*1000, by = 1)
      Xchunk <- X[,ind,drop=F]
      Xchunk <- as.data.frame(Xchunk) %>% dplyr::mutate_all(as.factor)
      model.matrix( ~ . , data = Xchunk )
    })
  }
  ## last chunk
  ind <- seq( 1+res*1000, res*1000+div, by = 1)
  Xchunk <- X[,ind,drop=F]
  Xchunk <- as.data.frame(Xchunk) %>% dplyr::mutate_all(as.factor)
  resX2 <- model.matrix( ~ . , Xchunk )
  cbind(resX, resX2)
}
