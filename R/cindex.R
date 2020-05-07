cindex_concordance <- function(x, ...) {
  time <- x[,1]
  status <- x[,2]
  mu <- x[,3]
  if(ncol(x) == 4) {
    w <- x[,4]
  } else {
    w <- rep(1, nrow(x))
    w <- w/sum(w)
  }
  y <- survival::Surv(time, status)
  cindex <- survival:::survConcordance.fit(y = y, x = mu, weight = w)
  agree <- cindex[1]
  disagree <- cindex[2]
  tied <- cindex[3]
  (agree + tied/2) / (agree + disagree + tied)
  # (agree-disagree)/(agree+disagree)
}

cindex_gamma <- function(x, ...) {
  time <- x[,1]
  status <- x[,2]
  mu <- x[,3]
  if(ncol(x) == 4) {
    w <- x[,4]
  } else {
    w <- rep(1, nrow(x))
    w <- w/sum(w)
  }
  y <- survival::Surv(time, status)
  cindex <- survival:::survConcordance.fit(y = y, x = mu, weight = w)
  agree <- cindex[1]
  disagree <- cindex[2]
  tied <- cindex[3]
  #(agree + tied/2) / (agree + disagree + tied)
  (agree-disagree)/(agree+disagree)
}


cindex_harrell <- function( x, ...) {
  if(!require(survcomp)){
    stop("survcomp package needed for function cindex_harrell.")
  }
  time <- x[,1]
  status <- x[,2]
  mu <- x[,3]
  if(ncol(x) == 4) {
    w <- x[,4]
  } else {
    w <- rep(1, nrow(x))
    w <- w/sum(w)
  }

  survcomp::concordance.index(x = mu, surv.time = time, surv.event = status, weights = w, ...)$c.index
}


get_cost <- function(x){
  x$loglik[2]
}
