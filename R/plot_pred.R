#' Plot Brier Score
#'
#' Plot Brier Score
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
plot_brier <- function(brier_score){
  Brier_Score = brier_score$AppErr[[2]]
  Time = brier_score$time
  plot(Brier_Score ~ Time, type= "l")
}

#' Plot ROC Curves
#'
#' Take various ROC curves and plot, a similar function can be found in the caret Tutorial.
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
#' @author Tarj

roc.plot2 <- function(...,x.var="false positive rate",y.var="true positive rate",x.lab=x.var,y.lab=y.var,ident=TRUE){
  fitList <- list(...)

  roc.data <- plyr::ldply(seq_along(fitList),function(i){
    roc.data <- data.frame(
      x = 1 -unlist(fitList[[i]]$ROC["spec"]),
      y = unlist(fitList[[i]]$ROC["sens"]),
      model=paste0(paste(deparse(attr(fitList[[i]], "predction.of.model") ),collapse = "\n"), " (AUC = ", round(fitList[[i]]$AUC[1], 2),
              ", lower = ",round(fitList[[i]]$AUC[3], 2),
              ", upper = ", round(fitList[[i]]$AUC[4], 2), ")" ) )

    arrange(roc.data,y)
  })

  p <- ggplot(roc.data, aes(x, y)) + theme_bw() +
    geom_line(aes(colour=model)) +
    guides(col = guide_legend(ncol = 1,title=NULL)) + theme(legend.position="bottom") +
    labs(title = "time dependent ROC")
  if(ident) p <- p + geom_abline(slope=1,alpha=0.2)
  p <- p + scale_x_continuous(x.lab)
  p <- p + scale_y_continuous(y.lab)
  p
}

#' Plot Brier Score
#'
#' Plot Brier Score
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
plot_brier2 <- function(...,x.var="Time",y.var="Brier Score",x.lab=x.var,y.lab=y.var,ident= FALSE){
  brierList <- list(...)
  ibrier <- lapply(brierList, function(i) pec::crps(i, models = "matrix"))

  brier.data <- plyr::ldply(seq_along(brierList),function(i){
    roc.data <- data.frame(
      x = unlist(brierList[[i]]$time) ,
      y = unlist(brierList[[i]]$AppErr[[2]]),
      model=paste0(paste(deparse(attr(brierList[[i]], "predction.of.model") ),collapse = "\n"), " (integrated Brier Score = ", round(ibrier[[i]][1], 2) ,")" ) )

  })

  p <- ggplot(brier.data, aes(x, y)) + theme_bw() +
    geom_line(aes(colour=model)) +
    guides(col = guide_legend(ncol = 1,title=NULL)) +
    theme(legend.position="bottom")+
    labs(title = "Brier score", subtitles = paste0("Time = ", deparse(round(brierList[[1]]$maxtime)) ))
  if(ident) p <- p + geom_abline(slope=1,alpha=0.2)
  p <- p + scale_x_continuous(x.lab)
  p <- p + scale_y_continuous(y.lab)
  p

}

