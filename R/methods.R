#' Contrast two cost functions
#' @param cost1 cost function 1
#' @param cost2 cost function 2
#' @param nzv near zero variance.
#' @export
plot_contrast_cost <- function(cost1, cost2, nzv, xlab, ylab, title, subtitle = "") {
  if(missing(nzv)){
    nzv <- rep(F, length(cost1))
  }
  plot_data <- data.frame(x = cost1,
                          y = cost2,
                          nzv = nzv) %>%
    dplyr::mutate(nzv = ifelse(nzv,
                               paste("Non variable feature : " ,
                                     sum(nzv), sep ="" ),
                               paste("Variable feature : ",
                                     sum(!nzv), sep = "" ) ) )

  p <-  ggplot(plot_data,
               aes(x, y, color = nzv)) +
    geom_point(size = 0.5, alpha = .9) +
    ggpubr::theme_pubr() +
    labs(x = xlab, y = ylab, color = "") +
    ggtitle(title, subtitle)
  p
}
