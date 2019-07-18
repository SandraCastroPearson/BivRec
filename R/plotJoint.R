
#' Plot Joint CDF after non-parametric analysis
#'
#' @description
#' This function plots joint cdf for recurrent event data.
#' @noRd
#' @param x must be an object of \code{bivrecNP} class.
#'
#' @return A contour plot of joint cdf.
#'
#' @importFrom stats ftable
#' @examples
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' bdat <- is.bivrecSurv(bivrec_data)
#' npresult <- bivrecNP(bdat,ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),conditional = FALSE, given.interval=c(0, 10))
#' plotJoint(npresult)
#'
plotJoint <- function(x) {
  if (!is.bivrecNP(x)) stop("Object must be a bivrecNP class")

  forplot <- x$joint_cdf

  #####OLD MAY RE-USE LATTER: Wald CI and plot
  # rgl::plot3d(forplot[,1], forplot[,2], forplot[,3], col = "black", xlab = "x",
  #        main = "Joint cdf", ylab ="y", zlab = expression(P(X^0 <= x, Y^0 <= y)),  expand = 1.1, cex.lab = 1.5)
  # for (i in 1:nrow(forplot)) {
  #   rgl::rgl.lines(forplot[i,1], forplot[i,2], as.matrix(forplot[i,5:6]), col="red")
  # }


  forplot <- forplot[1:3]
  colnames(forplot) <- c("X", "Y", "Cumm.Prob")
  myx <- as.factor(forplot$X)
  myy <- as.factor(forplot$Y)
  lx <- length(levels(myx))
  forplot2 <- matrix(ftable(forplot, row.vars = 1, col.vars = 2), nrow=lx)
  rownames(forplot2) <- levels(myx)
  colnames(forplot2) <- levels(myy)
  for (i in 1:lx) {
      index <- which(forplot$X==as.numeric(rownames(forplot2)[i]))
      forplot2[i,] = forplot$Cumm.Prob[index]
  }
  graphics::filled.contour(x=as.numeric(levels(myx)), y= as.numeric(levels(myy)),
                           forplot2, color.palette = grDevices::heat.colors, cex.main=1.5,
                           xlab="x", ylab="y", main = expression(P(X^0 <= x, Y^0 <= y)))

}


