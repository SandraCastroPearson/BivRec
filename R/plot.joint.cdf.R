
#' Plot Joint CDF after non-parametric analysis
#'
#' @description
#' This function plots joint cdf for recurrent event data.  Called from biv.rec.np(). No user interface.
#'
#' @param bivrec.nonparam.result List with joing.cdf, formula, data. Passed from biv.rec.np()
#' @param CI Confidence level for CI. Passed from biv.rec.np()
#'
#' @return A 3D scatter plot of joint cdf with confidence interval.
#'
#' @importFrom stats ftable
#' @importClassesFrom grDevices heat.colors
#' @keywords internal
#'
plot.joint.cdf <- function(bivrec.nonparam.result, CI) {

  forplot <- bivrec.nonparam.result$cdf

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
                           forplot2, color.palette = heat.colors, cex.main=1.5,
                           xlab="x", ylab="y", main = expression(P(X^0 <= x, Y^0 <= y)))

}


