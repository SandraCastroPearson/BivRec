
#' Plot Joint CDF after non-parametric analysis
#'
#' @description
#' This function plots joint cdf for recurrent event data.  Called from biv.rec.np(). No user interface.
#'
#' @param bivrec.nonparam.result List with joing.cdf, formula, data. Passed from biv.rec.np()
#' @param CI Confidence level for CI. Passed from biv.rec.np()
#'
#' @importFrom rgl plot3d
#' @importFrom rgl rgl.lines
#' @return A 3D scatter plot of joint cdf with confidence interval.
#'
#' @keywords internal
#'
plot.joint.cdf <- function(bivrec.nonparam.result, CI) {

  forplot <- bivrec.nonparam.result$cdf
  #####Wald CI and plot

  plot3d(forplot[,1], forplot[,2], forplot[,3], col = "black", xlab = "x",
         main = "Joint cdf", ylab ="y", zlab = expression(P(X^0 <= x, Y^0 <= y)),  expand = 1.1, cex.lab = 1.5)
  for (i in 1:nrow(forplot)) {
    rgl.lines(forplot[i,1], forplot[i,2], as.matrix(forplot[i,5:6]), col="red")
  }

}
