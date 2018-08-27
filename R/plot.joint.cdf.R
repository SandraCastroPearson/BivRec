
#' Plot Joint CDF after non-parametric analysis
#'
#' @description
#' This function plots joint cdf for recurrent event data.  Called from biv.rec.np(). No user interface.
#'
#' @param bivrec.nonparam.result List with joing.cdf, formula, data. Passed from biv.rec.np()
#' @param CI Confidence level for CI. Passed from biv.rec.np()
#'
#' @importFrom rgl open3d
#' @importFrom rgl plot3d
#' @importFrom rgl rgl.lines
#' @return A 3D scatter plot of joint cdf with confidence interval.
#'
#' @keywords internal
#'
plot.joint.cdf <- function(bivrec.nonparam.result, CI) {

  forplot <- bivrec.nonparam.result$cdf
  formula <- bivrec.nonparam.result$formula
  data <- bivrec.nonparam.result$data

  #####Wald CI and plot
  conf.lev = 1 - ((1-CI)/2)
  forplot$lower <- forplot[,3] - qnorm(conf.lev)*forplot[,4]
  forplot$upper <- forplot[,3] + qnorm(conf.lev)*forplot[,4]
  index <- which(forplot$lower<0)
  forplot[index, -1] <- forplot[index[1]-1, -1]

  open3d()
  plot3d(forplot[,1], forplot[,2], forplot[,3], col = "black", xlab = "x",
         main = "Joint cdf", ylab ="y", zlab = TeX('$P(X^0 \\leq x, Y^0\\leq y)$'),  expand = 1.09)
  for (i in 1:nrow(forplot)) {
    rgl.lines(forplot[i,1], forplot[i,2], as.matrix(forplot[i,5:6]), col="red")
  }


}
