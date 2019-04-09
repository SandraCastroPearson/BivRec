
#' Plot graphs (joint cdf, marginal and conditional (if user wants)) after non-parametric analysis
#'
#' @description
#' This function plots joint cdf, marginal survival and conditional for recurrent event data.  Called from biv.rec.np(). No user interface.
#'
#' @param bivrec.nonparam.result List with joint.cdf, formula, data. Passed from biv.rec.np()
#' @param CI Confidence level for CI. Passed from biv.rec.np()
#'
#' @return A 3D scatter plot of joint cdf with confidence interval.
#'
#' @importFrom stats ftable
#' @keywords internal
#'
plot.bivrecNP <- function(x) {
  if (!is.bivrecNP(x)) stop("Object must be a bivrecNP class")
  forplot <- x$joint.cdf
  
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

#marginal survival
  mgplot <- x$marginal.survival[1:3]
  #formula <- bivrec.nonparam.result$formula
  data <- x$df
  
  #variables <- all.vars(formula)
  xij <- x$df$xij
  mx <- round(max(xij), digits = 0)
  str_mx <- substring(as.character(mx), 1, nchar(as.character(mx))-1)
  str_mx <- paste(as.numeric(str_mx)+1, 0, sep="")
  mx <- round(as.numeric(str_mx), digits=1)
  mgplot <- rbind(c(0, 1, 0), mgplot, c(mx, 0, mgplot[nrow(mgplot),3]))
  
  #####95% Wald CI and plot
  conf.lev = 1 - ((1-x$CI)/2)
  mgplot$lower <- mgplot[,2] - qnorm(conf.lev)*mgplot[,3]
  mgplot$upper <- mgplot[,2] + qnorm(conf.lev)*mgplot[,3]
  index <- which(mgplot$lower<0)
  mgplot[index, -1] <- mgplot[index[1]-1, -1]
  
  #par(mfrow=c(1,2))
  #mar = c(0,0,0,0)
  #plot for marginal survival
  plot(mgplot$Time, mgplot$Marginal.Survival, type = "l", xlab = "Type I Gap Times (x)", ylab = "Marginal Survival",
       yaxp  = c(0, 1, 10), xaxp  = round(c(0, mx, 10), digits=1), main = expression(1 - P(X^0 <= x)))
  graphics::lines(mgplot$Time, mgplot$lower, lty = 2)
  graphics::lines(mgplot$Time, mgplot$upper, lty = 2)
  #plot for joint cdf 
  graphics::filled.contour(x=as.numeric(levels(myx)), y= as.numeric(levels(myy)),
                           forplot2, color.palette = grDevices::heat.colors, cex.main=1.5,
                           xlab="x", ylab="y", main = expression(P(X^0 <= x, Y^0 <= y)))
  
  
}

#have to figure out how to get heatmap and plot(ggplot?) side by side.
