#######################    plotJoint     ########################
#' Plot Joint CDF After Nonparametric Analysis
#'
#' @description
#' This function plots the joint CDF for a \code{bivrecNP} object.
#'
#' @param object An object of \code{bivrecNP} class.
#' @importFrom stats ftable
#' @keywords internal
#'

#@return A contour plot of joint cdf.

plotJoint <- function(object) {

  x = object

  if (!inherits(x, "bivrecNP")) stop("Object must be a bivrecNP class")

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

########################    plotMarg     ########################
#' Plot Marginal Survival After Nonparametric Analysis
#'
#' @description
#' This function plots the marginal survival for a \code{bivrecNP} object.
#'
#' @param object An object of \code{bivrecNP} class.
#' @keywords internal

#@return A plot of marginal survival vs. first gap time with confidence interval.

plotMarg <- function(object) {
  x <- object

  if (!inherits(x, "bivrecNP")) stop("Object must be a bivrecNP class")
  xij <- x$xij
  forplot <- x$marginal_survival[1:3]
  #formula <- bivrec.nonparam.result$formula

  #variables <- all.vars(formula)
  mx <- round(max(xij), digits = 0)
  str_mx <- substring(as.character(mx), 1, nchar(as.character(mx))-1)
  str_mx <- paste(as.numeric(str_mx)+1, 0, sep="")
  mx <- round(as.numeric(str_mx), digits=1)
  forplot <- rbind(c(0, 1, 0), forplot, c(mx, 0, forplot[nrow(forplot),3]))

  ##### Wald CI and plot
  conf.lev = 1 - ((1-x$level)/2)
  forplot$lower <- forplot[,2] - qnorm(conf.lev)*forplot[,3]
  forplot$upper <- forplot[,2] + qnorm(conf.lev)*forplot[,3]
  index <- which(forplot$lower<0)
  forplot[index, -1] <- forplot[index[1]-1, -1]

  plot(forplot$Time, forplot$Marginal_Survival, type = "l", xlab = "Type I Gap Times (x)",
       ylab = "Marginal Survival", yaxp  = c(0, 1, 10),
       xaxp  = round(c(0, mx, 10), digits=1), main = expression(1 - P(X^0 <= x))
  )
  graphics::lines(forplot$Time, forplot$lower, lty = 2)
  graphics::lines(forplot$Time, forplot$upper, lty = 2)

}

########################    plotCond     ########################
#' Plot Conditional CDF After Nonparametric Analysis
#'
#' @description
#' This function plots conditional cdf for a \code{bivrecNP} object.
#'
#' @param object An object of \code{bivrecNP} class where the analysis has specified conditional = TRUE.
#' @importFrom stats ftable
#' @keywords internal

#@return A plot of conditional cdf in the given interval.

plotCond <- function(object) {
  x=object
  cond <-x$conditional_cdf$conditional
  plot(cond$Time, cond[,5], type="l", lty = 2, xlab = "Type II Gap Times (y)",
       ylab = "Conditional Probability", xlim=c(0, round(max(x$conditional_cdf$ygrid), digits=1)),
       ylim=c(0, round(max(cond[,5]), digits=1)),
       main=substitute(paste("P(", Y^0 <= y, "|", X^0 %in% "[", gi1, ",", gi2, "])"),
                       list(gi1 = x$given.interval[1], gi2 = x$given.interval[2]))
  )
  graphics::lines(cond$Time, cond[,4], lty = 2)
  graphics::lines(cond$Time, cond$Conditional.Probability,lty = 1)
}

########################    plot.bivrecNP     ########################
#' Plot Results of Non-Parametric Analysis of Bivariate Recurrent Events
#'
#' @description
#' This function plots all the estimated functions (joint CDF, marginal survival and/or conditioncal cdf) from a \code{bivrecNP}() object in one step.
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#'
#' @param x An object of class \code{bivrecNP}.
#' @param y Either empty or NULL
#' @param main Optional string with plot title. Default is no title.
#' @param xlab Optional string with label for horizontal axis.
#' @param ylab Optional string with label for vertical axis.
#' @param type Optional vector of strings to label type 1 and type 2 gap times. Default is c("Type 1", "Type 2").
#' @param ... Additional arguments to be passed to graphical methods if needed.
#'
#' @export
#'

plot.bivrecNP <-function(x, y=NULL, type = NULL,
                         main = NULL, xlab = NULL, ylab = NULL, ...){

  if (!inherits(x, "bivrecNP")) stop("Object must be a bivrecNP class")

  cond=x$conditional #boolean saying if conditional is in bivrecNP object

  if (cond==FALSE){
    par(mar=c(5,4,4,2)+0.1)
    plotJoint(x)
    par(mar=c(5,4,4,2)+0.1)
    plotMarg(x)
  }
  else {
    plotJoint(x)
    par(mar=c(5,4,4,2)+0.1, mfrow=c(1,2))
    plotMarg(x)
    plotCond(x)
    par(mfrow=c(1, 1))
  }
}
