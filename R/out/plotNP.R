
#' Plot graphs (joint cdf, marginal and conditional cdf (if user wants)) after non-parametric analysis
#'
#' @description
#' This function plots joint cdf, marginal survival and conditional cdf for recurrent event data.  Called from biv.rec.np(). No user interface.
#'
#' @param object an object of class \code{bivrecNP}
#' @param main for labels
#' @param xlab for labels
#' @param ylab for labels
#' @param type1 for labels
#' @param type2 for labels
#' @param data only use for formula object
#'
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom graphics legend
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @export
#'
plot.bivrecNP <-function(object, main=NULL, xlab=NULL, ylab=NULL, type1=NULL, type2=NULL, data=NULL) {
  if (!is.bivrecNP(x)) stop("Object must be a bivrecNP class")
  cond=x$conditional #boolean saying if conditional is in bivrecNP object
  if (cond==FALSE){
    plotJoint(x)
    plotMarg(x)
  }
  else {
    plotJoint(x)
    #par(mfrow=c(1,2))
    plotMarg(x)
    plotCond(x)
  }
}


