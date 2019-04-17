
#' Plot graphs (joint cdf, marginal and conditional (if user wants)) after non-parametric analysis
#'
#' @description
#' This function plots joint cdf, marginal survival and conditional for recurrent event data.  Called from biv.rec.np(). No user interface.
#'
#' @param x an object of class \code{bivrecNP}

#'
#' @return A 2 (or 3) plots
#'
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom graphics legend
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @export
#'
plot.bivrecNP <-function(x) {
  if (!is.bivrecNP(x)) stop("Object must be a bivrecNP class")
  cond=x$conditional #boolean saying if conditional is in bivrecNP object 
  if (cond==FALSE){
    plotJoint(x)
    plotMarg(x)
  }
  else {
    plotJoint(x)
    par(mfrow=c(1,2))
    plotMarg(x)
    plotCond(x)
  }
}
  

