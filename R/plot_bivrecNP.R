########################    plot.bivrecNP     ########################
#' Plot Results of Non-Parametric Analysis of Bivariate Recurrent Events
#'
#' @description
#' This function plots the joint CDF, marginal survival and conditioncal cdf obtained from \code{bivrecNP}().
#' For examples see \code{bivrecNP}.
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#'
#' @param x an object of class \code{bivrecNP}.
#' @param y either empty or NULL
#' @param main Optional string with plot title. Default is no title.
#' @param xlab Optional string with label for horizontal axis. Default is "Gap Times".
#' @param ylab Optional string with label for vertical axis. Default is "Individual".
#' @param type Optional vector of strings to label type 1and type 2 gap times. Default is c("Type 1", "Type 2").
#' @param ... arguments to be passed to graphical methods as needed.
#'
#' @export
#'

plot.bivrecNP <-function(x, y=NULL, type = NULL,
                         main = NULL, xlab = NULL, ylab = NULL, ...){

  if (!inherits(x, "bivrecNP")) stop("Object must be a bivrecNP class")

  #check arguments for labels
  if (missing(type)) {type=c("Type 1","Type 2")}
  if (missing(xlab)) {xlab="Gap Times"}
  if (missing(ylab)) {ylab="Individual"}
  if (missing(main)) {main=""}

  args = c(main, xlab, ylab, type)
  cond=x$conditional #boolean saying if conditional is in bivrecNP object

  if (cond==FALSE){
    plotJoint(x)
    plotMarg(x)
  }
  else {
    plotJoint(x)
    par(mar=c(1,1,1,1))
    par(mfrow=c(1,2))
    plotMarg(x)
    plotCond(x)
    par(mfrow=c(1, 1))
  }
}
