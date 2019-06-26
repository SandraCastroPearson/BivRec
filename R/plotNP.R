#' Plot Joint cdf, Marginal and Conditional cdf (if user desires) after non-parametric analysis
#'
#' @description
#' This function plots joint cdf and marginal survival by default and conditional cdf (if user wants) for recurrent event data.
#'
#' @param x an object of class \code{bivrecNP}
#'
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom graphics legend
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @export
#'
#' @examples
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' bdat <- is.bivrecSurv(bivrec_data)
#' npresult <- bivrecNP(bdat,ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),conditional = FALSE, given.interval=c(0, 10))
#' plot(npresult)
#'
#' \dontrun{
#' #This is an example with longer runtime (it runs the conditional graph)
# npresult2 <- bivrecNP(bdat,ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),conditional = TRUE, given.interval=c(0, 10))
# plot(npresult2)
#'}
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
    plotMarg(x)
    plotCond(x)
  }
}


