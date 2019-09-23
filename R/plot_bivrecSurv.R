########################    plot.bivrecSurv     ########################
#' Plot Bivariate Alternating Recurrent Series
#'
#' @description
#' This function plots bivariate recurrent event gap times from a bivrecSurv object.
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#'
#' @param x An object of class \code{bivrecSurv}.
#' @param y either empty or NULL
#' @param main Optional string with plot title. Default is no title.
#' @param xlab Optional string with label for horizontal axis. Default is "Gap Times".
#' @param ylab Optional string with label for vertical axis. Default is "Individual".
#' @param type Optional vector of strings to label type 1and type 2 gap times. Default is c("Type 1", "Type 2").
#' @param ... arguments to be passed to graphical methods as needed.
#'
#' @export
#'
#' @examples
#'# Simulate bivariate alternating recurrent event data
#' library(BivRec)
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' plot(x = with(bivrec_data, bivrecSurv(id, epi, xij, yij, d1, d2)), main="Example")
#'

plot.bivrecSurv <- function(x, y=NULL, type = NULL,
                            main = NULL, xlab = NULL, ylab = NULL, ...){
  if (!inherits(x, "bivrecSurv")) stop("Object must be a bivrecSurv class")
  object <- x

  #check arguments for labels
  if (missing(type)) {type=c("Type 1","Type 2")}
  if (missing(xlab)) {xlab="Time"}
  if (missing(ylab)) {ylab="Individual"}
  if (missing(main)) {main=""}

  args = c(main, xlab, ylab, type)

  ##EXTRACT VECTORS FOR PLOTTING FUNCTION
  parameters <- object$data4Creg[-(5:7)]
  colnames(parameters) <- c("id", "episode", "xij", "yij", "ci")
  ctimes <- object$data4Lreg$ctime
  nsubject <- object$data4Lreg$n
  temp <- NULL

  basicplot(parameters, ctimes, nsubject, temp, args, a = 1/2, b = 1/3)

}
