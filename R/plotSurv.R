#' #' Produce Bivariate Alternating Recurrent Series Plot
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
#' @param object must be an object of \code{bivrecSurv} class.
#' @param main optional string with plot title (default is no title)
#' @param xlab optional string with label for horizontal axis (default is "Gap Times")
#' @param ylab optional string with label for vertical axis (default is "Individual")
#' @param type1 optional string to label type 1 gap times (default is "Type 1")
#' @param type2 optional string to label type 2 gap times (default is "Type 2")
#'
#' @examples
#' library(BivRec)
#'
#' #Simulate data
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#'
#' #Create a bivrecSurv object
#' bivrecSurvObject <- with(bivrec_data, bivrecSurv(id, epi, xij, yij, d1, d2))
#'
#' #Plot
#' plot(bivrecSurvObject)
#'
#' @export
#'
plot.bivrecSurv <- function(object, main, xlab, ylab, type1, type2){
  #check arguments for labels
  if (missing(main)) {main=""}
  if (missing(xlab)) {xlab="Gap Times"}
  if (missing(ylab)) {ylab="Individual"}
  if (missing(type1)) {type1="Type 1"}
  if (missing(type2)) {type2="Type 2"}

  args = c(main, xlab, ylab, type1, type2)

  ##EXTRACT VECTORS FOR PLOTTING FUNCTION
  parameters <- object$data4Creg[-(5:7)]
  colnames(parameters) <- c("id", "episode", "xij", "yij", "ci")
  ctimes <- object$data4Lreg$ctime
  nsubject <- object$data4Lreg$n
  temp <- NULL

  basicplot(parameters, ctimes, nsubject, temp, args, a = 1/2, b = 1/3)

}
