#' #' Produce Bivariate Alternating Recurrent Series Plot
#'
#' @description
#' This function plots bivariate recurrent event gap times by levels of the categorical covariates that were included in a semi-parametric regression using either Lee et al or Chang methods.
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#'
#' @param object must be an object of \code{bivrecReg} class.
#' @param main optional string with plot title (default is no title)
#' @param xlab optional string with label for horizontal axis (default is "Gap Times")
#' @param ylab optional string with label for vertical axis (default is "Individual")
#' @param type1 optional string to label type 1 gap times (default is "Type 1")
#' @param type2 optional string to label type 2 gap times (default is "Type 2")
#' @param data only for formula
#'
#' @examples
#' \dontrun{
#'# Simulate bivariate alternating recurrent event data
#' library(BivRec)
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#'
#' # Apply Lee C, Huang CY, Xu G, Luo X (2017) method
#' lee_reg <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
#'                     data = bivrec_data, method="Lee.et.al")
#'
#' #Plot bivrecReg object
#' plot(lee_reg)
#' }
#' @export
plot.bivrecReg <- function(object, main, xlab, ylab, type1, type2, data=NULL) {
    newobject = object$formula
    newdata = object$data$original
    plot.formula(newobject, newdata, main, xlab, ylab, type1, type2)
    }

