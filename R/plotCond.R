
#' Plot conditional CDF after non-parametric analysis
#'
#' @description
#' This function plots conditional cdf for recurrent event data.
#'
#' @param x must be an object of \code{bivrecNP} class where the analysis has specified conditional = TRUE.
#'
#' @return A plot of conditional cdf in the given interval.
#'
#' @importFrom stats ftable
#' @examples
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' bdat <- is.bivrecSurv(bivrec_data)
#' npresult <- bivrecNP(bdat,ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),conditional = TRUE, given.interval=c(0, 10))
#' plotCond(npresult)
#'
plotCond <- function(x) {
  if (!is.bivrecNP(x)) stop("Object must be a bivrecNP class")
  cond <-x$conditional.cdf$conditional
  plot(cond$Time, cond[,5], type="l", lty = 2, xlab = "Type II Gap Times (y)", ylab = "Conditional Probability",
  xlim=c(0, round(max(x$conditional.cdf$ygrid), digits=1)),
  ylim=c(0, round(max(cond[,5]), digits=1)),
  main=substitute(
  paste("P(", Y^0 <= y, "|", X^0 %in% "[", gi1, ",", gi2, "])"),
  list(gi1 = x$given.interval[1], gi2 = x$given.interval[2]))
  )
  graphics::lines(cond$Time, cond[,4], lty = 2)
  graphics::lines(cond$Time, cond$Conditional.Probability,lty = 1)
}



