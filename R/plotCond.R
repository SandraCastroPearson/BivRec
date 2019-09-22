
#' Plot conditional CDF after non-parametric analysis
#'
#' @description
#' This function plots conditional cdf for recurrent event data. For example see \code{bivrecNP}.
#'
#' @param x must be an object of \code{bivrecNP} class where the analysis has specified conditional = TRUE.
#'
#' @return A plot of conditional cdf in the given interval.
#'
#' @importFrom stats ftable
#' @keywords internal

plotCond <- function(x) {
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



